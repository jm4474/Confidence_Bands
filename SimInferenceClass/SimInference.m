classdef SimInference
    %This Matlab class provides some functions to implement simultaneous
    %confidence bands in different problems (be it time series or cross-sectional)
    
    properties
    end
    
    methods(Static=true)
    %% 1)Transformation of Interest:    
        function  theta=theta(mu,h)
            %This function computes the transformation of interest
            %(theta in the notation used by Montiel-Olea and Plagborg-Moller)
            %-----
            %INPUT
            %-----
            %a) mu in R^p x I
            %b) h is a function that gives a 1 x k vector
            %-----
            %OUTPUT
            %-----
            %a) theta is I x k
            
            %----------------------------------------------
            I=size(mu,2);
            k=size(h(mu(:,1)),2);
            theta=zeros(I,k);
            for ix=1:I
                theta(ix,:)=h(mu(:,ix));
            end                        
        end
        
    %% 2)Standard Errors based on resampling:
        function [thetadraws,rs_stderrors,Vartheta]=rs_stderrors(mudraws,h)
            %This function estimates the covariance matrix of 
            %thetahat based on a vector of mdraws. The name of the function
            %"rs_stderrors" is an abbreviation for "resampled standard errors"
            %(mudraws should be p times I)
            %------
            %INPUT
            %------
            %a)mudraws is p x I (I is the number of draws)
            %b)h is the function of interest
            %------
            %OUTPUT
            %a)rs_stderrors: estimate of the std. errors for each 
            %coordinate based on resampling (1 x k)
            %b) Vartheta: Covariance matrix of theta k x k 
            %------
            thetadraws=SimInference.theta(mudraws,h);
            rs_stderrors=std(thetadraws,0,1);
            Vartheta=cov(thetadraws');
        end
    %% 3)Typical Critical Values for Rectangular Bands:
        function [Pwise,Sidak,Bonferroni,thetaproj,muproj]...
            =critvalues(confidence,p,k)
        %This function reports typical critical values for 
        %simultaneous inference
        Pwise=norminv(1-((1-confidence)/2),0,1);
        Sidak=chi2inv(confidence^(1/k),1)^.5;  
        Bonferroni=norminv(1-((1-confidence)/(2*(k))),0,1);
        thetaproj=chi2inv(confidence,k)^.5;
        muproj=chi2inv(confidence,p)^.5;
%         if k > p
%             disp('Note that the thetaproj critical value is not valid if k>p.')
%             disp('(See Inoue and Kilian (2016) in JOE)')
%         else
%         end
        
        end
        %% 4)(Plug-in) Supt Critical Value for Rectangular Bands:
        function [suptplugin]=suptcritval_plugin(confidence,Vartheta,I)
            %This function reporst the supt critical value, based on MVN
            %draws.
            d=diag(Vartheta).^(.5);   % Extracts the diagonal of Vartheta
            d_nonzero = (d>eps); % Elements with nonzero variance
            Corrmat=bsxfun(@ldivide, d(d_nonzero), bsxfun(@rdivide, Vartheta(d_nonzero,d_nonzero), d(d_nonzero)')); %Matrix of Correlations
            [Corrmat_V,Corrmat_D] = eig(Corrmat);
            Corrmat_chol = sqrt(Corrmat_D)*Corrmat_V';
            
            t=abs(randn(I,size(Corrmat_chol,1))*Corrmat_chol); %limit of tstats
            suptplugin=quantile(max(t,[],2),confidence,1);                 
        end
        %% 5) 2-D plots of Rectangular Bands based on Theoretical C-values
        function [Rband]=Rbands(cval,stderror,thetahat)
            %INPUT
            %-----
            %a) cval is a scalar
            %b) stderror is a 1 x k vector
            %c) thetahat is the point estimator of theta (1 x k vector)
            %-----
            %OUTPUT
            %-----
            %Rband is a 2 times K matrix
            %-----
            
            Rband(1,:)=thetahat-(cval*stderror);
            Rband(2,:)=thetahat+(cval*stderror);
        end
        %% 6) "Credibility" of rectangular bands based on resampled values of mu
        function [cred]=credibility_Rbands(thetadraws,Rband)
            %INPUT
            %------
            %a)mudraws is p x I (I is the number of draws)
            %b)h is the function of interest
            %c)Rband is a 2 times K matrix (the first row is the lower bound)
            %------
            %OUTPUT
            %a)cred is a scalar reporting the number of draws
            %that fall into a particular confidence band
            k=size(thetadraws,2);
            aux1=sum(bsxfun(@le,thetadraws,Rband(2,:)).*...
                bsxfun(@ge,thetadraws,Rband(1,:)),2);
            cred=mean(aux1==k);
        end
        
        %% 7) Confidence Band based on the quantiles of mudraws
        function [quantileRband]=quantile_Rbands(thetadraws,confidence)
        quantileRband(1,:)=quantile(thetadraws,(1-confidence)/2,1);
        quantileRband(2,:)=quantile(thetadraws,1-((1-confidence)/2),1);
        end
        
        %% 8) "Calibrated" Confidence Band
        function [Calibrated,x,Cred]=calibrated_Rbands(thetadraws,confidence)
            k=size(thetadraws,2);
            f=@(x) SimInference.credibility_Rbands(thetadraws,SimInference.quantile_Rbands(thetadraws,max(min(x,1),0)))-confidence;
            x0=.5*(1-(1-confidence)/k) + .5*confidence;
            [x,fval]=fzero(f,x0);
            Calibrated=SimInference.quantile_Rbands(thetadraws,x);
            Cred=fval+confidence;
        end
        
        %% 9) Wrapper function for returning any number of desired plug-in bands
        function the_bands = bands_plugin(thetahat, Vartheta, p, band_list, I, confidence)
            % Returns a cell array of 2D bands
            % band_list is a cell array containing any combination of: 'Pwise', 'Sidak', 'supt', 'Bonferroni', 'thetaproj', 'muproj'
            numbands = length(band_list);
            the_bands = cell(1,numbands);
            k = size(Vartheta,1);
            stderror = sqrt(diag(Vartheta))';
            cv = struct;
            [cv.Pwise,cv.Sidak,cv.Bonferroni,cv.thetaproj,cv.muproj]=SimInference.critvalues(confidence,p,k);
            for j=1:numbands
                if strcmp(band_list{j}, 'supt')
                    the_cv = SimInference.suptcritval_plugin(confidence,Vartheta,I);
                else
                    the_cv = cv.(band_list{j});
                end
                the_bands{j} = SimInference.Rbands(the_cv,stderror,thetahat);
            end
        end
        
    end
         
end
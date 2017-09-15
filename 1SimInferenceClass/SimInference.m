classdef SimInference
% -----------------------------------------------------------------------------
% This Matlab class provides generic functions to implement simultaneous
% confidence bands in different problems (be it time series or cross-sectional)
%
% This version: September 15th, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ------------------------------------------------------------------------------
    
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
            %a)  mu:  Matrix collecting I different values of 
            %         the parameter mu                 (p x I)
            %b)   h:  matlab function that transforms mu 
            %         in R^p into the parameter theta. 
            %         (the output of h is a 1 x k vector)
            %-----
            %OUTPUT
            %-----
            %a)theta: parameter of interest                (I x k)
            %         Each row contains one evaluation of 
            %         the parameter of interest.            
            %----------------------------------------------
            I      = size(mu,2);
            
            k      = size(h(mu(:,1)),2);
            
            theta  = zeros(I,k);
            
            for ix = 1:I
                
                theta(ix,:)...
                   = h(mu(:,ix));
                
            end    
            
        end
        
    %% 2)Standard Errors based on resampling:
    
        function [thetadraws,rs_stderrors,Sigma]=rs_stderrors(mudraws,h)
            %This function estimates the covariance matrix of 
            %thetahat based on a vector of mdraws. The name of the function
            %"rs_stderrors" is an abbreviation for "resampled standard errors"
            %(mudraws should be p times I)
            %------
            %INPUT
            %------
            %a)mudraws: Matrix collecting I different values of 
            %           the parameter mu                   (p x I)
            %b)      h: matlab function that transforms mu 
            %           in R^p into the parameter theta. 
            %           (the output of h is a 1 x k vector)
            %------
            %OUTPUT
            %a)rs_stderrors: estimate of the std. errors for each 
            %                coordinate based on resampling (1 x k)
            %b) Sigma:       Covariance matrix of theta     (k x k) 
            %------
            
            thetadraws   = SimInference.theta(mudraws,h);
            
            rs_stderrors = std(thetadraws,0,1);
            
            Sigma        = cov(thetadraws');
            
        end
        
        
    %% 3)Critical Values for Rectangular Bands:
    
        function [Pwise,Sidak,Bonferroni,thetaproj,muproj]...
            =critvalues(confidence,p,k)
            %This function reports typical critical values for 
            %simultaneous inference (Pointwise, Sidak, Bonferroni,
            %thetaproj,muproj)
            %INPUT
            %------
            %a)confidence: confidence level                 (1 x 1)
            %b)         p: dimension of the parameter mu    (1 x 1)
            %c)         k: dimension of the parameter theta (1 x 1)   
            %------
            %OUTPUT
            %a)     Pwise: pointwise critical value         (1 x 1)
            %b)     Sidak: Sidak's critical value           (1 x 1)
            %c)Bonferroni: Bonferroni critical value        (1 x 1) 
            %d) thetaproj: Critical value corresponding to  (1 x 1)
            %              the projection of theta
            %e)    muproj: Critical value corresponding to  (1 x 1)
            %              the projection of mu
            %------
        Pwise      = norminv(1-((1-confidence)/2),0,1);
        
        Sidak      = chi2inv(confidence^(1/k),1)^.5;  
        
        Bonferroni = norminv(1-((1-confidence)/(2*(k))),0,1);
        
        thetaproj  = chi2inv(confidence,k)^.5;
        
        muproj     = chi2inv(confidence,p)^.5;
        
        end
        
        
        %% 4)(Plug-in) Sup-t Critical Value for Rectangular Bands:
        
        function [suptplugin]=suptcritval_plugin(confidence,Sigma,I)
            %This function reports the sup-t critical value, based on
            %multivariate normal draws.
            %INPUT
            %------
            %a)confidence: confidence level                 (1 x 1)
            %b)     Sigma: cov matrix of theta              (k x k)
            %              (this matrix can be singular)
            %c)         I: number of draws to compute the 
            %              supt critical value              (1 x 1)   
            %------
            %OUTPUT
            %a)suptplugin: supt critical value              (1 x 1)
            %              (we use plug-in, as we the input
            %              Sigma in this function is a consistent
            %              estimator of the true covariance matrix)
            %------
            d           = diag(Sigma).^(.5);     % Extracts the diagonal 
                                                 % of Vartheta
                                               
            d_nonzero   = (d>eps);               % Elements with nonzero 
                                                 % variance
                                                 % This line implicitly 
                                                 % drops out all the zero
                                                 % variance components of
                                                 % Signa
                                               
            Corrmat     = bsxfun(@ldivide, ...
                                 d(d_nonzero), ...
                                 bsxfun(@rdivide,...
                                 Sigma(d_nonzero,d_nonzero),...
                                 d(d_nonzero)')); %Matrix of Correlations
                             
            % Compute square root using eigendecomposition
            [Corrmat_V,Corrmat_D] ...
                         = eig(Corrmat);
            
            Corrmat_sqrt = sqrt(Corrmat_D)*Corrmat_V';
            
            t            = abs(randn(I,size(Corrmat_sqrt,1))*Corrmat_sqrt); %limit of tstats
            
            suptplugin   = quantile(max(t,[],2),confidence,1);
            
        end
        
        
        %% 5) 2-D plots of Rectangular Bands based on Theoretical C-values
        
        function [Rband]=Rbands(cval,stderror,thetahat)
            %This function simply constructs a rectangular confidence band
            %of the form thetahat +- cval*stderror.
            %INPUT
            %-----
            %a)     cval:  critical value                 (1 x 1)
            %b) stderror:  std error                      (1 x 1)
            %c) thetahat:  estimator of theta             (1 x k)
            %-----
            %OUTPUT
            %-----
            %      Rband:  is a 2 times K matrix
            %-----
            
            Rband(1,:) = thetahat-(cval*stderror);
            
            Rband(2,:) = thetahat+(cval*stderror);
            
        end
        
        %% 6) "Credibility" of rectangular bands based on resampled values of mu
        
        function [cred]=credibility_Rbands(thetadraws,Rband)
            %This function computes the "credibility" of a rectangular band
            %based on resampled values of mu (be it bootstrap draws or
            %posterior draws). Note that we are using the term "credibility"
            %(see code below) simply to refer to the share of draws of mu 
            %that fall inside a given rectangular band. 
            %INPUT
            %------
            %a) thetadraws: bootstrap or posterior draws of theta (k x I)
            %b)      Rband: rectangular band                      (2 x K)
            %               (the first row is the lower bound)
            %------
            %OUTPUT
            %a)       cred: is a scalar reporting the number of 
            %               draws that fall into a particular 
            %               confidence band.                      (1 x 1)
            %-----
            k    = size(thetadraws,2);
            
            aux1 = sum(bsxfun(@le,thetadraws,Rband(2,:)).*...
                bsxfun(@ge,thetadraws,Rband(1,:)),2);
            
            cred = mean(aux1==k);
        end
        
        
        %% 7) Equal-tailed Confidence Band based on the quantiles of thetadraws
        
        function [quantileRband]=quantile_Rbands(thetadraws,confidence)
            %This function computes an equal-tailed confidence band based
            %on the quantiles of thetadraws
            %INPUT
            %------
            %a)    thetadraws: bootstrap or posterior draws of theta (k x I)
            %b)    confidence: marginal confidence level             (1 x 1)
            %------
            %OUTPUT
            %a) quantileRband: equal-tailed confidence band         (2 x k)
            %               draws that fall into a particular 
            %-----
            
        quantileRband(1,:) = quantile(thetadraws,(1-confidence)/2,1);
        
        quantileRband(2,:) = quantile(thetadraws,1-((1-confidence)/2),1);
        end
        
        
        %% 8) "Calibrated" Confidence Band
        
        function [Calibrated,x,Cred]=calibrated_Rbands(thetadraws,confidence)
            %This function "caibrates" the marginal confidence/credibility
            %level of the confidence/credible band to achieve a desired 
            %simulatenous confidence/credibility level
            %INPUT
            %------
            %a)    thetadraws: bootstrap or posterior draws of theta (k x I)
            %b)    confidence: target simultaneous 
            %                  confidence/credibility level          (1 x 1)
            %------
            %OUTPUT
            %a)    Calibrated: Calibrated confidence/credible band   (2 x k)
            %b)             x: Calibrated marginal 
            %                  confidence/credibility level          (1 x 1)
            %c)          Cred: Calibrated simulatenous 
            %                  confidence/credibility level          (1 x 1)
            %-----
            k          = size(thetadraws,2);
            
            f          = @(x) ...
                         SimInference.credibility_Rbands(thetadraws,...
                         SimInference.quantile_Rbands(thetadraws,max(min(x,1),0)))...
                         -confidence;
                     
            x0         = .5*(1-(1-confidence)/k) + .5*confidence; 
                                          %Initial condition (marginal level)
              
            [x,fval]   = fzero(f,x0);     %Calibrate the marginal conf/cred     
            
            Calibrated = SimInference.quantile_Rbands(thetadraws,x);
                                          %Compute the simulatenous
                                          %conf/cred
                                          
            Cred       = fval+confidence; %If the numerical procedure worked,
                                          %fval should be close to zero
                                          %(this cred should be close to
                                          %the target conf/cred level)
                                         
        end
        
        %% 9) Wrapper function for returning desired plug-in bands
        
        function the_bands = bands_plugin(thetahat, Vartheta, p, band_list, I, confidence)
            % This functions returns a cell array of 2D bands
            %INPUT
            %------
            %a)    thetahat: point estimator of theta              (k x 1)
            %b)    Vartheta: estimated cov matrix of thetahat      (k x k)
            %c)           p: dimension of mu                       (1 x 1)
            %d)   band_list: is a cell array containing any 
            %                combination of: 'Pwise', 'Sidak', 
            %                'supt', 'Bonferroni', 'thetaproj', 
            %                'muproj'                              (1 x l)
            %e)           I: number of draws to compute suptplugin (1 x 1)
            %    confidence: nominal confidence level              (1 x 1)
            %------
            %OUTPUT
            %a)   the bands: cell array containing the cbands      (1 x l)
            %-----
            % 
            numbands  = length(band_list);      % numbands = l
            
            the_bands = cell(1,numbands);
            
            k         = size(Vartheta,1);
            
            stderror  = sqrt(diag(Vartheta))';
            
            cv        = struct;                % Initialize structure
            
            [cv.Pwise,...
             cv.Sidak,...
             cv.Bonferroni,...
             cv.thetaproj,...
             cv.muproj]...
                     = SimInference.critvalues(confidence,p,k);
                 
            for j    = 1:numbands
                
                if strcmp(band_list{j}, 'supt')
                    
                    the_cv ...
                     = SimInference.suptcritval_plugin(confidence,...
                                                       Vartheta,I);
                
                else
                    
                    the_cv ...
                     = cv.(band_list{j});
                end
                
                the_bands{j} = SimInference.Rbands(the_cv,...
                                                   stderror,...
                                                   thetahat);
                
            end
        end
        
    end
         
end
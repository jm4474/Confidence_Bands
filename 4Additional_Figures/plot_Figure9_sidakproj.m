clear;

% Figure comparing Sidak and mu-projection
% "Comparison of ?idák and ?-projection critical values"
% (Figure 9 in the paper)

% This script takes approximately .3 seconds to run 
% on a Macbook Pro 
%         @2.4 GHz Intel Core i7 (8 GB 1600 MHz DDR3)
% running Matlab R2016b.

% This version: September 15, 2017
% J. L. Montiel Olea & M. Plagborg-Moller

%% Settings

alphas          = [0.05 0.1 0.32];              % Significance levels to plot

plot_colors     = {[0 0 0]; [1 0 0]; [0 0 1]};  % Line colors for each signif.level 

plot_linestyles = {'-'; '--'; ':'};             % Line styles for each signif. level

ps              = 1:15;                         % Values of dimension p to plot

%% Compute "break-even" value of p

xs              = zeros(length(alphas),length(ps)); 
                                               % Will contain values x st. chi_{1,x} = chi_{p,1-alpha}

the_legend      = cell(length(alphas),1);      % Will contain figure legend

for a=1:length(alphas)                         % For every signif_level...
    
    for j=1:length(ps) % For every dimension p...
        
        % Find root of this function (normalized to be on the "scale" 1)
        f         = @(x) (chi2inv(1-alphas(a),ps(j)) - chi2inv(x,1))/ps(j);
        
        % Store root
        xs(a,j)   = fzero(f, [1-alphas(a) 1-eps]);
        
    end
    
    the_legend{a} = sprintf('%s%3.2f', '\alpha=', alphas(a)); % Add signif. level to legend
    
end

% Compute k s.t. (1-alpha)^(1/k) = x
ks                = bsxfun(@rdivide, log(1-alphas)', log(xs));

%% Plot

figure('Unit', 'normalize', 'Position', [0.2 0.2 0.6 0.6]);

% Plot break-even k against p on semilog scale
h = semilogy(ps, ks, 'LineWidth', 2);

% Linestyles, labels, legend
set(h, {'color'}, plot_colors);

set(h, {'linestyle'}, plot_linestyles);

title('$k$ such that $\chi_{(1-\alpha)^{1/k}}=\chi_{p,1-\alpha}$', 'Interpreter', 'Latex');

xlabel('$p$', 'Interpreter', 'Latex');

ylabel('$k$ (log scale)', 'Interpreter', 'Latex');

legend(the_legend, 'Location', 'SouthEast');

set(gca, 'FontSize', 14);

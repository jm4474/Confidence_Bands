function plot_compare(Thetas, Theta_bands, the_xlabel, the_ylabel, the_xticks, the_legend, varargin)
% ------------------------------------------------------------
% Labeled plot comparing different point estimates and bands
%
% Inputs:
% - Thetas:                 cell array of point estimates
% - Theta_bands:            cell array of band upper and lower bounds
% - the_xlabel:             x axis label string
% - the_ylabel:             y axis label string
% - the_xticks:             x axis ticks
% - the_legend:             legend object
% - varargin (optional):    cell array of line styles for bands
%
% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ------------------------------------------------------------


%% Preliminaries

hori = length(Thetas{1})-1;
colormap lines;


%% Store IRFs and band bounds

irfs = nan(hori+1,length(Theta_bands));
lowers = zeros(hori+1,length(Theta_bands));
uppers = lowers;
for m=1:length(Theta_bands) % For each band
    if m<=length(Thetas)
        irfs(:,m) = Thetas{m};
    end
    lowers(:,m) = Theta_bands{m}(1,:);
    uppers(:,m) = Theta_bands{m}(2,:);
end


%% Draw plot

lines = plot_IRF_band(irfs, lowers, uppers, varargin{:});

grid on;
xlim([min(the_xticks) max(the_xticks)]);
xlabel(the_xlabel);
ylabel(the_ylabel);
set(gca, 'XTick', the_xticks);
set(gca, 'FontSize', 14);
legend(lines, the_legend, 'Location', 'Best', 'FontSize', 14); % Apply legend to last subplot

end
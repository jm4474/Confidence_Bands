function lines = plot_IRF_band(IRFs, lower, upper, varargin)
% ------------------------------------------------------------
% Plot IRFs and confidence bands
%
% Inputs:
% - IRFs:       matrix with each IRF along columns
% - lower:      matrix with lower bounds of bands
% - upper:      matrix with upper bounds of bands
% - varargin    (optional): cell array of line styles for bands
% Outputs:
% - lines:      vector of IRF line objects (for use in legends, etc.)
%
% This version: August 24, 2017
% J. L. Montiel Olea & M. Plagborg-Moller
% ------------------------------------------------------------


%% Plots
cmap = colormap;
hori = size(IRFs,1)-1;

if isempty(varargin)
    linestyles = repmat({'-'}, 1, size(IRFs,2));
else
    linestyles = varargin{1};
end

lines = [];
hold on;
for i=1:size(IRFs,2)
   plot(0:hori, IRFs(:,i), 'Color', cmap(i,:), 'LineWidth', 3);
%    plot(0:hori, IRFs(:,i), '.-', 'Color', cmap(i,:), 'LineWidth', 3, 'MarkerSize', 20);
   plot(0:hori, lower(:,i), 'Color', cmap(i,:), 'LineStyle', linestyles{i});
   l = plot(0:hori, upper(:,i), 'Color', cmap(i,:), 'LineStyle', linestyles{i});
   lines = [lines l];
end
line(xlim, [0 0], 'Color', 'k');
hold off;

end
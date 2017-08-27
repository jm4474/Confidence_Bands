function plot_compare(beta_hats, bands_plot, control_labels, bands_legend, bands_linestyle, bands_linewidth)


%% Plot point estimates

k = length(beta_hats);
plot(1:k, beta_hats, 'xk', 'LineWidth', 2, 'MarkerSize', 10);
the_xlim = [0.5 k+0.5];                 % Horizontal axis limits
line(the_xlim, [0 0], 'Color', 'k');    % Draw horizontal axis
hold on;


%% Plot bands

colormap lines;
cmap = colormap;
cmap = [0 0 0; cmap]; % Insert black color first in colormap

lines = [];

for i=1:length(bands_plot) % For every band...
    
    % Plot upper boundary of band
    the_plot = plot(1:k, bands_plot{i}(1,:), 'Color', cmap(i,:), 'LineStyle', bands_linestyle{i}, 'LineWidth', bands_linewidth(i));
    lines = [lines the_plot]; % Save line for drawing legend below
    
    % Plot lower boundary
    plot(1:k, bands_plot{i}(2,:), 'Color', cmap(i,:), 'LineStyle', bands_linestyle{i}, 'LineWidth', bands_linewidth(i));
    
end

hold off;


%% Legend, labels, ticks

legend(lines, bands_legend, 'Location', 'NorthEast', 'FontSize', 16);
ylabel('log points');
set(gca,'XGrid','on');
set(gca, 'XTick', 1:k);
set(gca, 'XTickLabel', control_labels);
set(gca, 'FontSize', 16);
xlim(the_xlim); % Enforce horizontal axis limits

end
clear;

% Plot 2-D confidence band illustration
% This version: August 27, 2017
% J. L. Montiel Olea & M. Plagborg-Moller


%% Settings

rho             = 0.9;          % Correlation
var_theta2      = 0.25;         % Variance of thetahat_2
alpha           = 0.1;          % Significance level
supt_numdraws   = 1e6;          % Number of draws for sup-t crit val
axislim_theta1  = [-2.5 2.5];   % Axis limits, theta_1
axislim_theta2  = [-1.5 1.5];   % Axis limits, theta_2
ticks_theta1    = -2:2;         % Axis ticks, theta_1
ticks_theta2    = -1:0.5:1;     % Axis ticks, theta_2
rng(20170217);                  % Seed random number generator

% Styles of different bands
regions = {'Wald', 'Pointwise', 'Sup-t', 'Sidak', 'Bonferroni', '\theta-projection'}; % Labels
colors =  {'k',    [0 0.7 0],   'r',     'b',     [1 0.8 0],    [0.8 0 0.8]};         % Colors
lstyles = {'-',    '-',         '-',     '-',     '--',         '-'};                 % Linestyles


%% Compute critical values

critval = zeros(1,length(regions));
critval(2) = norminv((1 + 1-alpha)/2);                  % Pointwise
critval(3) = quantile(max(abs(mvnrnd(zeros(1,2), ...
                             [1 rho; rho 1], ...
                             supt_numdraws)), ...
                          [], 2), ...
                      1-alpha);                         % Sup-t
critval(4) = norminv((1+(1-alpha)^(1/2))/2);            % Sidak
critval(5) = norminv((1 + 1-alpha/2)/2);                % Bonferroni
critval(6) = sqrt(chi2inv(1-alpha,2));                  % Theta-proj


%% Compute Wald ellipse

theta1vals = linspace(-critval(6), critval(6), 500); % theta_1 values
ellipse = sqrt(var_theta2)*real(0.5*bsxfun(@plus, ...
                                           2*rho*theta1vals, ...
                                           [1; -1]*sqrt((2*rho*theta1vals).^2 - 4*(theta1vals.^2 - (1-rho^2)*critval(6)^2)) ...
                                           ) ...
                                ); % Lower and upper theta_2 values in ellipse

                            
%% Plot figure

figure('Units', 'normalize', 'Position', [0.2 0.2 0.6 0.6]);

figs = gobjects(1, length(regions)); % Will collect all figure objects for legend

% Ellipse
figs(1) = plot(theta1vals, ellipse(1,:), 'Color', 'k', 'LineWidth', 2);
hold on;
plot(theta1vals, ellipse(2,:), 'Color', 'k', 'LineWidth', 2);

% Rectangles
for j=2:length(regions)
    figs(j) = plot(NaN, NaN, 'Color', colors{j}, 'LineStyle', lstyles{j}, 'LineWidth', 2); % Empty figure for legend
    rectangle('Position', [-1 -sqrt(var_theta2) 2 2*sqrt(var_theta2)]*critval(j), ...
              'EdgeColor', colors{j}, ...
              'LineStyle', lstyles{j}, ...
              'LineWidth', 2); % Actual rectangle to plot
end
hold off;

% Display axes
line(axislim_theta1, [0 0], 'Color', 'k');
line([0 0], axislim_theta1, 'Color', 'k');

% Axes settings, labels, legend
axis equal;
xlim(axislim_theta1);
ylim(axislim_theta2);
xlabel('$\theta_1$', 'Interpreter', 'Latex');
ylabel('$\theta_2$', 'Interpreter', 'Latex');
set(gca, 'FontSize', 16);
set(gca, 'XTick', ticks_theta1);
set(gca, 'YTick', ticks_theta2);
legend(figs, regions, 'Location', 'EastOutside', 'FontSize', 16);

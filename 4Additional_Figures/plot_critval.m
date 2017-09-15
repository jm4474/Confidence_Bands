clear;

% Critical value comparison figure
% This version: August 27, 2017
% J. L. Montiel Olea & M. Plagborg-Moller


%% Settings

alphas      = [0.05 0.1 0.32];  % Significance levels
colors      = {'k', 'r', 'b'};  % Color for each signif. level
linestyles  = {'-', '--', ':'}; % Linestyle for each signif. level
ks          = 1:50;             % Values of dimension k
font_size   = 12;               % Figure font size
title_size  = 16;               % Figure title size


%% Compute critical values

% Pointwise
pw_critval = norminv(1-alphas'/2);                                  % (1-alpha) quantile of |N(0,1)|

% Sidak
sidak_critval = norminv((1+bsxfun(@power, 1-alphas', 1./ks))/2);    % (1-alpha)^(1/k) quantile of |N(0,1)|

% Bonferroni
bonf_critval = norminv(1-bsxfun(@rdivide, alphas', ks)/2);          % (1-alpha/k) quantile of |N(0,1)|

% theta-projection
proj_critval = zeros(length(alphas), length(ks));
for k=1:length(ks)
    proj_critval(:,k) = sqrt(chi2inv(1-alphas', ks(k)));            % (1-alpha) quantile of sqrt(chi^2(k))
end


%% Plot

figure('Unit', 'normalize', 'Position', [0.25 0.05 0.5 0.9]);

% Sidak vs. pointwise
subplot(3,1,1);
hold on;
for a=1:length(alphas)
    plot(ks, sidak_critval(a, :)/pw_critval(a), 'Color', colors{a}, 'LineStyle', linestyles{a}, 'LineWidth', 2);
end
hold off;
set(gca, 'FontSize', font_size);
title('Sidak vs. pointwise: $\chi_{(1-\alpha)^{1/k}}/\chi_{1-\alpha}$', 'Interpreter', 'Latex', 'FontSize', title_size);

% Bonferroni vs. Sidak
subplot(3,1,2);
hold on;
for a=1:length(alphas)
    plot(ks, bonf_critval(a,:)./sidak_critval(a, :), 'Color', colors{a}, 'LineStyle', linestyles{a}, 'LineWidth', 2);
end
hold off;
set(gca, 'FontSize', font_size);
title('Bonferroni vs. Sidak: $\chi_{1-\alpha/k}/\chi_{(1-\alpha)^{1/k}}$', 'Interpreter', 'Latex', 'FontSize', title_size);

% theta-projection vs. Sidak
subplot(3,1,3);
hold on;
the_legend = cell(1,length(alphas));
for a=1:length(alphas)
    plot(ks, proj_critval(a,:)./sidak_critval(a, :), 'Color', colors{a}, 'LineStyle', linestyles{a}, 'LineWidth', 2);
    the_legend{a} = sprintf('%s%3.2f', '\alpha = ', alphas(a));
end
hold off;
xlabel('$k$', 'Interpreter', 'Latex');
set(gca, 'FontSize', font_size);
title('$\theta$-projection vs. Sidak: $\chi_{k,1-\alpha}/\chi_{(1-\alpha)^{1/k}}$', 'Interpreter', 'Latex', 'FontSize', title_size);

% Add legend to last subplot
legend(the_legend, 'Location', 'SouthEast', 'FontSize', font_size);
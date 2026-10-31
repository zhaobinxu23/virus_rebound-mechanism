%% plot_Figure6.m
% Draw Figure 6 from the two saved fitting results.
% Panel d solves the cumulative probability equation
%   P(80 <= t_start <= 270) = 0.208
% for t_start ~ N(mu_t, sigma_t^2) truncated to t_start >= 0,
% with mu_t = 300.
%
% Required files:
%   parameter_estimation_rebound_rates.mat
%   estimate_tstart_sigma_p15fixed.mat

clear; clc; close all;

%% ---------- Load saved results ----------
file1 = 'parameter_estimation_rebound_rates.mat';
file2 = 'estimate_tstart_sigma_p15fixed.mat';

if ~isfile(file1)
    error('Cannot find %s. Run estimate_rebound_rates.m first.', file1);
end
if ~isfile(file2)
    error('Cannot find %s. Run estimate_tstart_sigma_p15_fixed.m first.', file2);
end

S1 = load(file1);
S2 = load(file2);

p15_grid          = S1.p15_grid(:);
rebound_untreated = logical(S1.rebound_untreated(:));
p15_crit          = S1.p15_crit;
mu_log            = S1.mu_log;
sigma_log         = S1.sigma_log;

t_start_grid  = S2.t_start_grid(:);
rebound_flags = logical(S2.rebound_flags(:));

%% ---------- Window, mean, target ----------
win_a  = 80;
win_b  = 270;
mu_t   = 300;
target = 0.208;

fprintf('p15_crit  = %.4f\n', p15_crit);
fprintf('sigma_log = %.4f\n', sigma_log);
fprintf('Window    = [%.0f, %.0f]\n', win_a, win_b);
fprintf('mu_t      = %.1f\n', mu_t);
fprintf('Target    = %.3f\n', target);

%% ---------- Solve P(win_a <= t_start <= win_b) = target ----------
% Cumulative probability as a function of sigma_t
cum_prob = @(s) truncated_window_mass(win_a, win_b, mu_t, s);

% Scan sigma on a log grid
sigma_grid = logspace(0, 3, 500);
mass_grid  = arrayfun(cum_prob, sigma_grid);

[max_mass, idx_max] = max(mass_grid);
fprintf('Max achievable cumulative P = %.4f at sigma = %.2f\n', ...
    max_mass, sigma_grid(idx_max));

if max_mass < target
    warning(['Max P = %.4f < target %.3f. No sigma_t exists for this ' ...
             'window and mean.'], max_mass, target);
    roots = [];
    sigma_primary = NaN;
else
    roots = [];
    for i = 1:length(sigma_grid)-1
        f1 = mass_grid(i)   - target;
        f2 = mass_grid(i+1) - target;
        if f1*f2 < 0
            r = fzero(@(s) cum_prob(s) - target, ...
                      [sigma_grid(i), sigma_grid(i+1)]);
            roots = [roots; r];
        end
    end
    if isempty(roots)
        sigma_primary = NaN;
    else
        roots = sort(roots);
        % Select smaller root as primary (more compact, biologically plausible)
        sigma_primary = roots(1);
        fprintf('Solutions: sigma_t = %s\n', mat2str(roots', 4));
        fprintf('Primary (smaller) sigma_t = %.2f, variance = %.2f\n', ...
            sigma_primary, sigma_primary^2);
    end
end

%% ---------- Plot ----------
set(groot,'DefaultAxesFontName','Arial');
set(groot,'DefaultAxesFontSize',9);
set(groot,'DefaultLineLineWidth',1.3);
set(groot,'DefaultTextFontName','Arial');

figure('Color','w','Units','centimeters','Position',[2 2 22 16]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% ---------- Panel a: untreated rebound vs p15 ----------
nexttile;
semilogx(p15_grid, double(rebound_untreated), 'ko-', ...
    'MarkerSize',5, 'MarkerFaceColor','k');
hold on;
xline(p15_crit, 'r--', sprintf('p15_{crit} = %.2f', p15_crit), ...
    'LineWidth',1.0);
xlabel('Immune competence p15');
ylabel('Untreated rebound (0/1)');
title('a  Untreated rebound vs p15','FontWeight','normal');
ylim([-0.1 1.1]); box off;

% ---------- Panel b: fitted lognormal distribution ----------
nexttile;
p15_fine = logspace(log10(min(p15_grid)), log10(max(p15_grid)), 400);
pdf_fine = exp(-0.5*((log10(p15_fine)-mu_log)./sigma_log).^2) ./ ...
           (sigma_log*sqrt(2*pi)*p15_fine*log(10));
plot(p15_fine, pdf_fine, 'k-','LineWidth',1.6); hold on;
xline(p15_crit, 'r--', sprintf('p15_{crit} = %.2f', p15_crit), ...
    'LineWidth',1.0);
xline(10^mu_log, 'b:', sprintf('median = %.1f', 10^mu_log), ...
    'LineWidth',1.0);
xlabel('p15'); ylabel('Probability density');
title(sprintf('b  Lognormal fit, \\sigma_{log} = %.3f', sigma_log), ...
    'FontWeight','normal');
box off;

% ---------- Panel c: treated rebound window ----------
nexttile;
stairs(t_start_grid, double(rebound_flags), 'k-','LineWidth',1.3); hold on;
yline(0.5, 'Color', [0.6 0.6 0.6]);
xline(win_a, 'r--', sprintf('%d', win_a), 'LineWidth',1.0);
xline(win_b, 'r--', sprintf('%d', win_b), 'LineWidth',1.0);
xlabel('Treatment start time (model time)');
ylabel('Treated rebound (0/1)');
title('c  Treated rebound window (p15 = 5)','FontWeight','normal');
ylim([-0.1 1.1]); box off;
xlim([t_start_grid(1) t_start_grid(end)]);

% ---------- Panel d: sigma_t fit ----------
nexttile;
semilogx(sigma_grid, mass_grid, 'k-','LineWidth',1.6); hold on;
yline(target, 'r--', sprintf('%.1f%% target', 100*target), ...
    'LineWidth',1.0);
if ~isempty(roots)
    for k = 1:length(roots)
        if k == 1
            label = sprintf('\\sigma_t = %.1f (primary)', roots(k));
        else
            label = sprintf('\\sigma_t = %.1f', roots(k));
        end
        xline(roots(k), 'b:', label, 'LineWidth',1.2);
    end
end
xlabel('\sigma_t (model time)');
ylabel('Cumulative P(85 \leq t_{start} \leq 270)');
title(sprintf('d  \\mu_t = %.0f, window [%d, %d]', mu_t, win_a, win_b), ...
    'FontWeight','normal');
box off; xlim([sigma_grid(1) sigma_grid(end)]);

%% ---------- Export ----------
exportgraphics(gcf, 'Figure6.png', 'Resolution', 300);
exportgraphics(gcf, 'Figure6.pdf', 'ContentType', 'vector');
savefig(gcf, 'Figure6.fig');

fprintf('Figure 6 saved: Figure6.png / Figure6.pdf / Figure6.fig\n');

%% ---------- Local function: cumulative probability ----------
function p = truncated_window_mass(a, b, mu, sigma)
    % P(a <= X <= b) for X ~ N(mu,sigma^2), truncated to X >= 0
    % = [Phi((b-mu)/sigma) - Phi((a-mu)/sigma)] / [1 - Phi(-mu/sigma)]
    if sigma <= 0
        p = 0; return;
    end
    Z = 1 - normcdf(0, mu, sigma);
    if Z < 1e-12
        p = 0; return;
    end
    p = (normcdf(b, mu, sigma) - normcdf(a, mu, sigma)) / Z;
    p = max(min(p, 1), 0);
end
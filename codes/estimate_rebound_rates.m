%% estimate_rebound_rates.m
% Two-stage probabilistic estimation using clinical rebound rates
% Targets: untreated rebound rate = 1.8%, treated rebound rate = 20.8%
% Reference: Edelstein et al., Ann Intern Med 2023
%
% Stage 1: p15 ~ lognormal, median fixed at 5; estimate sigma_log
%          from untreated rebound rate = 1.8%.
% Stage 2: fix p15 distribution; estimate kD (p.para(23))
%          from treated rebound rate = 20.8%.
% Stage 3: Monte Carlo verification.

clear; clc;

%% ---------- Load base model ----------
run('main_consider_drug_dose.m');   % defines p, y0, tspan, options

p_base = p;
p_base.para(15) = 5;      % baseline immune competence (median)
p_base.para(23) = 0.01;   % baseline kD, will be estimated

tspan = [0 10000];
options = odeset('RelTol',1e-4,'AbsTol',1e-8, ...
                 'NonNegative',1:numel(y0));

% Disable extinction switch during untreated scan to avoid artificial clearing
p_untreated = p_base;
p_untreated.use_extinction_switch = false;
p_untreated.extinction_V_threshold = 1e-20;
p_untreated.extinction_I_total_threshold = 1e-20;

%% ---------- Treated protocol ----------
% NOTE: model time unit is not yet mapped to days.
% Placeholder: treat_start=1, interval=0.5, 10 doses.
% After time mapping, replace these by the correct model-time values.
treat_start   = 100;
dose_interval = 100;
n_doses       = 10;
dose_amount   = 1e7;   % middle of scanned dose range

protocol_treated.start_time = treat_start;
protocol_treated.single_dose = dose_amount;
protocol_treated.num_doses = n_doses;
protocol_treated.dose_interval = dose_interval;
protocol_treated.dose_times = treat_start + (0:n_doses-1)*dose_interval;
protocol_treated.dose_amounts = dose_amount * ones(1,n_doses);
protocol_treated.last_dose_time = protocol_treated.dose_times(end);
protocol_treated.treatment_duration = (n_doses-1)*dose_interval;

%% ============================================================
% Stage 1: untreated rebound rate vs p15
% ============================================================
% p15_grid = logspace(log10(0.2), log10(20), 25);
p15_grid = logspace(log10(2), log10(5), 10);
rebound_untreated = false(size(p15_grid));

fprintf('Stage 1: scanning p15 for untreated rebound...\n');
for i = 1:length(p15_grid)
    p_try = p_untreated;
    p_try.para(15) = p15_grid(i);
    p_try.drug_dose_times = [];
    p_try.drug_dose_amounts = [];

    odeFun = @(t,y) sys_ode_modified_new(t,y,p_try);
    try
        [t, y] = ode15s(odeFun, tspan, y0, options);
        N = p_try.N;
        V = max(y(:,N+2),0);
        rebound_untreated(i) = detect_rebound_untreated(t, V);
    catch ME
        warning('p15=%g failed: %s', p15_grid(i), ME.message);
        rebound_untreated(i) = false;
    end
    fprintf('  p15=%.3g: rebound=%d\n', p15_grid(i), rebound_untreated(i));
end

%% Stage 1b: fit sigma_log of p15 distribution
p15_med = 5;
target_untreated = 0.018;

idx_true = find(rebound_untreated);
if isempty(idx_true)
    error(['No untreated rebound in p15 range. ' ...
           'Cannot fit 1.8%%. Check rebound definition or model regime.']);
elseif all(rebound_untreated)
    error('All p15 values rebound. Cannot fit 1.8%%.');
else
    p15_crit = max(p15_grid(rebound_untreated));
    fprintf('Critical p15 for untreated rebound: %.3g\n', p15_crit);
end

mu_log = log10(p15_med);
z_target = norminv(1 - target_untreated);   % ≈ 2.096
% If no Statistics Toolbox: z_target = sqrt(2)*erfinv(2*(1-target)-1);
sigma_log = (mu_log - log10(p15_crit)) / z_target;
fprintf('Fitted sigma_log(p15) = %.4f\n', sigma_log);
fprintf('Check P(p15<p15_crit) = %.4f (target %.4f)\n', ...
    normcdf((log10(p15_crit)-mu_log)/sigma_log), target_untreated);

% %% ============================================================
% % Stage 2: treated rebound rate vs kD
% % ============================================================
% kD_grid = logspace(log10(0.002), log10(0.05), 15);
% target_treated = 0.208;
% 
% % Numerical integration over p15 distribution using quantile midpoints
% n_quad = 15;
% q = ((1:n_quad) - 0.5) / n_quad;          % midpoint quantiles
% p15_quad = 10.^(mu_log + sigma_log*norminv(q));
% weights = ones(n_quad,1)/n_quad;
% 
% treated_rate = zeros(size(kD_grid));
% fprintf('Stage 2: scanning kD for treated rebound...\n');
% 
% for ik = 1:length(kD_grid)
%     kD = kD_grid(ik);
%     rebound_prob = 0;
%     for iq = 1:n_quad
%         p_try = p_base;
%         p_try.para(15) = p15_quad(iq);
%         p_try.para(23) = kD;
%         p_try.drug_dose_times = protocol_treated.dose_times;
%         p_try.drug_dose_amounts = protocol_treated.dose_amounts;
% 
%         sim = simulate_one_drug_protocol(p_try, y0, tspan, options, ...
%                                          protocol_treated);
%         if sim.success && sim.metrics.rebound_flag
%             rebound_prob = rebound_prob + weights(iq);
%         end
%     end
%     treated_rate(ik) = rebound_prob;
%     fprintf('  kD=%.4g: treated rebound rate=%.3f\n', kD, treated_rate(ik));
% end
% 
% %% Stage 2b: interpolate kD to match 20.8%
% if all(treated_rate < target_treated) || all(treated_rate > target_treated)
%     warning('Target %.3f outside scanned range [%.3f, %.3f]. Extend kD grid.', ...
%         target_treated, min(treated_rate), max(treated_rate));
%     kD_fit = NaN;
% else
%     kD_fit = interp1(treated_rate, kD_grid, target_treated, 'linear');
%     fprintf('Fitted kD = %.4g (model time units)\n', kD_fit);
%     fprintf('Corresponding half-life = %.4g model time units\n', log(2)/kD_fit);
% end
% 
% %% ============================================================
% % Stage 3: Monte Carlo verification
% % ============================================================
% n_mc = 300;
% rng(1);   % reproducibility
% p15_samples = 10.^(mu_log + sigma_log*randn(n_mc,1));
% rebound_mc = false(n_mc,1);
% 
% fprintf('Stage 3: Monte Carlo with %d samples...\n', n_mc);
% parfor i = 1:n_mc
%     p_try = p_base;
%     p_try.para(15) = p15_samples(i);
%     p_try.para(23) = kD_fit;
%     p_try.drug_dose_times = protocol_treated.dose_times;
%     p_try.drug_dose_amounts = protocol_treated.dose_amounts;
% 
%     sim = simulate_one_drug_protocol(p_try, y0, tspan, options, ...
%                                      protocol_treated);
%     rebound_mc(i) = sim.success && sim.metrics.rebound_flag;
% end
% fprintf('MC treated rebound rate = %.3f (target %.3f)\n', ...
%     mean(rebound_mc), target_treated);
% 
% %% ---------- Save ----------
save('parameter_estimation_rebound_rates.mat', ...
    'p15_grid','rebound_untreated','p15_crit','mu_log','sigma_log');

%% ============================================================
% Local function: untreated rebound detection
% ============================================================
function flag = detect_rebound_untreated(t, V)
    % Second peak after first nadir:
    % peak2 >= 10 * max(nadir,1) AND peak2 >= 100.
    t = t(:); V = max(V(:),0);
    [~, idx_peak1] = max(V);
    if idx_peak1 >= length(V)-2
        flag = false; return;
    end
    V_after = V(idx_peak1:end);
    [min_after, idx_min_rel] = min(V_after);
    if idx_min_rel >= length(V_after)-2
        flag = false; return;
    end
    V_post = V_after(idx_min_rel:end);
    peak2 = max(V_post);
    flag = (peak2 >= 10*max(min_after,1)) && (peak2 >= 100);
end
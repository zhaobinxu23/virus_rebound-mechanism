%% estimate_tstart_sigma_p15_fixed.m
% p15 fixed at 5
% Scan t_start, simulate rebound (deterministic 0/1 for fixed p15)
% Given mu_t, solve sigma_t such that P(t_start in rebound window) = 0.208
%
% Time scale: 200 time units = 1 day
% User-specified mean treatment start: mu_t = 400 (2 days)

clear; clc;
run('main_consider_drug_dose.m');

p_base = p;
p_base.para(15) = 5;      % FIXED
p_base.para(23) = 0.01;   % fixed kD

tspan = [0 4000];
options = odeset('RelTol',1e-4,'AbsTol',1e-8,'NonNegative',1:numel(y0));

% Protocol (fixed)
dose_amount   = 1e7;
dose_interval = 100;   % 12 h
n_doses       = 10;    % 5 days BID

% Scan range: extended to cover mu_t = 400
t_start_grid = 0:5:500;

% Simulate
rebound_flags = false(size(t_start_grid));
fprintf('Scanning t_start 0:5:%d (p15 = 5)...\n', t_start_grid(end));
for i = 1:length(t_start_grid)
    protocol = build_protocol(t_start_grid(i), dose_amount, dose_interval, n_doses);
    p_try = p_base;
    p_try.drug_dose_times = protocol.dose_times;
    p_try.drug_dose_amounts = protocol.dose_amounts;
    sim = simulate_one_drug_protocol(p_try, y0, tspan, options, protocol);
    rebound_flags(i) = sim.success && sim.metrics.rebound_flag;
    if mod(i,10)==0 || i==1
        fprintf('  t_start=%5.0f: rebound=%d\n', t_start_grid(i), rebound_flags(i));
    end
end

% Identify rebound windows
windows = find_rebound_windows(t_start_grid, rebound_flags);
fprintf('\nRebound windows:\n');
if isempty(windows)
    error('No rebound windows found. Cannot fit sigma_t.');
end
for k = 1:size(windows,1)
    fprintf('  [%.0f, %.0f], width = %.0f\n', ...
        windows(k,1), windows(k,2), windows(k,2)-windows(k,1));
end

% Given mu_t, solve sigma_t
mu_t   = 300;    % user-specified mean treatment start
target = 0.208;
dt     = t_start_grid(2) - t_start_grid(1);

% Grid of sigma
sigma_grid = logspace(0, 3.5, 300);
mass_grid  = arrayfun(@(s) window_mass(windows, mu_t, s, dt), sigma_grid);
[max_mass, idx_max] = max(mass_grid);
fprintf('\nMax achievable P(rebound) with mu_t=%.0f: %.4f (at sigma=%.1f)\n', ...
    mu_t, max_mass, sigma_grid(idx_max));

if max_mass < target
    warning(['Max achievable P(rebound) = %.4f < target %.4f for mu_t = %.0f. ' ...
             'No sigma_t exists. Options: (a) set mu_t inside or near the rebound window ' ...
             '(e.g. mu_t = 150-200), (b) widen the window by changing dose/protocol, ' ...
             '(c) increase dose.'], max_mass, target, mu_t);
    sigma_fit = NaN;
    crossings = [];
else
    % Find all crossings of mass_grid through target
    crossings = [];
    for i = 1:length(sigma_grid)-1
        f1 = mass_grid(i)   - target;
        f2 = mass_grid(i+1) - target;
        if f1*f2 < 0
            sc = fzero(@(s) window_mass(windows, mu_t, s, dt) - target, ...
                       [sigma_grid(i), sigma_grid(i+1)]);
            crossings = [crossings; sc];
        end
    end
    if isempty(crossings)
        sigma_fit = NaN;
    else
        % If multiple crossings, pick the largest sigma (widest, most conservative)
        sigma_fit = max(crossings);
        fprintf('Crossings at sigma = %s\n', mat2str(crossings', 4));
        fprintf('Selected sigma_t = %.2f, variance = %.2f\n', ...
            sigma_fit, sigma_fit^2);
    end
end

% Plot
figure('Color','w','Units','centimeters','Position',[2 2 22 8]);

subplot(1,2,1);
stairs(t_start_grid, double(rebound_flags), 'k-','LineWidth',1.2); hold on;
xlabel('t_{start} (model time)');
ylabel('Rebound indicator');
title('Rebound window (p15 = 5)');
ylim([-0.1 1.1]); box off;

subplot(1,2,2);
semilogx(sigma_grid, mass_grid, 'k-','LineWidth',1.2); hold on;
yline(target, 'r--','20.8% target');
if ~isnan(sigma_fit)
    xline(sigma_fit, 'b:','Fitted \sigma_t');
end
xlabel('\sigma_t');
ylabel('P(rebound)');
title(sprintf('\\mu_t = %.0f', mu_t));
box off;

% Save
save('estimate_tstart_sigma_p15fixed.mat', ...
    't_start_grid','rebound_flags','windows','mu_t','sigma_fit', ...
    'target','sigma_grid','mass_grid','crossings');

%% ============================================================
% Local functions
% ============================================================
function protocol = build_protocol(t_start, dose_amount, dose_interval, n_doses)
    protocol.start_time = t_start;
    protocol.single_dose = dose_amount;
    protocol.num_doses = n_doses;
    protocol.dose_interval = dose_interval;
    protocol.dose_times = t_start + (0:n_doses-1)*dose_interval;
    protocol.dose_amounts = dose_amount * ones(1,n_doses);
    protocol.last_dose_time = protocol.dose_times(end);
    protocol.treatment_duration = (n_doses-1)*dose_interval;
end

function windows = find_rebound_windows(t_grid, flags)
    windows = [];
    in_window = false;
    start_idx = 0;
    for i = 1:length(flags)
        if flags(i) && ~in_window
            start_idx = i;
            in_window = true;
        elseif ~flags(i) && in_window
            windows = [windows; t_grid(start_idx), t_grid(i-1)];
            in_window = false;
        end
    end
    if in_window
        windows = [windows; t_grid(start_idx), t_grid(end)];
    end
end

function mass = window_mass(windows, mu, sigma, dt)
    % P(N(mu, sigma^2) truncated to [0, inf) falls in union of windows)
    Z = 1 - normcdf(0, mu, sigma);
    if Z < 1e-12
        mass = 0; return;
    end
    mass = 0;
    for k = 1:size(windows,1)
        a = windows(k,1) - dt/2;
        b = windows(k,2) + dt/2;
        mass = mass + normcdf(b, mu, sigma) - normcdf(a, mu, sigma);
    end
    mass = mass / Z;
end


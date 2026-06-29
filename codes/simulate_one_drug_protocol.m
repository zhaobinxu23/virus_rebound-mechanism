function sim = simulate_one_drug_protocol(p_base, y0, tspan, ode_options, protocol)
% simulate_one_drug_protocol
% Run one drug protocol simulation and calculate rebound metrics.

%% =========================
% 1. Update drug protocol
% ==========================
p = p_base;

p.drug_dose_times = protocol.dose_times;
p.drug_dose_amounts = protocol.dose_amounts;

% Keep compatibility with older code
p.drug_on_time = protocol.dose_times(1);
p.drug_off_time = protocol.dose_times(end);

if ~isfield(p, 'drug_internal_dt')
    p.drug_internal_dt = 1;
end

if ~isfield(p, 'vin_cap')
    p.vin_cap = 1e50;
end

if ~isfield(p, 'vin_exp_cap')
    p.vin_exp_cap = 700;
end

%% =========================
% 2. Run ODE
% ==========================
odeFun = @(t, y) sys_ode_modified_new(t, y, p);

try
   [t, y, extinctionInfo] = solve_with_extinction_restart( ...
    odeFun, tspan, y0, ode_options, p);
catch ME
    warning('Simulation failed: start=%g, dose=%g, n=%d\n%s', ...
        protocol.start_time, protocol.single_dose, ...
        protocol.num_doses, ME.message);

    sim.t = [];
    sim.y = [];
    sim.p = p;
    sim.protocol = protocol;
    sim.metrics = failed_metrics();
    sim.success = false;
    sim.error_message = ME.message;
    return;
end

%% =========================
% 3. Extract variables
% ==========================
N = p.N;

T = y(:, 1);
I_total = sum(y(:, 2:N+1), 2);
V = y(:, N+2);
Tc = y(:, N+3);

T = max(T, 0);
I_total = max(I_total, 0);
V = max(V, 0);
Tc = max(Tc, 0);

%% =========================
% 4. Calculate rebound metrics
% ==========================
metrics = calculate_rebound_metrics(t, V, protocol);

%% =========================
% 5. Output
% ==========================
sim.t = t;
sim.y = y;
sim.T = T;
sim.I_total = I_total;
sim.V = V;
sim.Tc = Tc;
sim.p = p;
sim.protocol = protocol;
sim.metrics = metrics;
sim.success = true;
sim.error_message = '';
sim.extinctionInfo = extinctionInfo;
end

%% ========================================================================
% Local helper functions
% ========================================================================

function metrics = calculate_rebound_metrics(t, V, protocol)

% Default scalar outputs. Never return [].
metrics.peak_V_all = NaN;
metrics.time_peak_V_all = NaN;
metrics.min_V_after_treatment = NaN;
metrics.time_min_V_after_treatment = NaN;
metrics.rebound_peak_V = NaN;
metrics.time_rebound_peak = NaN;
metrics.rebound_amplitude = NaN;
metrics.rebound_flag = false;
metrics.auc_V = NaN;
metrics.auc_log10V = NaN;
metrics.final_V = NaN;

t = t(:);
V = V(:);
V = max(real(V), 0);

validIdx = isfinite(t) & isfinite(V);

if ~any(validIdx)
    return;
end

t = t(validIdx);
V = V(validIdx);

if isempty(t) || isempty(V)
    return;
end

metrics.final_V = V(end);

% Overall peak
[metrics.peak_V_all, idxPeakAll] = max(V);
metrics.time_peak_V_all = t(idxPeakAll);

% AUC
if numel(t) >= 2
    metrics.auc_V = trapz(t, V);
    metrics.auc_log10V = trapz(t, log10(V + 1));
else
    metrics.auc_V = 0;
    metrics.auc_log10V = 0;
end

% Protocol safety
if isfield(protocol, 'start_time')
    startTime = protocol.start_time;
else
    startTime = t(1);
end

if isfield(protocol, 'last_dose_time')
    lastDose = protocol.last_dose_time;
else
    lastDose = startTime;
end

if isfield(protocol, 'dose_interval')
    doseInterval = protocol.dose_interval;
else
    doseInterval = 0;
end

% Minimum viral load after treatment begins
idxAfterTreatment = t >= startTime;

if any(idxAfterTreatment)
    V_after = V(idxAfterTreatment);
    t_after = t(idxAfterTreatment);

    [metrics.min_V_after_treatment, idxMin] = min(V_after);
    metrics.time_min_V_after_treatment = t_after(idxMin);
else
    [metrics.min_V_after_treatment, idxMin] = min(V);
    metrics.time_min_V_after_treatment = t(idxMin);
end

% Rebound observation window
% Prefer one dose interval after last dose.
reboundStart = lastDose + doseInterval;
idxPost = t >= reboundStart;

% If no time point exists after reboundStart, fall back to after last dose.
if ~any(idxPost)
    idxPost = t >= lastDose;
end

% If still empty, fall back to after treatment start.
if ~any(idxPost)
    idxPost = t >= startTime;
end

% If still empty, use the full trajectory.
if ~any(idxPost)
    idxPost = true(size(t));
end

V_post = V(idxPost);
t_post = t(idxPost);

if isempty(V_post)
    return;
end

[metrics.rebound_peak_V, idxRebound] = max(V_post);
metrics.time_rebound_peak = t_post(idxRebound);

metrics.rebound_amplitude = metrics.rebound_peak_V - metrics.min_V_after_treatment;
metrics.rebound_amplitude = max(metrics.rebound_amplitude, 0);

% Rebound criterion
% 条件 1：后期峰值至少比治疗后最低值高 10 倍
% 条件 2：后期峰值本身超过 100
baselineForFold = max(metrics.min_V_after_treatment, 1);
foldCriterion = metrics.rebound_peak_V >= 10 * baselineForFold;
absoluteCriterion = metrics.rebound_peak_V >= 100;

metrics.rebound_flag = foldCriterion && absoluteCriterion;

end

function [tAll, yAll, extinctionInfo] = solve_with_extinction_restart( ...
    odeFun, tspan, y0, ode_options, p)

extinctionInfo.triggered = false;
extinctionInfo.time = NaN;

tStart = tspan(1);
tEnd = tspan(end);

eventOptions = odeset(ode_options, ...
    'Events', @(t, y) extinction_event(t, y, p));

[t1, y1, te, ye] = ode15s(odeFun, [tStart tEnd], y0, eventOptions);

if isempty(te)
    tAll = t1;
    yAll = y1;
    return;
end

eventTime = te(end);
eventState = ye(end, :)';

extinctionInfo.triggered = true;
extinctionInfo.time = eventTime;

% After extinction, infected cells, free virus, and virus-antibody
% complexes must be removed before restarting the ODE.
yRestart = enforce_extinction_state(eventState, p);
yRestart = max(yRestart, 0);

restartOptions = odeset(ode_options, 'Events', []);

[t2, y2] = ode15s(odeFun, [eventTime tEnd], yRestart, restartOptions);

tAll = [t1; t2(2:end)];
yAll = [y1; y2(2:end, :)];

% Enforce zero virus-related states for all post-extinction output rows.
idxPost = tAll >= eventTime;
yAll(idxPost, :) = enforce_extinction_state_matrix(yAll(idxPost, :), p);
yAll = max(yAll, 0);

end

function state = enforce_extinction_state(state, p)

N = p.N;

% Infected-cell age compartments.
state(2:N+1) = 0;

% Free extracellular virus.
state(N+2) = 0;

% x0 starts after T, I_1...I_N, V and Tc.
% Full state: y = [T; I_1...I_N; V; Tc; x0]
x0Offset = N + 3;

% x0(801:1200):
% 801:900    IgM-BCR-V complex
% 901:1000   IgM-V complex
% 1001:1100  IgG-BCR-V complex
% 1101:1200  IgG-V complex
virusComplexIdx = x0Offset + (801:1200);
virusComplexIdx = virusComplexIdx(virusComplexIdx <= numel(state));

state(virusComplexIdx) = 0;

end

function states = enforce_extinction_state_matrix(states, p)

if isempty(states)
    return;
end

N = p.N;

% Infected-cell age compartments.
states(:, 2:N+1) = 0;

% Free extracellular virus.
states(:, N+2) = 0;

% x0 starts after T, I_1...I_N, V and Tc.
% Full state: y = [T; I_1...I_N; V; Tc; x0]
x0Offset = N + 3;

% x0(801:1200):
% 801:900    IgM-BCR-V complex
% 901:1000   IgM-V complex
% 1001:1100  IgG-BCR-V complex
% 1101:1200  IgG-V complex
virusComplexIdx = x0Offset + (801:1200);
virusComplexIdx = virusComplexIdx(virusComplexIdx <= size(states, 2));

states(:, virusComplexIdx) = 0;

end

function [value, isterminal, direction] = extinction_event(~, y, p)

N = p.N;

I_total = sum(max(y(2:N+1), 0));
V = max(y(N+2), 0);

V_threshold = p.extinction_V_threshold;
I_threshold = p.extinction_I_total_threshold;

value = max(I_total / I_threshold, V / V_threshold) - 1;
isterminal = 1;
direction = -1;

end
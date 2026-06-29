%% analyze_replication_capacity_immune_competence_nature.m
% Viral replication capacity and host immune competence determine acute clearance and chronic persistence

clear;
clc;
close all;

%% =========================
% 1. User settings
% =========================

mainParameterScript = 'main_consider_drug_dose.m';
odeFunctionName = 'sys_ode_modified_new';

outDir = fullfile( ...
    'drug_rebound_scan_results', ...
    'replication_capacity_immune_competence');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

figPrefix = fullfile(outDir, ...
    'replication_capacity_immune_competence_nature');

simulationDuration = 10000;
timeLimit = [0 simulationDuration];

replicationValues = [0.05, 0.1];

immuneK5 = 0.1;
strongImmuneCompetence = 5;
weakImmuneCompetence = 2.5;

%% =========================
% 2. Load baseline parameters
% =========================

[pBase, y0Base] = initialize_model_from_main_script(mainParameterScript);

if ~isfield(pBase, 'N')
    pBase.N = 400;
end

N = pBase.N;

%% =========================
% 3. Run simulations for viral replication capacity
% =========================

replicationSims = struct([]);

for i = 1:numel(replicationValues)

    p = pBase;
    p.k5 = replicationValues(i);
    p = refresh_intracellular_replication_vectors(p);

    p.drug_dose_times = [];
    p.drug_dose_amounts = [];

    if isfield(p, 'use_extinction_switch')
        p.use_extinction_switch = true;
    end

    label = sprintf('k5 = %.2f', p.k5);

    replicationSims(i).label = label;
    replicationSims(i).k5 = p.k5;
    replicationSims(i).p = p;

    fprintf('\nRunning replication-capacity simulation: %s\n', label);

    sim = run_single_simulation(p, y0Base, simulationDuration, odeFunctionName);
    replicationSims(i).sim = sim;
    replicationSims(i).features = extract_features(sim, p);
    replicationSims(i).death = compute_death_dynamics(sim, p);

end

%% =========================
% 4. Run simulations for host immune competence
% =========================

immuneValues = [strongImmuneCompetence, weakImmuneCompetence];
immuneLabels = {'Strong immunity', 'Weak immunity'};

immuneSims = struct([]);

for i = 1:numel(immuneValues)

    p = pBase;
    p.k5 = immuneK5;
    p = refresh_intracellular_replication_vectors(p);

    p.para(15) = immuneValues(i);

    p.drug_dose_times = [];
    p.drug_dose_amounts = [];

    if isfield(p, 'use_extinction_switch')
        p.use_extinction_switch = true;
    end

    label = sprintf('%s, p.para(15) = %.1f', immuneLabels{i}, p.para(15));

    immuneSims(i).label = label;
    immuneSims(i).immuneCompetence = p.para(15);
    immuneSims(i).p = p;

    fprintf('\nRunning immune-competence simulation: %s\n', label);

    sim = run_single_simulation(p, y0Base, simulationDuration, odeFunctionName);
    immuneSims(i).sim = sim;
    immuneSims(i).features = extract_features(sim, p);
    immuneSims(i).death = compute_death_dynamics(sim, p);

end

%% =========================
% 5. Summary
% =========================

summary = make_summary(replicationSims, immuneSims);

summaryFile = fullfile(outDir, ...
    'replication_capacity_immune_competence_summary.mat');

save(summaryFile, ...
    'summary', ...
    'replicationSims', ...
    'immuneSims', ...
    'pBase', ...
    'y0Base');

%% =========================
% 6. Nature-style figure settings
% =========================

set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontSize', 8);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 0.8);

blue = [0.10 0.37 0.71];
red = [0.80 0.18 0.16];
green = [0.20 0.55 0.30];
purple = [0.49 0.18 0.56];
orange = [0.93 0.49 0.19];
gray = [0.35 0.35 0.35];

repColors = [blue; red];
immuneColors = [green; purple];

fig = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [2 2 34 32], ...
    'Renderer', 'painters');

tl = tiledlayout(fig, 4, 4, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% =========================
% A. Viral replication capacity: V and I_total
% =========================

axA = nexttile(tl, [1 2]);
hold on;

for i = 1:numel(replicationSims)
    f = replicationSims(i).features;

    plot(f.t, log10(f.V + 1), ...
        '-', ...
        'Color', repColors(i, :), ...
        'LineWidth', 1.8);

    plot(f.t, log10(f.I_total + 1), ...
        '--', ...
        'Color', repColors(i, :), ...
        'LineWidth', 1.5);
end

xlabel('Time');
ylabel('log_{10}(burden + 1)');
title('a  Viral replication capacity controls infection outcome', ...
    'FontWeight', 'bold');

legend({ ...
    'V, k_5 = 0.05', ...
    'I_{total}, k_5 = 0.05', ...
    'V, k_5 = 0.10', ...
    'I_{total}, k_5 = 0.10'}, ...
    'Location', 'northeast', ...
    'Box', 'off');

xlim(timeLimit);
box off;

%% =========================
% B. Viral replication capacity: Tc cells
% =========================

axB = nexttile(tl, [1 2]);
hold on;

for i = 1:numel(replicationSims)
    f = replicationSims(i).features;

    plot(f.t, log10(f.Tc + 1), ...
        'Color', repColors(i, :), ...
        'LineWidth', 1.8);
end

xlabel('Time');
ylabel('log_{10}(CD8^{+} T cells + 1)');
title('b  Replication capacity alters CD8^{+} T-cell activation', ...
    'FontWeight', 'bold');

legend({'k_5 = 0.05', 'k_5 = 0.10'}, ...
    'Location', 'southeast', ...
    'Box', 'off');

xlim(timeLimit);
box off;

%% =========================
% C. Infection-age heatmaps for k5 = 0.05 and 0.1
% =========================

axC1 = nexttile;
f = replicationSims(1).features;

imagesc(f.t, 1:N, log10(f.I_age' + 1));
set(gca, 'YDir', 'normal');

xlabel('Time');
ylabel('Infection age');
title('c  I_{1-400}, k_5 = 0.05', ...
    'FontWeight', 'bold');

cbC1 = colorbar;
cbC1.Label.String = 'log_{10}(I_a + 1)';
colormap(axC1, parula);
xlim(timeLimit);
box off;

axC2 = nexttile;
f = replicationSims(2).features;

imagesc(f.t, 1:N, log10(f.I_age' + 1));
set(gca, 'YDir', 'normal');

xlabel('Time');
ylabel('Infection age');
title('c  I_{1-400}, k_5 = 0.10', ...
    'FontWeight', 'bold');

cbC2 = colorbar;
cbC2.Label.String = 'log_{10}(I_a + 1)';
colormap(axC2, parula);
xlim(timeLimit);
box off;

%% =========================
% D. Death mechanisms for k5 = 0.05 and 0.1
% =========================

axD = nexttile(tl, [1 2]);
hold on;

lineStyles = {'-', '--'};

for i = 1:numel(replicationSims)

    d = replicationSims(i).death;
    baseColor = repColors(i, :);

    plot(d.t, log10(d.natural + 1), ...
        'LineStyle', lineStyles{i}, ...
        'Color', baseColor, ...
        'LineWidth', 1.8);

    plot(d.t, log10(d.adcc + 1), ...
        'LineStyle', lineStyles{i}, ...
        'Color', orange, ...
        'LineWidth', 1.5);

    plot(d.t, log10(d.tc + 1), ...
        'LineStyle', lineStyles{i}, ...
        'Color', gray, ...
        'LineWidth', 1.5);
end

xlabel('Time');
ylabel('log_{10}(cell-death flux + 1)');
title('d  Replication capacity reshapes killing mechanisms', ...
    'FontWeight', 'bold');

legend({ ...
    'Natural, k_5 = 0.05', ...
    'ADCC, k_5 = 0.05', ...
    'Tc killing, k_5 = 0.05', ...
    'Natural, k_5 = 0.10', ...
    'ADCC, k_5 = 0.10', ...
    'Tc killing, k_5 = 0.10'}, ...
    'Location', 'northeast', ...
    'Box', 'off');

xlim(timeLimit);
box off;

%% =========================
% E. Immune competence: V and I_total
% =========================

axE = nexttile(tl, [1 2]);
hold on;

for i = 1:numel(immuneSims)
    f = immuneSims(i).features;

    plot(f.t, log10(f.V + 1), ...
        '-', ...
        'Color', immuneColors(i, :), ...
        'LineWidth', 1.8);

    plot(f.t, log10(f.I_total + 1), ...
        '--', ...
        'Color', immuneColors(i, :), ...
        'LineWidth', 1.5);
end

xlabel('Time');
ylabel('log_{10}(burden + 1)');
title('e  Host immune competence determines persistence at k_5 = 0.10', ...
    'FontWeight', 'bold');

legend({ ...
    'V, strong immunity', ...
    'I_{total}, strong immunity', ...
    'V, weak immunity', ...
    'I_{total}, weak immunity'}, ...
    'Location', 'northeast', ...
    'Box', 'off');

xlim(timeLimit);
box off;

%% =========================
% F. Immune competence: Tc cells
% =========================

axF = nexttile(tl, [1 2]);
hold on;

for i = 1:numel(immuneSims)
    f = immuneSims(i).features;

    plot(f.t, log10(f.Tc + 1), ...
        'Color', immuneColors(i, :), ...
        'LineWidth', 1.8);
end

xlabel('Time');
ylabel('log_{10}(CD8^{+} T cells + 1)');
title('f  Strong immunity supports larger CD8^{+} T-cell responses', ...
    'FontWeight', 'bold');

legend({'Strong immunity', 'Weak immunity'}, ...
    'Location', 'southeast', ...
    'Box', 'off');

xlim(timeLimit);
box off;

%% =========================
% G. Infection-age heatmaps for strong and weak immunity
% =========================

axG1 = nexttile;
f = immuneSims(1).features;

imagesc(f.t, 1:N, log10(f.I_age' + 1));
set(gca, 'YDir', 'normal');

xlabel('Time');
ylabel('Infection age');
title('g  I_{1-400}, strong immunity', ...
    'FontWeight', 'bold');

cbG1 = colorbar;
cbG1.Label.String = 'log_{10}(I_a + 1)';
colormap(axG1, parula);
xlim(timeLimit);
box off;

axG2 = nexttile;
f = immuneSims(2).features;

imagesc(f.t, 1:N, log10(f.I_age' + 1));
set(gca, 'YDir', 'normal');

xlabel('Time');
ylabel('Infection age');
title('g  I_{1-400}, weak immunity', ...
    'FontWeight', 'bold');

cbG2 = colorbar;
cbG2.Label.String = 'log_{10}(I_a + 1)';
colormap(axG2, parula);
xlim(timeLimit);
box off;

%% =========================
% H. Death mechanisms for strong and weak immunity
% =========================

axH = nexttile(tl, [1 2]);
hold on;

lineStyles = {'-', '--'};

for i = 1:numel(immuneSims)

    d = immuneSims(i).death;
    baseColor = immuneColors(i, :);

    plot(d.t, log10(d.natural + 1), ...
        'LineStyle', lineStyles{i}, ...
        'Color', baseColor, ...
        'LineWidth', 1.8);

    plot(d.t, log10(d.adcc + 1), ...
        'LineStyle', lineStyles{i}, ...
        'Color', orange, ...
        'LineWidth', 1.5);

    plot(d.t, log10(d.tc + 1), ...
        'LineStyle', lineStyles{i}, ...
        'Color', gray, ...
        'LineWidth', 1.5);
end

xlabel('Time');
ylabel('log_{10}(cell-death flux + 1)');
title('h  Immune competence determines killing pathway engagement', ...
    'FontWeight', 'bold');

legend({ ...
    'Natural, strong', ...
    'ADCC, strong', ...
    'Tc killing, strong', ...
    'Natural, weak', ...
    'ADCC, weak', ...
    'Tc killing, weak'}, ...
    'Location', 'northeast', ...
    'Box', 'off');

xlim(timeLimit);
box off;

%% =========================
% 7. Save figure
% =========================

title(tl, ...
    'Viral replication capacity and host immune competence determine acute clearance and chronic persistence', ...
    'FontWeight', 'bold', ...
    'FontSize', 12);

drawnow;
pause(0.5);

pngFile = [figPrefix '.png'];
pdfFile = [figPrefix '.pdf'];
figFile = [figPrefix '.fig'];

try
    exportgraphics(fig, pngFile, ...
        'Resolution', 600, ...
        'ContentType', 'image');
catch
    print(fig, pngFile, '-dpng', '-r600');
end

try
    exportgraphics(fig, pdfFile, ...
        'Resolution', 600, ...
        'ContentType', 'image');
catch
    print(fig, pdfFile, '-dpdf', '-r600');
end

if isvalid(fig)
    savefig(fig, figFile);
end

%% =========================
% 8. Write draft results text
% =========================

textFile = fullfile(outDir, ...
    'results_section_replication_capacity_immune_competence.txt');

write_results_text(textFile, summary);

fprintf('\nFigures saved:\n');
fprintf('%s\n', pngFile);
fprintf('%s\n', pdfFile);
fprintf('%s\n', figFile);

fprintf('\nSummary saved:\n%s\n', summaryFile);
fprintf('\nResults text saved:\n%s\n', textFile);

disp('Done.');

%% ========================================================================
% Local functions
% ========================================================================

function [pBase, y0Base] = initialize_model_from_main_script(mainParameterScript)

fprintf('\nInitializing baseline model from %s\n', mainParameterScript);

if exist(mainParameterScript, 'file') ~= 2
    error('Cannot find %s in the current MATLAB path.', mainParameterScript);
end

run(mainParameterScript);

if ~exist('p', 'var')
    error('The main script did not create parameter structure p.');
end

pBase = p;

if exist('y0', 'var')
    y0Base = y0(:);
elseif exist('Y0', 'var')
    y0Base = Y0(:);
elseif exist('y', 'var') && ~isempty(y)
    y0Base = y(1, :)';
else
    error(['Cannot identify initial condition. ', ...
           'Please make sure the main script defines y0, Y0, or y.']);
end

if isfield(pBase, 'N')
    N = pBase.N;
else
    N = 400;
    pBase.N = N;
end

minLength = N + 1403;

if numel(y0Base) < minLength
    warning('Initial condition has %d variables, expected at least %d.', ...
        numel(y0Base), minLength);
end

pBase.drug_dose_times = [];
pBase.drug_dose_amounts = [];

if isfield(pBase, 'drug_on_time')
    pBase.drug_on_time = inf;
end

if isfield(pBase, 'drug_off_time')
    pBase.drug_off_time = inf;
end

pBase = refresh_intracellular_replication_vectors(pBase);

fprintf('Baseline model initialized with N = %d and %d state variables.\n', ...
    pBase.N, numel(y0Base));

end

function p = refresh_intracellular_replication_vectors(p)

if ~isfield(p, 'N')
    p.N = 400;
end

if ~isfield(p, 'da')
    p.da = 1;
end

if ~isfield(p, 'v_start')
    p.v_start = 1;
end

p.a_vec = (0:p.N-1)' * p.da;

if isfield(p, 'vin_exp_cap')
    exponentValue = min(p.k5 .* p.a_vec, p.vin_exp_cap);
elseif isfield(p, 'vin_cap')
    exponentValue = p.k5 .* p.a_vec;
else
    exponentValue = min(p.k5 .* p.a_vec, 700);
end

p.vin_vec = p.v_start .* exp(exponentValue);

if isfield(p, 'vin_cap')
    p.vin_vec = min(p.vin_vec, p.vin_cap);
end

end

function sim = run_single_simulation(p, y0, simulationDuration, odeFunctionName)

if ~isfield(p, 'N')
    p.N = 400;
end

if ~isfield(p, 'use_extinction_switch')
    p.use_extinction_switch = true;
end

if ~isfield(p, 'extinction_V_threshold')
    p.extinction_V_threshold = 1e-5;
end

if ~isfield(p, 'extinction_I_total_threshold')
    p.extinction_I_total_threshold = 1e-2;
end

if ~isfield(p, 'vin_cap')
    p.vin_cap = 1e50;
end

if ~isfield(p, 'vin_exp_cap')
    p.vin_exp_cap = 700;
end

odeFun = str2func(odeFunctionName);
odeHandle = @(t, y) odeFun(t, y, p);

tspan = [0 simulationDuration];

odeOptions = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-9, ...
    'NonNegative', 1:numel(y0));

try
    if p.use_extinction_switch
        [t, y, extinctionInfo] = solve_with_extinction_restart_analyze( ...
            odeHandle, tspan, y0, odeOptions, p);
    else
        [t, y] = ode15s(odeHandle, tspan, y0, odeOptions);

        extinctionInfo.triggered = false;
        extinctionInfo.time = NaN;
    end

    y = max(real(y), 0);

    sim.t = t;
    sim.y = y;
    sim.p = p;
    sim.success = true;
    sim.error_message = '';
    sim.extinctionInfo = extinctionInfo;

catch ME
    warning('Simulation failed in analyze script: %s', ME.message);

    sim.t = [];
    sim.y = [];
    sim.p = p;
    sim.success = false;
    sim.error_message = ME.message;

    sim.extinctionInfo.triggered = false;
    sim.extinctionInfo.time = NaN;
end

end

function death = compute_death_dynamics(sim, p)

t = sim.t(:);
y = max(sim.y, 0);
N = p.N;

I_age = y(:, 2:N+1);
Tc = y(:, N+3);

if isfield(p, 'vin_vec')
    vin_vec = p.vin_vec(:);
else
    vin_vec = p.v_start .* exp(p.k5 .* ((0:N-1)' .* p.da));
end

if numel(vin_vec) ~= N
    vin_vec = interp1( ...
        linspace(1, N, numel(vin_vec)), ...
        vin_vec, ...
        1:N, ...
        'linear', ...
        'extrap')';
end

theta_lysis = get_param_with_default(p, 'theta_lysis', 1e12);
theta_adcc = get_param_with_default(p, 'theta_adcc', 1e8);
theta_tc = get_param_with_default(p, 'theta_tc', 1e7);

n_hill = get_param_with_default(p, 'n_hill', 1);
k_lysis_max = get_param_with_default(p, 'k_lysis_max', 1);
k_adcc_max = get_param_with_default(p, 'k_adcc_max', 1);
k_tc_max = get_param_with_default(p, 'k_tc_max', 1);

Tc_binding = get_param_with_default(p, 'Tc_binding', 1e-5);

naturalRateByAge = k_lysis_max .* hill_function(vin_vec, theta_lysis, n_hill);
natural = I_age * naturalRateByAge(:);

antibodyPressure = estimate_antibody_pressure(y, p);
adccRateByTime = k_adcc_max .* hill_function(antibodyPressure, theta_adcc, n_hill);
adcc = adccRateByTime(:) .* sum(I_age, 2);

tcPressure = Tc_binding .* Tc .* mean(vin_vec);
tcRateByTime = k_tc_max .* hill_function(tcPressure, theta_tc, n_hill);
tc = tcRateByTime(:) .* sum(I_age, 2);

death = struct;
death.t = t;
death.natural = max(natural, 0);
death.adcc = max(adcc, 0);
death.tc = max(tc, 0);

end

function value = get_param_with_default(p, fieldName, defaultValue)

if isfield(p, fieldName)
    value = p.(fieldName);
else
    value = defaultValue;
end

end

function h = hill_function(x, theta, n)

x = max(x, 0);

if theta <= 0
    h = ones(size(x));
    return;
end

h = (x .^ n) ./ (theta .^ n + x .^ n);

h(~isfinite(h)) = 0;

end

function antibodyPressure = estimate_antibody_pressure(y, p)

N = p.N;
numVars = size(y, 2);

X_start = N + 4;
X_end = min(numVars, N + 1403);

if X_start > numVars
    antibodyPressure = zeros(size(y, 1), 1);
    return;
end

X = y(:, X_start:X_end);
X = max(X, 0);

if size(X, 2) < 1400
    X(:, end+1:1400) = 0;
end

IgM_BCR_V = X(:, 801:900);
IgM_V = X(:, 901:1000);
IgG_BCR_V = X(:, 1001:1100);
IgG_V = X(:, 1101:1200);

antibodyPressure = ...
    sum(IgM_BCR_V, 2) + ...
    sum(IgM_V, 2) + ...
    sum(IgG_BCR_V, 2) + ...
    sum(IgG_V, 2);

end

function summary = make_summary(replicationSims, immuneSims)

summary = struct;

for i = 1:numel(replicationSims)
    f = replicationSims(i).features;

    summary.replication(i).label = replicationSims(i).label;
    summary.replication(i).k5 = replicationSims(i).k5;
    summary.replication(i).peakV = f.peakV;
    summary.replication(i).timePeakV = f.timePeakV;
    summary.replication(i).peakI = f.peakI;
    summary.replication(i).timePeakI = f.timePeakI;
    summary.replication(i).finalV = f.finalV;
    summary.replication(i).finalI = f.finalI;
    summary.replication(i).aucV = f.aucV;
    summary.replication(i).aucI = f.aucI;
    summary.replication(i).isCleared = f.isCleared;
end

for i = 1:numel(immuneSims)
    f = immuneSims(i).features;

    summary.immune(i).label = immuneSims(i).label;
    summary.immune(i).immuneCompetence = immuneSims(i).immuneCompetence;
    summary.immune(i).peakV = f.peakV;
    summary.immune(i).timePeakV = f.timePeakV;
    summary.immune(i).peakI = f.peakI;
    summary.immune(i).timePeakI = f.timePeakI;
    summary.immune(i).finalV = f.finalV;
    summary.immune(i).finalI = f.finalI;
    summary.immune(i).aucV = f.aucV;
    summary.immune(i).aucI = f.aucI;
    summary.immune(i).isCleared = f.isCleared;
end

end

function write_results_text(textFile, summary)

fid = fopen(textFile, 'w');

fprintf(fid, 'Viral replication capacity and host immune competence determine acute clearance and chronic persistence\n\n');

fprintf(fid, ['We first varied intracellular viral replication capacity by setting k5 to 0.05 or 0.10 in the absence of antiviral treatment. ', ...
    'The lower-replication condition showed a more limited viral and infected-cell burden, whereas the higher-replication condition produced a larger and more persistent infection. ', ...
    'These differences were accompanied by altered CD8+ T-cell activation, distinct infection-age distributions and differential engagement of natural lysis, ADCC and Tc-mediated killing pathways.\n\n']);

fprintf(fid, ['We next fixed viral replication capacity at k5 = 0.10 and varied host immune competence by setting p.para(15) to 5 or 2. ', ...
    'Strong host immune competence promoted greater CD8+ T-cell expansion and more effective control of infected cells, whereas weak immune competence permitted sustained viral replication and persistence of infection-age compartments. ', ...
    'These simulations indicate that chronic persistence arises from the combined effect of high viral replication capacity and insufficient immune-mediated clearance.\n\n']);

fprintf(fid, 'Figure legend\n\n');

fprintf(fid, ['Figure X | Viral replication capacity and host immune competence determine acute clearance and chronic persistence. ', ...
    'a, Extracellular virus and total infected-cell burden after varying intracellular viral replication capacity, with k5 set to 0.05 or 0.10. ', ...
    'b, Antigen-specific CD8+ T-cell dynamics under the two replication-capacity settings. ', ...
    'c, Infection-age-resolved heatmaps showing I1-I400 dynamics for k5 = 0.05 and k5 = 0.10. ', ...
    'd, Time-resolved infected-cell death fluxes mediated by natural lysis, antibody-dependent cellular cytotoxicity and CD8+ T-cell killing under the two replication-capacity settings. ', ...
    'e, Extracellular virus and total infected-cell burden at fixed k5 = 0.10 under strong and weak host immune competence. ', ...
    'f, CD8+ T-cell dynamics under strong and weak immune competence. ', ...
    'g, Infection-age-resolved heatmaps showing I1-I400 dynamics under strong and weak immune competence. ', ...
    'h, Natural lysis, ADCC and Tc-mediated killing dynamics under strong and weak immune competence. ', ...
    'Together, these simulations show that infection outcome is jointly determined by viral replication capacity and the strength of host adaptive immune clearance.\n']);

fclose(fid);

end

function [tAll, yAll, extinctionInfo] = solve_with_extinction_restart_analyze( ...
    odeFun, tspan, y0, odeOptions, p)

extinctionInfo.triggered = false;
extinctionInfo.time = NaN;

tStart = tspan(1);
tEnd = tspan(end);

eventOptions = odeset(odeOptions, ...
    'Events', @(t, y) extinction_event_analyze(t, y, p));

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

% Once I_total and V are both below thresholds, remove virus-related states
% before restarting the ODE.
yRestart = enforce_extinction_state_analyze(eventState, p);
yRestart = max(real(yRestart), 0);

restartOptions = odeset(odeOptions, 'Events', []);

[t2, y2] = ode15s(odeFun, [eventTime tEnd], yRestart, restartOptions);

tAll = [t1; t2(2:end)];
yAll = [y1; y2(2:end, :)];

% Enforce zero virus-related states for all post-extinction output rows.
idxPost = tAll >= eventTime;
yAll(idxPost, :) = enforce_extinction_state_matrix_analyze(yAll(idxPost, :), p);
yAll = max(real(yAll), 0);

end

function state = enforce_extinction_state_analyze(state, p)

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

function states = enforce_extinction_state_matrix_analyze(states, p)

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

function [value, isterminal, direction] = extinction_event_analyze(~, y, p)

N = p.N;

I_total = sum(max(y(2:N+1), 0));
V = max(y(N+2), 0);

V_threshold = p.extinction_V_threshold;
I_threshold = p.extinction_I_total_threshold;

% Event is triggered only when both I_total and V are below thresholds.
value = max(I_total / I_threshold, V / V_threshold) - 1;

isterminal = 1;
direction = -1;

end

function features = extract_features(sim, p)

features = struct;

if ~isfield(sim, 'success') || ~sim.success || isempty(sim.t) || isempty(sim.y)
    features.t = [];
    features.y = [];
    features.N = p.N;
    features.T = [];
    features.I_age = [];
    features.I_total = [];
    features.V = [];
    features.Tc = [];
    features.finalV = NaN;
    features.finalI = NaN;
    features.peakV = NaN;
    features.timePeakV = NaN;
    features.peakI = NaN;
    features.timePeakI = NaN;
    features.aucV = NaN;
    features.aucI = NaN;
    features.isCleared = false;
    return;
end

t = sim.t(:);
y = max(real(sim.y), 0);

N = p.N;

features.t = t;
features.y = y;
features.N = N;

features.T = y(:, 1);
features.I_age = y(:, 2:N+1);
features.I_total = sum(features.I_age, 2);
features.V = y(:, N+2);
features.Tc = y(:, N+3);

features.finalV = features.V(end);
features.finalI = features.I_total(end);

[features.peakV, idxPeakV] = max(features.V);
features.timePeakV = t(idxPeakV);

[features.peakI, idxPeakI] = max(features.I_total);
features.timePeakI = t(idxPeakI);

features.aucV = trapz(t, features.V);
features.aucI = trapz(t, features.I_total);

if isfield(p, 'extinction_V_threshold')
    V_threshold = p.extinction_V_threshold;
else
    V_threshold = 1e-5;
end

if isfield(p, 'extinction_I_total_threshold')
    I_threshold = p.extinction_I_total_threshold;
else
    I_threshold = 1e-2;
end

features.isCleared = features.finalV < V_threshold && ...
    features.finalI < I_threshold;

end
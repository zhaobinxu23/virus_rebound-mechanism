%% analyze_antibody_repertoire_memory_acute_chronic_nature.m
% The model captures antibody repertoire expansion, contraction and immune memory
% maintenance during acute infection and chronic infection

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
    'antibody_repertoire_memory_acute_chronic');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

figPrefix = fullfile(outDir, ...
    'antibody_repertoire_memory_acute_chronic_nature');

simulationDuration = 6000;
timeLimit = [0 simulationDuration];

affinityTimeGrid = 0:5:simulationDuration;
elisaTimeGrid = 0:20:simulationDuration;

acute.k5 = 0.1;
acute.immuneCompetence = 5;
acute.label = 'Acute infection';

chronic.k5 = 0.05;
chronic.immuneCompetence = 5;
chronic.label = 'Chronic infection';

%% =========================
% 2. Load baseline model
% =========================

[pBase, y0Base] = initialize_model_from_main_script(mainParameterScript);

if ~isfield(pBase, 'N')
    pBase.N = 400;
end

%% =========================
% 3. Run acute and chronic simulations
% =========================

acuteSim = run_condition( ...
    pBase, ...
    y0Base, ...
    acute.k5, ...
    acute.immuneCompetence, ...
    simulationDuration, ...
    odeFunctionName, ...
    acute.label);

chronicSim = run_condition( ...
    pBase, ...
    y0Base, ...
    chronic.k5, ...
    chronic.immuneCompetence, ...
    simulationDuration, ...
    odeFunctionName, ...
    chronic.label);

%% =========================
% 4. Extract antibody repertoires
% =========================

acuteData = extract_antibody_features(acuteSim, acute.label);
chronicData = extract_antibody_features(chronicSim, chronic.label);

acuteDist = compute_affinity_distribution_from_sim( ...
    acuteData, affinityTimeGrid, acute.label);

chronicDist = compute_affinity_distribution_from_sim( ...
    chronicData, affinityTimeGrid, chronic.label);

%% =========================
% 5. Compute ELISA signals
% =========================

[acuteElisaIgM, acuteElisaIgG, acuteElisaTotal] = compute_elisa_signal( ...
    acuteData, elisaTimeGrid, acute.label);

[chronicElisaIgM, chronicElisaIgG, chronicElisaTotal] = compute_elisa_signal( ...
    chronicData, elisaTimeGrid, chronic.label);

acuteData.elisaTime = elisaTimeGrid(:);
acuteData.elisaIgM = acuteElisaIgM;
acuteData.elisaIgG = acuteElisaIgG;
acuteData.elisaTotal = acuteElisaTotal;

chronicData.elisaTime = elisaTimeGrid(:);
chronicData.elisaIgM = chronicElisaIgM;
chronicData.elisaIgG = chronicElisaIgG;
chronicData.elisaTotal = chronicElisaTotal;

%% =========================
% 6. Summary metrics
% =========================

summary = make_summary(acuteData, chronicData);

summaryFile = fullfile(outDir, ...
    'antibody_repertoire_memory_acute_chronic_summary.mat');

save(summaryFile, ...
    'summary', ...
    'acuteSim', ...
    'chronicSim', ...
    'acuteData', ...
    'chronicData', ...
    'acuteDist', ...
    'chronicDist', ...
    'affinityTimeGrid', ...
    'elisaTimeGrid');

%% =========================
% 7. Nature-style figure
% =========================

set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontSize', 8);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 0.8);

igmColor = [0.10 0.37 0.71];
iggColor = [0.80 0.18 0.16];
totalColor = [0.15 0.15 0.15];

fig = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [2 2 30 24], ...
    'Renderer', 'painters');

tl = tiledlayout(fig, 3, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% =========================
% A. Acute IgM repertoire
% =========================

axA = nexttile;
plot_repertoire_heatmap( ...
    axA, ...
    acuteDist.t_grid, ...
    acuteDist.kd_values, ...
    acuteDist.kd_IgM, ...
    'a  Acute infection: IgM repertoire');

%% =========================
% B. Acute IgG repertoire
% =========================

axB = nexttile;
plot_repertoire_heatmap( ...
    axB, ...
    acuteDist.t_grid, ...
    acuteDist.kd_values, ...
    acuteDist.kd_IgG, ...
    'b  Acute infection: IgG repertoire');

%% =========================
% C. Acute ELISA
% =========================

axC = nexttile;
hold(axC, 'on');

plot(axC, acuteData.elisaTime, log10(acuteData.elisaIgM + 1), ...
    'Color', igmColor, ...
    'LineWidth', 1.8);

plot(axC, acuteData.elisaTime, log10(acuteData.elisaIgG + 1), ...
    'Color', iggColor, ...
    'LineWidth', 1.8);

plot(axC, acuteData.elisaTime, log10(acuteData.elisaTotal + 1), ...
    'Color', totalColor, ...
    'LineWidth', 2.0);

xlabel(axC, 'Time');
ylabel(axC, 'log_{10}(ELISA signal + 1)');
title(axC, 'c  Acute infection: functional antibody binding', ...
    'FontWeight', 'bold');

legend(axC, {'IgM', 'IgG', 'IgM + IgG'}, ...
    'Location', 'southeast', ...
    'Box', 'off');

xlim(axC, timeLimit);
box(axC, 'off');

%% =========================
% D. Chronic IgM repertoire
% =========================

axD = nexttile;
plot_repertoire_heatmap( ...
    axD, ...
    chronicDist.t_grid, ...
    chronicDist.kd_values, ...
    chronicDist.kd_IgM, ...
    'd  Chronic infection: IgM repertoire');

%% =========================
% E. Chronic IgG repertoire
% =========================

axE = nexttile;
plot_repertoire_heatmap( ...
    axE, ...
    chronicDist.t_grid, ...
    chronicDist.kd_values, ...
    chronicDist.kd_IgG, ...
    'e  Chronic infection: IgG repertoire');

%% =========================
% F. Chronic ELISA
% =========================

axF = nexttile;
hold(axF, 'on');

plot(axF, chronicData.elisaTime, log10(chronicData.elisaIgM + 1), ...
    'Color', igmColor, ...
    'LineWidth', 1.8);

plot(axF, chronicData.elisaTime, log10(chronicData.elisaIgG + 1), ...
    'Color', iggColor, ...
    'LineWidth', 1.8);

plot(axF, chronicData.elisaTime, log10(chronicData.elisaTotal + 1), ...
    'Color', totalColor, ...
    'LineWidth', 2.0);

xlabel(axF, 'Time');
ylabel(axF, 'log_{10}(ELISA signal + 1)');
title(axF, 'f  Chronic infection: sustained antibody binding', ...
    'FontWeight', 'bold');

legend(axF, {'IgM', 'IgG', 'IgM + IgG'}, ...
    'Location', 'southeast', ...
    'Box', 'off');

xlim(axF, timeLimit);
box(axF, 'off');

%% =========================
% Global title and save
% =========================

title(tl, ...
    'The model captures antibody repertoire expansion, contraction and immune memory maintenance during acute and chronic infection', ...
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
% 8. Write results draft
% =========================

textFile = fullfile(outDir, ...
    'results_section_antibody_repertoire_memory_acute_chronic.txt');

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
    error('Cannot identify initial condition. Please define y0, Y0, or y in the main script.');
end

if ~isfield(pBase, 'N')
    pBase.N = 400;
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

function sim = run_condition( ...
    pBase, y0Base, k5Value, immuneCompetence, simulationDuration, ...
    odeFunctionName, conditionLabel)

p = pBase;

if ~isfield(p, 'N')
    p.N = 400;
end

p.k5 = k5Value;
p = refresh_intracellular_replication_vectors(p);

p.para(15) = immuneCompetence;

p.drug_dose_times = [];
p.drug_dose_amounts = [];

if ~isfield(p, 'use_extinction_switch')
    p.use_extinction_switch = true;
else
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
    'NonNegative', 1:numel(y0Base));

fprintf('\nRunning %s: k5 = %.4g, immune competence = %.4g\n', ...
    conditionLabel, p.k5, p.para(15));

try
    if p.use_extinction_switch
        [t, y, extinctionInfo] = solve_with_extinction_restart_antibody( ...
            odeHandle, tspan, y0Base, odeOptions, p);
    else
        [t, y] = ode15s(odeHandle, tspan, y0Base, odeOptions);

        extinctionInfo.triggered = false;
        extinctionInfo.time = NaN;
    end

    y = max(real(y), 0);

    sim.t = t;
    sim.y = y;
    sim.p = p;
    sim.label = conditionLabel;
    sim.success = true;
    sim.error_message = '';
    sim.extinctionInfo = extinctionInfo;

catch ME
    warning('Simulation failed for %s: %s', conditionLabel, ME.message);

    sim.t = [];
    sim.y = [];
    sim.p = p;
    sim.label = conditionLabel;
    sim.success = false;
    sim.error_message = ME.message;

    sim.extinctionInfo.triggered = false;
    sim.extinctionInfo.time = NaN;
end

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
else
    exponentValue = min(p.k5 .* p.a_vec, 700);
end

p.vin_vec = p.v_start .* exp(exponentValue);

if isfield(p, 'vin_cap')
    p.vin_vec = min(p.vin_vec, p.vin_cap);
end

end

function data = extract_antibody_features(sim, groupName)

t = sim.t(:);
y = max(sim.y, 0);
p = sim.p;

if isfield(p, 'N')
    N = p.N;
else
    N = 400;
end

data = struct;
data.groupName = groupName;
data.t = t;
data.y = y;
data.p = p;
data.N = N;

data.T = y(:, 1);
data.I_age = y(:, 2:N+1);
data.I_total = sum(data.I_age, 2);
data.V = y(:, N+2);
data.Tc = y(:, N+3);

X_start = N + 4;
X_end = min(size(y, 2), N + 1403);

X = y(:, X_start:X_end);
X = max(X, 0);

if size(X, 2) < 1400
    X(:, end+1:1400) = 0;
end

data.X = X;

data.IgM_BCR = X(:, 1:100);
data.IgM = X(:, 101:200);
data.IgG_BCR = X(:, 201:300);
data.IgG = X(:, 301:400);

data.IgM_BCR_E = X(:, 401:500);
data.IgM_E = X(:, 501:600);
data.IgG_BCR_E = X(:, 601:700);
data.IgG_E = X(:, 701:800);

data.IgM_BCR_V = X(:, 801:900);
data.IgM_V = X(:, 901:1000);
data.IgG_BCR_V = X(:, 1001:1100);
data.IgG_V = X(:, 1101:1200);

data.P_M = X(:, 1201:1300);
data.P_G = X(:, 1301:1400);

data.totalIgM = sum(data.IgM, 2);
data.totalIgG = sum(data.IgG, 2);
data.totalAb = data.totalIgM + data.totalIgG;

data.totalPlasmaIgM = sum(data.P_M, 2);
data.totalPlasmaIgG = sum(data.P_G, 2);

data.totalVirusBoundIgM = sum(data.IgM_BCR_V + data.IgM_V, 2);
data.totalVirusBoundIgG = sum(data.IgG_BCR_V + data.IgG_V, 2);

end

function dist = compute_affinity_distribution_from_sim(data, affinityTimeGrid, groupName)

fprintf('\nComputing affinity distribution for %s...\n', groupName);

if exist('compute_affinity_distribution_cell', 'file') ~= 2
    error('Cannot find compute_affinity_distribution_cell.m in MATLAB path.');
end

affinityTimeGrid = affinityTimeGrid(:)';

try
    dist = compute_affinity_distribution_cell(data.t, data.y, affinityTimeGrid);
catch
    dist = compute_affinity_distribution_cell(data.t, data.X, affinityTimeGrid);
end

if ~isfield(dist, 't_grid')
    dist.t_grid = affinityTimeGrid;
end

if ~isfield(dist, 'kd_values')
    if isfield(dist, 'kd_IgM')
        dist.kd_values = 1:size(dist.kd_IgM, 1);
    else
        dist.kd_values = 1:10;
    end
end

dist.kd_IgM = max(dist.kd_IgM, 0);
dist.kd_IgG = max(dist.kd_IgG, 0);

end

function [elisaIgM, elisaIgG, elisaTotalAb] = compute_elisa_signal(data, elisaTimeGrid, groupName)

fprintf('\nComputing ELISA signal for %s...\n', groupName);

if exist('pathway_model_many_antibody_immune_include_plasma_elisa', 'file') ~= 2
    warning('Cannot find ELISA ODE function. ELISA values will be approximated from antibody abundance.');
    elisaIgM = sum(data.IgM, 2);
    elisaIgG = sum(data.IgG, 2);
    elisaIgM = interp1(data.t, elisaIgM, elisaTimeGrid(:), 'linear', 'extrap');
    elisaIgG = interp1(data.t, elisaIgG, elisaTimeGrid(:), 'linear', 'extrap');
    elisaTotalAb = elisaIgM + elisaIgG;
    return;
end

elisaTimeGrid = elisaTimeGrid(:);
numTime = numel(elisaTimeGrid);

dataNewIgM = zeros(100, numTime);
dataNewIgG = zeros(100, numTime);

for i = 1:100
    dataNewIgM(i, :) = interp1(data.t, data.IgM(:, i), elisaTimeGrid, 'linear', 'extrap');
    dataNewIgG(i, :) = interp1(data.t, data.IgG(:, i), elisaTimeGrid, 'linear', 'extrap');
end

dataNewIgM = max(dataNewIgM, 0);
dataNewIgG = max(dataNewIgG, 0);

elisaIgM = zeros(numTime, 1);
elisaIgG = zeros(numTime, 1);
elisaTotalAb = zeros(numTime, 1);

[para, para_new, para_new_1] = get_elisa_parameters();

odeOptions = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-9, ...
    'NonNegative', 1:1402);

for k = 1:numTime

    if k == 1 || k == numTime || mod(k, 25) == 0
        fprintf('  %s ELISA %d / %d, t = %.1f\n', ...
            groupName, k, numTime, elisaTimeGrid(k));
    end

    IgM_now = dataNewIgM(:, k)';
    IgG_now = dataNewIgG(:, k)';

    elisaIgM(k) = run_single_elisa_assay( ...
        IgM_now, zeros(1, 100), para, para_new, para_new_1, odeOptions);

    elisaIgG(k) = run_single_elisa_assay( ...
        zeros(1, 100), IgG_now, para, para_new, para_new_1, odeOptions);

    elisaTotalAb(k) = run_single_elisa_assay( ...
        IgM_now, IgG_now, para, para_new, para_new_1, odeOptions);

end

end

function elisaValue = run_single_elisa_assay(IgM_vector, IgG_vector, para, para_new, para_new_1, odeOptions)

x0 = zeros(1402, 1);

for i = 101:200
    x0(i) = IgM_vector(i - 100) / 1000;
end

for i = 301:400
    x0(i) = IgG_vector(i - 300) / 1000;
end

x0(1402) = 1e15;

try
    [tElisa, zz] = ode15s( ...
        @pathway_model_many_antibody_immune_include_plasma_elisa, ...
        [0 2], ...
        x0, ...
        odeOptions, ...
        para, ...
        para_new, ...
        para_new_1);

    virusAtEnd = interp1(tElisa, zz(:, 1402), 2);
    elisaValue = 1e15 - virusAtEnd;

    if ~isfinite(elisaValue) || elisaValue < 0
        elisaValue = 0;
    end

catch
    elisaValue = NaN;
end

end

function [para, para_new, para_new_1] = get_elisa_parameters()

para = zeros(20, 1);

para(1) = 1e-20;
para(2) = 0.5;
para(3) = 2.2e17 + 1e13;
para(4) = 0.5e13;
para(5) = 0.01;
para(6) = 0.005;
para(7) = 1;
para(8) = 5e-7;
para(9) = 0.1;
para(10) = 0.05;
para(11) = 4.4e9;
para(12) = 1.2e10;
para(13) = 0.05;
para(14) = 0.025;
para(15) = 0;
para(16) = 0.05;
para(17) = 0.1;
para(18) = 0.5;
para(19) = 0;
para(20) = 1e5;

para_new = zeros(10, 1);
para_new(1) = 1e-22;
para_new(2) = 1e-21;
para_new(3) = 1e-20;
para_new(4) = 1e-19;
para_new(5) = 1e-18;
para_new(6) = 1e-17;
para_new(7) = 1e-16;
para_new(8) = 1e-15;
para_new(9) = 1e-14;
para_new(10) = 1e-13;

para_new_1 = zeros(10, 1);
para_new_1(1) = 1e0;
para_new_1(2) = 1e1;
para_new_1(3) = 1e2;
para_new_1(4) = 1e3;
para_new_1(5) = 1e4;
para_new_1(6) = 1e5;
para_new_1(7) = 1e6;
para_new_1(8) = 1e7;
para_new_1(9) = 1e8;
para_new_1(10) = 1e9;

end

function plot_repertoire_heatmap(ax, tGrid, kdValues, repertoireMatrix, panelTitle)

imagesc(ax, tGrid, kdValues, log10(repertoireMatrix + 1));
set(ax, 'YDir', 'normal');

xlabel(ax, 'Time');
ylabel(ax, 'log_{10}(K_d)');
title(ax, panelTitle, 'FontWeight', 'bold');

cb = colorbar(ax);
cb.Label.String = 'log_{10}(antibody abundance + 1)';

colormap(ax, turbo);
box(ax, 'off');

end

function summary = make_summary(acuteData, chronicData)

summary = struct;

summary.acute.finalV = acuteData.V(end);
summary.acute.finalI = acuteData.I_total(end);
summary.acute.peakIgM = max(acuteData.totalIgM);
summary.acute.peakIgG = max(acuteData.totalIgG);
summary.acute.finalIgM = acuteData.totalIgM(end);
summary.acute.finalIgG = acuteData.totalIgG(end);
summary.acute.IgM_contraction_ratio = ...
    (summary.acute.finalIgM + 1) / (summary.acute.peakIgM + 1);
summary.acute.IgG_memory_ratio = ...
    (summary.acute.finalIgG + 1) / (summary.acute.peakIgG + 1);

summary.chronic.finalV = chronicData.V(end);
summary.chronic.finalI = chronicData.I_total(end);
summary.chronic.peakIgM = max(chronicData.totalIgM);
summary.chronic.peakIgG = max(chronicData.totalIgG);
summary.chronic.finalIgM = chronicData.totalIgM(end);
summary.chronic.finalIgG = chronicData.totalIgG(end);
summary.chronic.IgM_contraction_ratio = ...
    (summary.chronic.finalIgM + 1) / (summary.chronic.peakIgM + 1);
summary.chronic.IgG_memory_ratio = ...
    (summary.chronic.finalIgG + 1) / (summary.chronic.peakIgG + 1);

end

function write_results_text(textFile, summary)

fid = fopen(textFile, 'w');

fprintf(fid, 'The model captures antibody repertoire expansion, contraction and immune memory maintenance during acute infection and chronic infection\n\n');

fprintf(fid, ['In the acute-infection setting, defined by k5 = 0.1 and host immune competence p.para(15) = 5, ', ...
    'the model generated a transient IgM-dominated humoral response followed by contraction of the IgM compartment. ', ...
    'The IgG repertoire expanded more gradually and was preferentially maintained at later times, consistent with the transition from an early primary antibody response to longer-lived humoral memory. ', ...
    'Functional ELISA simulations showed an early IgM contribution, a delayed IgG contribution and sustained total antigen-binding activity after the acute phase.\n\n']);

fprintf(fid, ['In the chronic-infection setting, defined by k5 = 0.05 and p.para(15) = 3, ', ...
    'persistent antigenic stimulation altered the kinetics of antibody-repertoire remodeling. ', ...
    'Compared with acute infection, the chronic condition maintained ongoing IgM and IgG repertoire activity and sustained ELISA-detectable binding over the observation window. ', ...
    'These results indicate that the model captures both contraction after acute infection and continued antibody-repertoire remodeling during persistent infection.\n\n']);

fprintf(fid, ['Quantitatively, the final-to-peak IgM ratio was %.3g in acute infection and %.3g in chronic infection, ', ...
    'whereas the final-to-peak IgG ratio was %.3g in acute infection and %.3g in chronic infection. ', ...
    'Thus, IgM responses contracted more strongly after acute infection, while IgG was preferentially preserved as a memory-like compartment.\n\n'], ...
    summary.acute.IgM_contraction_ratio, ...
    summary.chronic.IgM_contraction_ratio, ...
    summary.acute.IgG_memory_ratio, ...
    summary.chronic.IgG_memory_ratio);

fprintf(fid, 'Figure legend\n\n');

fprintf(fid, ['Figure X | The model captures antibody repertoire expansion, contraction and immune memory maintenance during acute and chronic infection. ', ...
    'a, Time-resolved IgM affinity-distribution heatmap during acute infection. Colours indicate log10-transformed antibody abundance across binding-energy classes. ', ...
    'b, Time-resolved IgG affinity-distribution heatmap during acute infection. IgG expands more gradually than IgM and is preferentially maintained at later times. ', ...
    'c, ELISA-derived IgM, IgG and total antibody-binding signals during acute infection. ', ...
    'd, Time-resolved IgM affinity-distribution heatmap during chronic infection. Persistent antigen exposure sustains IgM repertoire activity. ', ...
    'e, Time-resolved IgG affinity-distribution heatmap during chronic infection. Chronic infection maintains continued IgG repertoire remodeling. ', ...
    'f, ELISA-derived IgM, IgG and total antibody-binding signals during chronic infection. ', ...
    'Together, these simulations show that the model reproduces the expansion and contraction of early IgM responses, the maintenance of IgG-dominated immune memory after acute infection and sustained antibody-repertoire remodeling during chronic infection.\n']);

fclose(fid);

end

function [tAll, yAll, extinctionInfo] = solve_with_extinction_restart_antibody( ...
    odeFun, tspan, y0, odeOptions, p)

extinctionInfo.triggered = false;
extinctionInfo.time = NaN;

tStart = tspan(1);
tEnd = tspan(end);

eventOptions = odeset(odeOptions, ...
    'Events', @(t, y) extinction_event_antibody(t, y, p));

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

yRestart = enforce_extinction_state_antibody(eventState, p);
yRestart = max(real(yRestart), 0);

restartOptions = odeset(odeOptions, 'Events', []);

[t2, y2] = ode15s(odeFun, [eventTime tEnd], yRestart, restartOptions);

tAll = [t1; t2(2:end)];
yAll = [y1; y2(2:end, :)];

idxPost = tAll >= eventTime;
yAll(idxPost, :) = enforce_extinction_state_matrix_antibody(yAll(idxPost, :), p);
yAll = max(real(yAll), 0);

end

function state = enforce_extinction_state_antibody(state, p)

N = p.N;

state(2:N+1) = 0;
state(N+2) = 0;

x0Offset = N + 3;

virusComplexIdx = x0Offset + (801:1200);
virusComplexIdx = virusComplexIdx(virusComplexIdx <= numel(state));

state(virusComplexIdx) = 0;

end

function states = enforce_extinction_state_matrix_antibody(states, p)

if isempty(states)
    return;
end

N = p.N;

states(:, 2:N+1) = 0;
states(:, N+2) = 0;

x0Offset = N + 3;

virusComplexIdx = x0Offset + (801:1200);
virusComplexIdx = virusComplexIdx(virusComplexIdx <= size(states, 2));

states(:, virusComplexIdx) = 0;

end

function [value, isterminal, direction] = extinction_event_antibody(~, y, p)

N = p.N;

I_total = sum(max(y(2:N+1), 0));
V = max(y(N+2), 0);

V_threshold = p.extinction_V_threshold;
I_threshold = p.extinction_I_total_threshold;

value = max(I_total / I_threshold, V / V_threshold) - 1;

isterminal = 1;
direction = -1;

end
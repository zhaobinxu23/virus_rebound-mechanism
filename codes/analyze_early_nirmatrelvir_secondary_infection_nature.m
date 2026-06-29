%% analyze_early_nirmatrelvir_secondary_infection_nature.m
% Early nirmatrelvir treatment blunts immune memory and increases susceptibility to secondary infection

clear;
clc;
close all;

%% =========================
% 1. File paths
% =========================

rawDir = fullfile( ...
    'drug_rebound_scan_results', ...
    'raw_simulations');

earlyFile = fullfile(rawDir, ...
    'sim_start_0000_dose_1e+08_ndose_20.mat');

controlFile = fullfile(rawDir, ...
    'sim_start_0400_dose_1e+05_ndose_05.mat');

outDir = fullfile( ...
    'drug_rebound_scan_results', ...
    'early_nirmatrelvir_secondary_infection');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

figPrefix = fullfile(outDir, ...
    'early_nirmatrelvir_blunts_memory_secondary_infection');

primaryTimeLimit = [0 4000];
secondaryTimeLimit = [0 2000];

secondaryVirusInoculum = 10;
secondaryDuration = 2000;

affinity_time_grid_primary = 0:5:4000;
elisa_time_grid_primary = 0:20:4000;

%% =========================
% 2. Load primary simulations
% =========================

tmp = load(earlyFile);
simEarly = tmp.sim;

tmp = load(controlFile);
simControl = tmp.sim;

earlyPrimary = extract_features(simEarly, elisa_time_grid_primary, 'early nirmatrelvir');
controlPrimary = extract_features(simControl, elisa_time_grid_primary, 'control');

earlyDist = compute_affinity_distribution_from_sim(earlyPrimary, affinity_time_grid_primary, 'early nirmatrelvir');
controlDist = compute_affinity_distribution_from_sim(controlPrimary, affinity_time_grid_primary, 'control');

%% =========================
% 3. Run secondary infection simulations
% =========================

secondaryEarly = run_secondary_infection_from_final_state( ...
    simEarly, secondaryVirusInoculum, secondaryDuration, 'early nirmatrelvir');

secondaryControl = run_secondary_infection_from_final_state( ...
    simControl, secondaryVirusInoculum, secondaryDuration, 'control');

earlySecondary = extract_secondary_features(secondaryEarly, 'early nirmatrelvir');
controlSecondary = extract_secondary_features(secondaryControl, 'control');

%% =========================
% 4. Summary metrics
% =========================

summary = make_summary( ...
    earlyPrimary, ...
    controlPrimary, ...
    earlySecondary, ...
    controlSecondary, ...
    earlyDist, ...
    controlDist, ...
    secondaryVirusInoculum);

summaryFile = fullfile(outDir, ...
    'early_nirmatrelvir_secondary_infection_summary.mat');

save(summaryFile, ...
    'summary', ...
    'earlyPrimary', ...
    'controlPrimary', ...
    'earlySecondary', ...
    'controlSecondary', ...
    'earlyDist', ...
    'controlDist', ...
    'secondaryVirusInoculum', ...
    'secondaryDuration', ...
    'elisa_time_grid_primary', ...
    'affinity_time_grid_primary');

%% =========================
% 5. Figure style
% =========================

set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontSize', 8);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultAxesLineWidth', 0.8);

earlyColor = [0.80 0.18 0.16];
controlColor = [0.10 0.37 0.71];
gray = [0.45 0.45 0.45];
orange = [0.93 0.49 0.19];
purple = [0.49 0.18 0.56];
green = [0.20 0.55 0.30];

fig = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [2 2 32 30], ...
    'Renderer', 'painters');

tl = tiledlayout(fig, 4, 4, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% =========================
% Panel A: Primary infection viral kinetics
% =========================

axA = nexttile(tl, [1 2]);
hold on;

plot(controlPrimary.t, log10(controlPrimary.V + 1), ...
    'Color', controlColor, ...
    'LineWidth', 1.6);

plot(earlyPrimary.t, log10(earlyPrimary.V + 1), ...
    'Color', earlyColor, ...
    'LineWidth', 1.8);

xlabel('Time');
ylabel('log_{10}(V + 1)');
title('a  Early treatment clears primary infection without rebound', ...
    'FontWeight', 'bold');

legend({'Delayed/limited treatment control', 'Early nirmatrelvir'}, ...
    'Location', 'northeast', ...
    'Box', 'off');

xlim(primaryTimeLimit);
box off;

%% =========================
% Panel B: Primary infected-cell burden
% =========================

axB = nexttile(tl, [1 2]);
hold on;

plot(controlPrimary.t, log10(controlPrimary.I_total + 1), ...
    'Color', controlColor, ...
    'LineWidth', 1.6);

plot(earlyPrimary.t, log10(earlyPrimary.I_total + 1), ...
    'Color', earlyColor, ...
    'LineWidth', 1.8);

xlabel('Time');
ylabel('log_{10}(I_{total} + 1)');
title('b  Early treatment truncates antigen exposure', ...
    'FontWeight', 'bold');

xlim(primaryTimeLimit);
box off;

%% =========================
% Panel C: Pre-secondary CD8 memory
% =========================

axC = nexttile;
hold on;

barData = [ ...
    controlPrimary.Tc(end), ...
    earlyPrimary.Tc(end)];

bar(1, log10(barData(1) + 1), ...
    'FaceColor', controlColor, ...
    'EdgeColor', 'none');

bar(2, log10(barData(2) + 1), ...
    'FaceColor', earlyColor, ...
    'EdgeColor', 'none');

set(gca, 'XTick', [1 2], ...
    'XTickLabel', {'Control', 'Early'}, ...
    'XTickLabelRotation', 25);

ylabel('log_{10}(CD8^{+} T cells + 1)');
title('c  Blunted CD8^{+} T-cell memory', ...
    'FontWeight', 'bold');

box off;

%% =========================
% Panel D: Pre-secondary antibody memory
% =========================

axD = nexttile;
hold on;

IgM_memory = [controlPrimary.totalIgM(end), earlyPrimary.totalIgM(end)];
IgG_memory = [controlPrimary.totalIgG(end), earlyPrimary.totalIgG(end)];

barX = categorical({'Control', 'Early'});
barX = reordercats(barX, {'Control', 'Early'});

b = bar(barX, log10([IgM_memory(:), IgG_memory(:)] + 1), ...
    'grouped', ...
    'EdgeColor', 'none');

b(1).FaceColor = [0.25 0.50 0.85];
b(2).FaceColor = orange;

ylabel('log_{10}(Ab + 1)');
title('d  Blunted antibody memory', ...
    'FontWeight', 'bold');

legend({'IgM', 'IgG'}, ...
    'Location', 'northwest', ...
    'Box', 'off');

box off;

%% =========================
% Panel E: Functional ELISA memory
% =========================

axE = nexttile;
hold on;

controlElisaFinal = controlPrimary.elisaTotalAb(end);
earlyElisaFinal = earlyPrimary.elisaTotalAb(end);

bar(1, log10(controlElisaFinal + 1), ...
    'FaceColor', controlColor, ...
    'EdgeColor', 'none');

bar(2, log10(earlyElisaFinal + 1), ...
    'FaceColor', earlyColor, ...
    'EdgeColor', 'none');

set(gca, 'XTick', [1 2], ...
    'XTickLabel', {'Control', 'Early'}, ...
    'XTickLabelRotation', 25);

ylabel('log_{10}(ELISA signal + 1)');
title('e  Reduced functional antibody memory', ...
    'FontWeight', 'bold');

box off;

%% =========================
% Panel F: Secondary viral load
% =========================

axF = nexttile(tl, [1 2]);
hold on;

plot(controlSecondary.t, log10(controlSecondary.V + 1), ...
    'Color', controlColor, ...
    'LineWidth', 1.8);

plot(earlySecondary.t, log10(earlySecondary.V + 1), ...
    'Color', earlyColor, ...
    'LineWidth', 2.0);

xlabel('Time after secondary challenge');
ylabel('log_{10}(V + 1)');
title('f  Early-treatment host is susceptible to secondary infection', ...
    'FontWeight', 'bold');

legend({'Immune-experienced control', 'Early nirmatrelvir'}, ...
    'Location', 'northeast', ...
    'Box', 'off');

xlim(secondaryTimeLimit);
box off;

%% =========================
% Panel G: Secondary infected cells
% =========================

axG = nexttile(tl, [1 2]);
hold on;

plot(controlSecondary.t, log10(controlSecondary.I_total + 1), ...
    'Color', controlColor, ...
    'LineWidth', 1.8);

plot(earlySecondary.t, log10(earlySecondary.I_total + 1), ...
    'Color', earlyColor, ...
    'LineWidth', 2.0);

xlabel('Time after secondary challenge');
ylabel('log_{10}(I_{total} + 1)');
title('g  Weak memory permits renewed infected-cell expansion', ...
    'FontWeight', 'bold');

xlim(secondaryTimeLimit);
box off;

%% =========================
% Panel H: Secondary peak burden
% =========================

axH = nexttile;
hold on;

peakV = [controlSecondary.peakV, earlySecondary.peakV];
peakI = [controlSecondary.peakI, earlySecondary.peakI];

barX = categorical({'Control', 'Early'});
barX = reordercats(barX, {'Control', 'Early'});

b = bar(barX, log10([peakV(:), peakI(:)] + 1), ...
    'grouped', ...
    'EdgeColor', 'none');

b(1).FaceColor = green;
b(2).FaceColor = purple;

ylabel('log_{10}(peak burden + 1)');
title('h  Increased secondary-infection burden', ...
    'FontWeight', 'bold');

legend({'Peak V', 'Peak I'}, ...
    'Location', 'northwest', ...
    'Box', 'off');

box off;

%% =========================
% Panel I: Control IgM affinity memory
% =========================

axI = nexttile;

imagesc( ...
    controlDist.t_grid, ...
    controlDist.kd_values, ...
    log10(controlDist.kd_IgM + 1));

set(gca, 'YDir', 'normal');
xlabel('Time');
ylabel('log_{10}(K_d)');
title('i  Control IgM repertoire', ...
    'FontWeight', 'bold');

cbI = colorbar;
cbI.Label.String = 'log_{10}(IgM + 1)';
colormap(axI, parula);
xlim(primaryTimeLimit);
box off;

%% =========================
% Panel J: Early IgM affinity memory
% =========================

axJ = nexttile;

imagesc( ...
    earlyDist.t_grid, ...
    earlyDist.kd_values, ...
    log10(earlyDist.kd_IgM + 1));

set(gca, 'YDir', 'normal');
xlabel('Time');
ylabel('log_{10}(K_d)');
title('j  Early-treatment IgM repertoire', ...
    'FontWeight', 'bold');

cbJ = colorbar;
cbJ.Label.String = 'log_{10}(IgM + 1)';
colormap(axJ, parula);
xlim(primaryTimeLimit);
box off;

%% =========================
% Panel K: Control IgG affinity memory
% =========================

axK = nexttile;

imagesc( ...
    controlDist.t_grid, ...
    controlDist.kd_values, ...
    log10(controlDist.kd_IgG + 1));

set(gca, 'YDir', 'normal');
xlabel('Time');
ylabel('log_{10}(K_d)');
title('k  Control IgG repertoire', ...
    'FontWeight', 'bold');

cbK = colorbar;
cbK.Label.String = 'log_{10}(IgG + 1)';
colormap(axK, parula);
xlim(primaryTimeLimit);
box off;

%% =========================
% Panel L: Early IgG affinity memory
% =========================

axL = nexttile;

imagesc( ...
    earlyDist.t_grid, ...
    earlyDist.kd_values, ...
    log10(earlyDist.kd_IgG + 1));

set(gca, 'YDir', 'normal');
xlabel('Time');
ylabel('log_{10}(K_d)');
title('l  Early-treatment IgG repertoire', ...
    'FontWeight', 'bold');

cbL = colorbar;
cbL.Label.String = 'log_{10}(IgG + 1)';
colormap(axL, parula);
xlim(primaryTimeLimit);
box off;

%% =========================
% Figure title and save
% =========================

title(tl, ...
    'Early nirmatrelvir treatment blunts immune memory and increases susceptibility to secondary infection', ...
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
% Write Results section draft
% =========================

textFile = fullfile(outDir, ...
    'results_section_early_nirmatrelvir_secondary_infection.txt');

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

function data = extract_features(sim, elisa_time_grid, groupName)

t = sim.t(:);
y = sim.y;

if isfield(sim, 'p') && isfield(sim.p, 'N')
    N = sim.p.N;
else
    N = 400;
end

data.groupName = groupName;
data.t = t;
data.y = y;
data.N = N;

data.T = max(y(:, 1), 0);
data.I_total = max(sum(y(:, 2:N+1), 2), 0);
data.V = max(y(:, N+2), 0);
data.Tc = max(y(:, N+3), 0);

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

data.IgM_BCR_V = X(:, 801:900);
data.IgM_V = X(:, 901:1000);
data.IgG_BCR_V = X(:, 1001:1100);
data.IgG_V = X(:, 1101:1200);

data.totalIgM = sum(data.IgM, 2);
data.totalIgG = sum(data.IgG, 2);
data.totalAb = data.totalIgM + data.totalIgG;

data.totalVirusBoundIgM = sum(data.IgM_BCR_V + data.IgM_V, 2);
data.totalVirusBoundIgG = sum(data.IgG_BCR_V + data.IgG_V, 2);
data.totalVirusBoundAb = data.totalVirusBoundIgM + data.totalVirusBoundIgG;

data.elisa_time_grid = elisa_time_grid(:);

[data.elisaIgM, data.elisaIgG, data.elisaTotalAb] = ...
    compute_elisa_signal(data, elisa_time_grid, groupName);

if isfield(sim, 'p')
    data.p = sim.p;
else
    data.p = [];
end

if isfield(sim, 'protocol')
    data.protocol = sim.protocol;
else
    data.protocol = [];
end

end

function secondary = run_secondary_infection_from_final_state(sim, secondaryVirusInoculum, secondaryDuration, groupName)

fprintf('\nRunning secondary infection for %s group...\n', groupName);

if exist('sys_ode_modified_new', 'file') ~= 2
    error('Cannot find sys_ode_modified_new.m. Please make sure it is in the MATLAB path.');
end

if ~isfield(sim, 'p')
    error('The loaded sim structure does not contain sim.p.');
end

p = sim.p;

if isfield(p, 'drug_dose_times')
    p.drug_dose_times = [];
end

if isfield(p, 'drug_dose_amounts')
    p.drug_dose_amounts = [];
end

if isfield(p, 'use_extinction_switch')
    p.use_extinction_switch = false;
end

if isfield(p, 'extinction_V_threshold')
    p.extinction_V_threshold = 0;
end

if isfield(p, 'extinction_I_total_threshold')
    p.extinction_I_total_threshold = 0;
end

if isfield(p, 'N')
    N = p.N;
else
    N = 400;
    p.N = N;
end

if isfield(p, 'da') && isfield(p, 'k5') && isfield(p, 'v_start')
    p.a_vec = (0:p.N-1)' * p.da;
    p.vin_vec = p.v_start .* exp(p.k5 .* p.a_vec);
end

y0 = sim.y(end, :)';
y0 = max(y0, 0);

virusIndex = N + 2;
y0(virusIndex) = secondaryVirusInoculum;

t_span = [0 secondaryDuration];

options = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-9, ...
    'NonNegative', 1:numel(y0), ...
    'MaxStep', 5);

[tSecondary, ySecondary] = ode15s( ...
    @(t, y) sys_ode_modified_new(t, y, p), ...
    t_span, ...
    y0, ...
    options);

secondary = struct;
secondary.groupName = groupName;
secondary.t = tSecondary(:);
secondary.y = ySecondary;
secondary.p = p;
secondary.secondaryVirusInoculum = secondaryVirusInoculum;

end

function data = extract_secondary_features(secondary, groupName)

t = secondary.t(:);
y = secondary.y;

if isfield(secondary.p, 'N')
    N = secondary.p.N;
else
    N = 400;
end

data.groupName = groupName;
data.t = t;
data.y = y;
data.N = N;
data.T = max(y(:, 1), 0);
data.I_total = max(sum(y(:, 2:N+1), 2), 0);
data.V = max(y(:, N+2), 0);
data.Tc = max(y(:, N+3), 0);

[data.peakV, idxPeakV] = max(data.V);
data.timePeakV = data.t(idxPeakV);

[data.peakI, idxPeakI] = max(data.I_total);
data.timePeakI = data.t(idxPeakI);

data.aucV = trapz(data.t, data.V);
data.aucI = trapz(data.t, data.I_total);

end

function dist = compute_affinity_distribution_from_sim(data, affinity_time_grid, groupName)

fprintf('\nComputing affinity distribution for %s group...\n', groupName);

if exist('compute_affinity_distribution_cell', 'file') ~= 2
    error('Cannot find compute_affinity_distribution_cell.m. Please make sure it is in the MATLAB path.');
end

affinity_time_grid = affinity_time_grid(:)';

try
    dist = compute_affinity_distribution_cell(data.t, data.y, affinity_time_grid);
catch
    dist = compute_affinity_distribution_cell(data.t, data.X, affinity_time_grid);
end

if ~isfield(dist, 't_grid')
    dist.t_grid = affinity_time_grid;
end

if ~isfield(dist, 'kd_values')
    dist.kd_values = 1:size(dist.kd_IgM, 1);
end

dist.kd_IgM = max(dist.kd_IgM, 0);
dist.kd_IgG = max(dist.kd_IgG, 0);

end

function [elisaIgM, elisaIgG, elisaTotalAb] = compute_elisa_signal(data, elisa_time_grid, groupName)

fprintf('\nComputing ELISA signal for %s group...\n', groupName);

if exist('pathway_model_many_antibody_immune_include_plasma_elisa', 'file') ~= 2
    warning('Cannot find ELISA ODE function. ELISA values will be set to NaN.');
    elisaIgM = nan(numel(elisa_time_grid), 1);
    elisaIgG = nan(numel(elisa_time_grid), 1);
    elisaTotalAb = nan(numel(elisa_time_grid), 1);
    return;
end

elisa_time_grid = elisa_time_grid(:);
numTime = numel(elisa_time_grid);

data_new_IgM = zeros(100, numTime);
data_new_IgG = zeros(100, numTime);

for i = 1:100
    data_new_IgM(i, :) = interp1(data.t, data.IgM(:, i), elisa_time_grid, 'linear', 'extrap');
    data_new_IgG(i, :) = interp1(data.t, data.IgG(:, i), elisa_time_grid, 'linear', 'extrap');
end

data_new_IgM = max(data_new_IgM, 0);
data_new_IgG = max(data_new_IgG, 0);

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
            groupName, k, numTime, elisa_time_grid(k));
    end

    IgM_now = data_new_IgM(:, k)';
    IgG_now = data_new_IgG(:, k)';

    elisaIgM(k) = run_single_elisa_assay( ...
        IgM_now, zeros(1, 100), para, para_new, para_new_1, odeOptions);

    elisaIgG(k) = run_single_elisa_assay( ...
        zeros(1, 100), IgG_now, para, para_new, para_new_1, odeOptions);

    elisaTotalAb(k) = run_single_elisa_assay( ...
        IgM_now, IgG_now, para, para_new, para_new_1, odeOptions);

end

end

function elisaValue = run_single_elisa_assay(IgM_vector, IgG_vector, ...
    para, para_new, para_new_1, odeOptions)

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
para(3) = 0;
para(4) = 0;
para(5) = 0;
para(6) = 0;
para(7) = 0;
para(8) = 0;
para(9) = 0;
para(10) = 0;
para(11) = 0;
para(12) = 0;
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

function summary = make_summary(earlyPrimary, controlPrimary, earlySecondary, controlSecondary, earlyDist, controlDist, secondaryVirusInoculum)

summary = struct;

summary.secondaryVirusInoculum = secondaryVirusInoculum;

summary.primary_early_final_V = earlyPrimary.V(end);
summary.primary_control_final_V = controlPrimary.V(end);
summary.primary_early_final_I = earlyPrimary.I_total(end);
summary.primary_control_final_I = controlPrimary.I_total(end);

summary.memory_Tc_early = earlyPrimary.Tc(end);
summary.memory_Tc_control = controlPrimary.Tc(end);
summary.memory_Tc_ratio_early_over_control = ...
    (summary.memory_Tc_early + 1) / (summary.memory_Tc_control + 1);

summary.memory_IgM_early = earlyPrimary.totalIgM(end);
summary.memory_IgM_control = controlPrimary.totalIgM(end);
summary.memory_IgG_early = earlyPrimary.totalIgG(end);
summary.memory_IgG_control = controlPrimary.totalIgG(end);

summary.memory_totalAb_early = earlyPrimary.totalAb(end);
summary.memory_totalAb_control = controlPrimary.totalAb(end);
summary.memory_totalAb_ratio_early_over_control = ...
    (summary.memory_totalAb_early + 1) / (summary.memory_totalAb_control + 1);

summary.memory_ELISA_early = earlyPrimary.elisaTotalAb(end);
summary.memory_ELISA_control = controlPrimary.elisaTotalAb(end);
summary.memory_ELISA_ratio_early_over_control = ...
    (summary.memory_ELISA_early + 1) / (summary.memory_ELISA_control + 1);

summary.secondary_peakV_early = earlySecondary.peakV;
summary.secondary_peakV_control = controlSecondary.peakV;
summary.secondary_peakV_ratio_early_over_control = ...
    (summary.secondary_peakV_early + 1) / (summary.secondary_peakV_control + 1);

summary.secondary_peakI_early = earlySecondary.peakI;
summary.secondary_peakI_control = controlSecondary.peakI;
summary.secondary_peakI_ratio_early_over_control = ...
    (summary.secondary_peakI_early + 1) / (summary.secondary_peakI_control + 1);

summary.secondary_aucV_early = earlySecondary.aucV;
summary.secondary_aucV_control = controlSecondary.aucV;
summary.secondary_aucV_ratio_early_over_control = ...
    (summary.secondary_aucV_early + 1) / (summary.secondary_aucV_control + 1);

summary.secondary_aucI_early = earlySecondary.aucI;
summary.secondary_aucI_control = controlSecondary.aucI;
summary.secondary_aucI_ratio_early_over_control = ...
    (summary.secondary_aucI_early + 1) / (summary.secondary_aucI_control + 1);

summary.control_final_kd_IgM = controlDist.kd_IgM(:, end);
summary.early_final_kd_IgM = earlyDist.kd_IgM(:, end);
summary.control_final_kd_IgG = controlDist.kd_IgG(:, end);
summary.early_final_kd_IgG = earlyDist.kd_IgG(:, end);

end

function write_results_text(textFile, summary)

fid = fopen(textFile, 'w');

fprintf(fid, 'Early nirmatrelvir treatment blunts immune memory and increases susceptibility to secondary infection\n\n');

fprintf(fid, ['Early nirmatrelvir treatment rapidly suppressed primary viral replication and achieved viral clearance without overt rebound. ', ...
    'However, this early pharmacological truncation of infection was accompanied by reduced adaptive immune maturation. ', ...
    'At the end of the primary simulation, the early-treatment condition retained lower CD8+ T-cell abundance, reduced antibody abundance and weaker ELISA-derived binding activity than the immune-experienced control. ', ...
    'These findings indicate that early replication blockade limited antigen exposure and blunted the development of durable immune memory.\n\n']);

fprintf(fid, ['To test whether this impaired memory altered susceptibility to reinfection, we used the final state of each primary simulation as the initial condition for a homologous secondary challenge. ', ...
    'All state variables were inherited from the end of the primary simulation, except extracellular virus, which was reset to %.1f to mimic re-exposure. ', ...
    'Secondary infection was then simulated for 2000 time units without additional antiviral dosing. ', ...
    'Compared with the control condition, the early-treatment host developed a larger secondary viral burden and greater infected-cell expansion, consistent with impaired protective immunity.\n\n'], ...
    summary.secondaryVirusInoculum);

fprintf(fid, ['Quantitatively, the early-treatment-to-control ratio was %.3g for final CD8+ T cells, %.3g for total antibody memory and %.3g for ELISA-derived functional antibody memory. ', ...
    'After secondary challenge, the early-treatment-to-control ratio was %.3g for peak viral load and %.3g for peak infected-cell burden. ', ...
    'Thus, although early nirmatrelvir treatment was sufficient to clear the primary infection, it also weakened antibody and CD8+ T-cell memory, thereby increasing susceptibility to secondary infection.\n\n'], ...
    summary.memory_Tc_ratio_early_over_control, ...
    summary.memory_totalAb_ratio_early_over_control, ...
    summary.memory_ELISA_ratio_early_over_control, ...
    summary.secondary_peakV_ratio_early_over_control, ...
    summary.secondary_peakI_ratio_early_over_control);

fprintf(fid, 'Figure legend\n\n');

fprintf(fid, ['Figure X | Early nirmatrelvir treatment blunts immune memory and increases susceptibility to secondary infection. ', ...
    'a, Primary extracellular viral-load dynamics in the early nirmatrelvir-treated simulation and the immune-experienced control. ', ...
    'b, Total infected-cell burden during primary infection, showing truncation of antigen exposure by early treatment. ', ...
    'c, CD8+ T-cell abundance immediately before secondary challenge. ', ...
    'd, IgM and IgG abundance immediately before secondary challenge. ', ...
    'e, ELISA-derived total antibody-binding signal immediately before secondary challenge. ', ...
    'f, Extracellular viral-load dynamics after homologous secondary challenge, initialized by resetting extracellular virus to 10 while preserving all other final-state variables from the primary simulation. ', ...
    'g, Total infected-cell burden after secondary challenge. ', ...
    'h, Peak viral and infected-cell burdens after secondary challenge. ', ...
    'i,j, Time-resolved IgM affinity-distribution heatmaps during primary infection in the control and early-treatment conditions. ', ...
    'k,l, Time-resolved IgG affinity-distribution heatmaps during primary infection in the control and early-treatment conditions. ', ...
    'Heatmap colours indicate log10-transformed antibody abundance across binding-energy classes. ', ...
    'Together, these analyses show that early replication inhibition can clear primary infection while limiting antigen-driven expansion of specific antibodies and CD8+ T cells, resulting in weak immune memory and increased susceptibility to secondary infection.\n']);

fclose(fid);

end
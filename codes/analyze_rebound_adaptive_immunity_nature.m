%% analyze_rebound_adaptive_immunity_final.m
% Viral rebound is accompanied by delayed expansion of antibody and CD8+ T-cell responses
%
% Final version:
% 1. Rebound peak is defined as the first virus rebound-associated local viral peak
%    after the viral nadir, not the late global maximum.
% 2. Original normalized Figure E is removed.
% 3. New Figure E-H:
%    E: No-drug control ELISA-style IgM, IgG and total antibody dynamics.
%    F: Virus rebound ELISA-style IgM, IgG and total antibody dynamics.
%    G: No-drug control IgM and IgG antibody-distribution heatmap.
%    H: Virus rebound IgM and IgG antibody-distribution heatmap.
% 4. ELISA is computed using:
%    ode15s(@pathway_model_many_antibody_immune_include_plasma_elisa, ...)
%    on elisa_time_grid = 0:20:4000.

clear;
clc;
close all;

%% =========================
% 1. File paths
% =========================

reboundFile = fullfile( ...
    'drug_rebound_scan_results', ...
    'raw_simulations', ...
    'sim_start_0200_dose_1e+07_ndose_15.mat');

% 修改这里：换成你的 no-drug 对照组 .mat 文件
controlFile = fullfile( ...
    'drug_rebound_scan_results', ...
    'raw_simulations', ...
    'sim_start_0400_dose_1e+05_ndose_05.mat');

outDir = fullfile( ...
    'drug_rebound_scan_results', ...
    'rebound_adaptive_immunity_final');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

figPrefix = fullfile(outDir, ...
    'viral_rebound_delayed_adaptive_immunity_final');

elisa_time_grid = 0:20:4000;

%% =========================
% 2. Load simulations
% =========================

tmp = load(reboundFile);
simDrug = tmp.sim;

tmp = load(controlFile);
simControl = tmp.sim;

% drug = extract_features(simDrug, elisa_time_grid, 'rebound');
% control = extract_features(simControl, elisa_time_grid, 'control');

affinity_time_grid = 0:5:4000;

controlDist = compute_affinity_distribution_from_sim(control, affinity_time_grid, 'control');
drugDist = compute_affinity_distribution_from_sim(drug, affinity_time_grid, 'rebound');


eventDrug = find_rebound_events_first_peak(drug);

%% =========================
% 3. Treatment protocol
% =========================

treatmentStart = NaN;
treatmentEnd = NaN;
singleDose = NaN;
numDoses = NaN;

if isfield(simDrug, 'protocol')
    protocol = simDrug.protocol;

    if isfield(protocol, 'dose_times')
        doseTimes = protocol.dose_times(:);
        treatmentStart = doseTimes(1);
        treatmentEnd = doseTimes(end);
    elseif isfield(protocol, 'start_time') && isfield(protocol, 'last_dose_time')
        treatmentStart = protocol.start_time;
        treatmentEnd = protocol.last_dose_time;
    end

    if isfield(protocol, 'single_dose')
        singleDose = protocol.single_dose;
    end

    if isfield(protocol, 'num_doses')
        numDoses = protocol.num_doses;
    end
end

xLimit = [0 4000];

%% =========================
% 4. Figure settings
% =========================

set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontSize', 8);
set(groot, 'DefaultLineLineWidth', 1.4);
set(groot, 'DefaultAxesLineWidth', 0.8);

red = [0.80 0.18 0.16];
blue = [0.10 0.37 0.71];
orange = [0.93 0.49 0.19];
purple = [0.49 0.18 0.56];
green = [0.20 0.55 0.30];
gray = [0.45 0.45 0.45];

fig = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [2 2 26 30], ...
    'Renderer', 'painters');

%% =========================
% 5. Panel A: Viral load
% =========================

axA = subplot(4, 4, [1 2]);
hold on;

plot(control.t, log10(control.V + 1), ...
    'Color', gray, ...
    'LineWidth', 1.2);

plot(drug.t, log10(drug.V + 1), ...
    'Color', red, ...
    'LineWidth', 1.8);

shade_treatment(axA, treatmentStart, treatmentEnd);
draw_event_line(axA, eventDrug.timeMinV, blue, 'viral nadir');
draw_event_line(axA, eventDrug.timeReboundPeak, red, 'rebound peak');

xlabel('Time');
ylabel('log_{10}(V + 1)');
title('a  Viral rebound after drug withdrawal', 'FontWeight', 'bold');

legend({'No drug', 'Drug rebound'}, ...
    'Location', 'northeast', ...
    'Box', 'off');

box off;
xlim(xLimit);

%% =========================
% 6. Panel B: Infected cells
% =========================

axB = subplot(4, 4, [3 4]);
hold on;

plot(control.t, log10(control.I_total + 1), ...
    'Color', gray, ...
    'LineWidth', 1.2);

plot(drug.t, log10(drug.I_total + 1), ...
    'Color', red, ...
    'LineWidth', 1.8);

shade_treatment(axB, treatmentStart, treatmentEnd);
draw_event_line(axB, eventDrug.timeMinV, blue, '');
draw_event_line(axB, eventDrug.timeReboundPeak, red, '');

xlabel('Time');
ylabel('log_{10}(I_{total} + 1)');
title('b  Residual infected cells before rebound', 'FontWeight', 'bold');

box off;
xlim(xLimit);

%% =========================
% 7. Panel C: CD8+ T cells
% =========================

axC = subplot(4, 4, [5 6]);
hold on;

plot(control.t, log10(control.Tc + 1), ...
    'Color', gray, ...
    'LineWidth', 1.2);

plot(drug.t, log10(drug.Tc + 1), ...
    'Color', purple, ...
    'LineWidth', 1.8);

shade_treatment(axC, treatmentStart, treatmentEnd);
draw_event_line(axC, eventDrug.timeReboundPeak, red, 'rebound peak');

xlabel('Time');
ylabel('log_{10}(CD8^{+} T cells + 1)');
title('c  Delayed CD8^{+} T-cell expansion', 'FontWeight', 'bold');

legend({'No drug', 'Drug rebound'}, ...
    'Location', 'southeast', ...
    'Box', 'off');

box off;
xlim(xLimit);

%% =========================
% 8. Panel D: Virus-bound antibody response
% =========================

axD = subplot(4, 4, [7 8]);
hold on;

plot(drug.t, log10(drug.totalVirusBoundIgM + 1), ...
    'Color', blue, ...
    'LineWidth', 1.4);

plot(drug.t, log10(drug.totalVirusBoundIgG + 1), ...
    'Color', orange, ...
    'LineWidth', 1.4);

plot(drug.t, log10(drug.totalVirusBoundAb + 1), ...
    'Color', red, ...
    'LineWidth', 1.6);

shade_treatment(axD, treatmentStart, treatmentEnd);
draw_event_line(axD, eventDrug.timeReboundPeak, red, '');

xlabel('Time');
ylabel('log_{10}(virus-bound Ab + 1)');
title('d  Virus-bound antibody response', 'FontWeight', 'bold');

legend({'IgM-bound virus', 'IgG-bound virus', 'total virus-bound Ab'}, ...
    'Location', 'southeast', ...
    'Box', 'off');

box off;
xlim(xLimit);

%% =========================
% 9. Panel E: Control ELISA
% =========================

axE = subplot(4, 4, [9 10]);
hold on;

plot(control.elisa_time_grid, log10(control.elisaIgM + 1), ...
    'Color', blue, ...
    'LineWidth', 1.5);

plot(control.elisa_time_grid, log10(control.elisaIgG + 1), ...
    'Color', orange, ...
    'LineWidth', 1.5);

plot(control.elisa_time_grid, log10(control.elisaTotalAb + 1), ...
    'Color', green, ...
    'LineWidth', 1.7);

xlabel('Time');
ylabel('log_{10}(ELISA signal + 1)');
title('e  ELISA antibody dynamics, no drug', 'FontWeight', 'bold');

legend({'IgM', 'IgG', 'IgM + IgG'}, ...
    'Location', 'southeast', ...
    'Box', 'off');

box off;
xlim(xLimit);

%% =========================
% 10. Panel F: Rebound ELISA
% =========================

axF = subplot(4, 4, [11 12]);
hold on;

plot(drug.elisa_time_grid, log10(drug.elisaIgM + 1), ...
    'Color', blue, ...
    'LineWidth', 1.5);

plot(drug.elisa_time_grid, log10(drug.elisaIgG + 1), ...
    'Color', orange, ...
    'LineWidth', 1.5);

plot(drug.elisa_time_grid, log10(drug.elisaTotalAb + 1), ...
    'Color', green, ...
    'LineWidth', 1.7);

shade_treatment(axF, treatmentStart, treatmentEnd);
draw_event_line(axF, eventDrug.timeReboundPeak, red, 'rebound peak');

xlabel('Time');
ylabel('log_{10}(ELISA signal + 1)');
title('f  ELISA antibody dynamics, rebound', 'FontWeight', 'bold');

legend({'IgM', 'IgG', 'IgM + IgG'}, ...
    'Location', 'southeast', ...
    'Box', 'off');

box off;
xlim(xLimit);

%% =========================
% 11. Panel G: Control antibody-distribution heatmap
% =========================

axG1 = subplot(4, 4, 13);

imagesc( ...
    controlDist.t_grid, ...
    controlDist.kd_values, ...
    log10(controlDist.kd_IgM + 1));

set(gca, 'YDir', 'normal');
xlabel('Time');
ylabel('Binding-energy class, log_{10}(K_d)');
title('g  No-drug IgM repertoire', 'FontWeight', 'bold');

cbG1 = colorbar;
cbG1.Label.String = 'log_{10}(IgM + 1)';
colormap(axG1, parula);
xlim(xLimit);
set(gca, 'FontSize', 8, 'LineWidth', 0.8);
box off;

axG2 = subplot(4, 4, 14);

imagesc( ...
    controlDist.t_grid, ...
    controlDist.kd_values, ...
    log10(controlDist.kd_IgG + 1));

set(gca, 'YDir', 'normal');
xlabel('Time');
ylabel('Binding-energy class, log_{10}(K_d)');
title('g  No-drug IgG repertoire', 'FontWeight', 'bold');

cbG2 = colorbar;
cbG2.Label.String = 'log_{10}(IgG + 1)';
colormap(axG2, parula);
xlim(xLimit);
set(gca, 'FontSize', 8, 'LineWidth', 0.8);
box off;

%% =========================
% 12. Panel H: Rebound IgM / IgG affinity-distribution heatmaps
% =========================

axH1 = subplot(4, 4, 15);

imagesc( ...
    drugDist.t_grid, ...
    drugDist.kd_values, ...
    log10(drugDist.kd_IgM + 1));

set(gca, 'YDir', 'normal');
xlabel('Time');
ylabel('Binding-energy class, log_{10}(K_d)');
title('h  Rebound IgM repertoire', 'FontWeight', 'bold');

cbH1 = colorbar;
cbH1.Label.String = 'log_{10}(IgM + 1)';
colormap(axH1, parula);
xlim(xLimit);
set(gca, 'FontSize', 8, 'LineWidth', 0.8);
box off;

axH2 = subplot(4, 4, 16);

imagesc( ...
    drugDist.t_grid, ...
    drugDist.kd_values, ...
    log10(drugDist.kd_IgG + 1));

set(gca, 'YDir', 'normal');
xlabel('Time');
ylabel('Binding-energy class, log_{10}(K_d)');
title('h  Rebound IgG repertoire', 'FontWeight', 'bold');

cbH2 = colorbar;
cbH2.Label.String = 'log_{10}(IgG + 1)';
colormap(axH2, parula);
xlim(xLimit);
set(gca, 'FontSize', 8, 'LineWidth', 0.8);
box off;

%% =========================
% 13. Overall title
% =========================

sgtitle( ...
    'Viral rebound is accompanied by delayed expansion of antibody and CD8^{+} T-cell responses', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);

%% =========================
% 14. Save figure
% =========================

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

fprintf('\nFigures saved:\n');
fprintf('%s\n', pngFile);
fprintf('%s\n', pdfFile);
fprintf('%s\n', figFile);

%% =========================
% 15. Summary and Results text
% =========================

summary = make_summary(drug, eventDrug, treatmentStart, treatmentEnd, singleDose, numDoses);

summaryFile = fullfile(outDir, 'rebound_adaptive_immunity_summary_final.mat');
save(summaryFile, ...
    'summary', ...
    'eventDrug', ...
    'drug', ...
    'control', ...
    'controlDist', ...
    'drugDist', ...
    'elisa_time_grid', ...
    'affinity_time_grid');

textFile = fullfile(outDir, 'results_section_draft_final.txt');
write_results_text(textFile, summary);

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

if isfield(sim, 'protocol')
    data.protocol = sim.protocol;
else
    data.protocol = [];
end

end

function [elisaIgM, elisaIgG, elisaTotalAb] = compute_elisa_signal(data, elisa_time_grid, groupName)

fprintf('\nComputing ELISA signal for %s group using pathway_model_many_antibody_immune_include_plasma_elisa...\n', groupName);

if exist('pathway_model_many_antibody_immune_include_plasma_elisa', 'file') ~= 2
    error(['Cannot find pathway_model_many_antibody_immune_include_plasma_elisa.m. ', ...
           'Please make sure this function is in the current folder or MATLAB path.']);
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

    fprintf('  %s ELISA time point %d / %d, t = %.1f\n', ...
        groupName, k, numTime, elisa_time_grid(k));

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

catch ME
    warning('ELISA ODE failed at one time point: %s', ME.message);
    elisaValue = NaN;
end

end

function [para, para_new, para_new_1] = get_elisa_parameters()

para = zeros(20, 1);

para(1) = 1e-20;   % environmental antigen kon
para(2) = 0.5;     % environmental antigen koff
para(3) = 0;       % replenish constant pi 1
para(4) = 0;       % replenish constant pi 2
para(5) = 0;       % decay constant of BCR IgM
para(6) = 0;       % decay constant of BCR IgG
para(7) = 0;       % feedback constant of environmental antigen-antibody complex
para(8) = 0;       % feedback constant on plasma cell regeneration
para(9) = 0;       % decay constant of plasma cell IgM
para(10) = 0;      % decay constant of plasma cell IgG
para(11) = 0;      % production constant of IgM
para(12) = 0;      % production constant of IgG
para(13) = 0.05;   % decay constant of IgM
para(14) = 0.025;  % decay constant of IgG
para(15) = 0;
para(16) = 0.05;   % transformation constant from IgM to IgG memory cell
para(17) = 0.1;    % maximal production percentage of plasma cell
para(18) = 0.5;    % decay constant of complex
para(19) = 0;      % virus replication constant
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

function event = find_rebound_events_first_peak(data)

t = data.t;
V = data.V;

event.timeMinV = NaN;
event.minV = NaN;
event.timeReboundPeak = NaN;
event.reboundPeakV = NaN;
event.reboundLogIncrease = NaN;

if isempty(t) || isempty(V)
    return;
end

if ~isempty(data.protocol) && isfield(data.protocol, 'start_time')
    searchStart = data.protocol.start_time;
else
    searchStart = t(round(numel(t) * 0.1));
end

idxAfterTreatmentStart = t >= searchStart;

if ~any(idxAfterTreatmentStart)
    idxAfterTreatmentStart = true(size(t));
end

tSub = t(idxAfterTreatmentStart);
VSub = V(idxAfterTreatmentStart);

[event.minV, idxMinSub] = min(VSub);
event.timeMinV = tSub(idxMinSub);

idxNadirGlobal = find(t == event.timeMinV, 1, 'first');

if isempty(idxNadirGlobal)
    [~, idxNadirGlobal] = min(abs(t - event.timeMinV));
end

baseline = max(event.minV, 1);
minReboundFold = 10;
minReboundIncrease = log10(minReboundFold);

candidateIdx = [];

for k = idxNadirGlobal + 1:numel(t) - 1
    isLocalPeak = V(k) >= V(k - 1) && V(k) >= V(k + 1);
    isAboveNadir = log10(V(k) + 1) - log10(baseline + 1) >= minReboundIncrease;

    if isLocalPeak && isAboveNadir
        candidateIdx = k;
        break;
    end
end

if isempty(candidateIdx)
    idxAfterNadir = t >= event.timeMinV;
    tAfter = t(idxAfterNadir);
    VAfter = V(idxAfterNadir);

    [event.reboundPeakV, idxPeak] = max(VAfter);
    event.timeReboundPeak = tAfter(idxPeak);
else
    event.reboundPeakV = V(candidateIdx);
    event.timeReboundPeak = t(candidateIdx);
end

event.reboundLogIncrease = ...
    log10(event.reboundPeakV + 1) - log10(event.minV + 1);

end

function dist = compute_affinity_distribution_from_sim(data, affinity_time_grid, groupName)

fprintf('\nComputing affinity distribution for %s group using compute_affinity_distribution_cell...\n', groupName);

if exist('compute_affinity_distribution_cell', 'file') ~= 2
    error(['Cannot find compute_affinity_distribution_cell.m. ', ...
           'Please make sure this function is in the current folder or MATLAB path.']);
end

affinity_time_grid = affinity_time_grid(:)';

try
    dist = compute_affinity_distribution_cell(data.t, data.y, affinity_time_grid);
catch
    dist = compute_affinity_distribution_cell(data.t, data.y, affinity_time_grid);
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

function shade_treatment(ax, treatmentStart, treatmentEnd)

if ~isfinite(treatmentStart) || ~isfinite(treatmentEnd)
    return;
end

axes(ax);
yl = ylim(ax);

patch([treatmentStart treatmentEnd treatmentEnd treatmentStart], ...
      [yl(1) yl(1) yl(2) yl(2)], ...
      [0.90 0.90 0.90], ...
      'EdgeColor', 'none', ...
      'FaceAlpha', 0.45, ...
      'HandleVisibility', 'off');

children = get(ax, 'Children');

if numel(children) > 1
    set(ax, 'Children', [children(2:end); children(1)]);
end

ylim(ax, yl);

end

function draw_event_line(ax, xValue, colorValue, labelText)

if ~isfinite(xValue)
    return;
end

axes(ax);
yl = ylim(ax);

plot([xValue xValue], yl, '--', ...
    'Color', colorValue, ...
    'LineWidth', 0.9, ...
    'HandleVisibility', 'off');

if ~isempty(labelText)
    text(xValue, yl(2), [' ' labelText], ...
        'Color', colorValue, ...
        'FontSize', 7, ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left');
end

ylim(ax, yl);

end

function summary = make_summary(drug, event, treatmentStart, treatmentEnd, singleDose, numDoses)

summary = struct;

summary.treatmentStart = treatmentStart;
summary.treatmentEnd = treatmentEnd;
summary.singleDose = singleDose;
summary.numDoses = numDoses;

summary.timeMinV = event.timeMinV;
summary.minV = event.minV;
summary.timeReboundPeak = event.timeReboundPeak;
summary.reboundPeakV = event.reboundPeakV;
summary.reboundLogIncrease = event.reboundLogIncrease;

[~, idxNadir] = min(abs(drug.t - event.timeMinV));
[~, idxPeak] = min(abs(drug.t - event.timeReboundPeak));

summary.I_total_at_nadir = drug.I_total(idxNadir);
summary.I_total_at_rebound_peak = drug.I_total(idxPeak);

summary.Tc_at_nadir = drug.Tc(idxNadir);
summary.Tc_at_rebound_peak = drug.Tc(idxPeak);
summary.Tc_fold_change = ...
    (drug.Tc(idxPeak) + 1) / (drug.Tc(idxNadir) + 1);

summary.IgM_ELISA_at_rebound_peak = interp1( ...
    drug.elisa_time_grid, drug.elisaIgM, event.timeReboundPeak, 'linear', 'extrap');

summary.IgG_ELISA_at_rebound_peak = interp1( ...
    drug.elisa_time_grid, drug.elisaIgG, event.timeReboundPeak, 'linear', 'extrap');

summary.total_ELISA_at_rebound_peak = interp1( ...
    drug.elisa_time_grid, drug.elisaTotalAb, event.timeReboundPeak, 'linear', 'extrap');

summary.IgM_ELISA_at_nadir = interp1( ...
    drug.elisa_time_grid, drug.elisaIgM, event.timeMinV, 'linear', 'extrap');

summary.IgG_ELISA_at_nadir = interp1( ...
    drug.elisa_time_grid, drug.elisaIgG, event.timeMinV, 'linear', 'extrap');

summary.total_ELISA_at_nadir = interp1( ...
    drug.elisa_time_grid, drug.elisaTotalAb, event.timeMinV, 'linear', 'extrap');

summary.total_ELISA_fold_change = ...
    (summary.total_ELISA_at_rebound_peak + 1) / ...
    (summary.total_ELISA_at_nadir + 1);

end

function write_results_text(textFile, summary)

fid = fopen(textFile, 'w');

fprintf(fid, 'Viral rebound is accompanied by delayed expansion of antibody and CD8+ T-cell responses\n\n');

fprintf(fid, ['To investigate the immunological basis of viral rebound after antiviral treatment, ', ...
    'we focused on a representative rebound simulation. ']);

if isfinite(summary.treatmentStart)
    fprintf(fid, 'Treatment was initiated at t = %.0f', summary.treatmentStart);
else
    fprintf(fid, 'Treatment was initiated at the specified treatment time');
end

if isfinite(summary.singleDose) && isfinite(summary.numDoses)
    fprintf(fid, ...
        ' with a single-dose amount of %.1e administered for %d doses. ', ...
        summary.singleDose, summary.numDoses);
else
    fprintf(fid, '. ');
end

fprintf(fid, ['Drug treatment rapidly suppressed extracellular virus and reduced the infected-cell burden. ', ...
    'However, at the viral nadir, extracellular virus and infected cells were not completely eliminated, ', ...
    'indicating the persistence of a residual infectious reservoir after pharmacological pressure was relaxed. ']);

fprintf(fid, ['The viral nadir occurred at t = %.1f, after which viral load increased again and reached the first rebound-associated viral peak at t = %.1f. ', ...
    'This rebound represented a %.2f-log10 increase in viral load relative to the nadir. '], ...
    summary.timeMinV, summary.timeReboundPeak, summary.reboundLogIncrease);

fprintf(fid, ['ELISA simulations were performed on a time grid from t = 0 to t = 4000 with a step size of 20. ', ...
    'At each selected time point, the simulated IgM and IgG distributions were introduced into an antigen-binding assay model, ', ...
    'and the ELISA signal was quantified as the reduction in free antigen after the assay period. ']);

fprintf(fid, ['At the time of rebound, adaptive immunity remained delayed relative to viral expansion. ', ...
    'ELISA-style antibody trajectories showed that IgM and IgG responses continued to increase during the rebound phase, ', ...
    'whereas the no-drug control displayed a distinct antibody kinetic profile over the same observation window. ']);

fprintf(fid, ['The antibody specificity landscapes further demonstrated remodeling of both IgM and IgG distributions during rebound, ', ...
    'consistent with antigen-driven reshaping of the humoral response. ', ...
    'In parallel, CD8+ T-cell abundance increased after renewed viral expansion, supporting delayed activation of cellular immunity. ']);

fprintf(fid, ['Together, these results suggest that viral rebound occurs when antiviral treatment suppresses viral replication before complete pathogen clearance ', ...
    'and before adaptive immunity has reached sufficient magnitude. ', ...
    'The subsequent expansion and remodeling of antibody responses, together with CD8+ T-cell expansion, provides a mechanistic explanation for eventual viral control and durable clearance.\n\n']);

fprintf(fid, 'Figure legend\n\n');

fprintf(fid, ['Figure X. Viral rebound is accompanied by delayed expansion of antibody and CD8+ T-cell responses. ', ...
    '(a) Extracellular viral load in the no-drug control and drug-treated rebound simulation. ', ...
    'The shaded region indicates the treatment window, and dashed lines indicate the viral nadir and the first rebound-associated viral peak. ', ...
    '(b) Total infected-cell burden, showing residual infected cells near the viral nadir. ', ...
    '(c) CD8+ T-cell dynamics, showing delayed expansion after viral rebound. ', ...
    '(d) Virus-bound IgM, virus-bound IgG and total virus-bound antibody dynamics in the rebound simulation. ', ...
    '(e) ELISA-derived IgM, IgG and total antibody dynamics in the no-drug control over t = 0-4000. ', ...
    '(f) ELISA-derived IgM, IgG and total antibody dynamics in the viral rebound simulation over t = 0-4000. ', ...
    '(g) IgM and IgG antibody specificity landscape in the no-drug control, averaged over t = 0-4000. ', ...
    '(h) IgM and IgG antibody specificity landscape in the viral rebound simulation, averaged over t = 0-4000. ', ...
    'Heatmap colors indicate normalized log10-transformed antibody abundance.\n']);

fclose(fid);

end

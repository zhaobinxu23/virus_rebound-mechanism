function plot_drug_rebound_parameter_scan_nature(Result, Representative, figDir)
% plot_drug_rebound_parameter_scan_nature
% Nature-style multi-panel figure for drug rebound parameter scan.

if ~exist(figDir, 'dir')
    mkdir(figDir);
end

%% =========================
% 1. Style
% ==========================
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 8);
set(groot, 'defaultTextFontSize', 8);
set(groot, 'defaultAxesLineWidth', 0.75);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');
set(groot, 'defaultLegendBox', 'off');

startTimes = unique(Result.start_time);
doseNumbers = unique(Result.num_doses);
singleDoses = unique(Result.single_dose);

%% =========================
% 2. Build heatmap matrices
% ==========================
% We show heatmaps for each start time:
% x-axis: number of doses
% y-axis: single dose
% color: log10 rebound peak

Heat = cell(numel(startTimes), 1);

for i = 1:numel(startTimes)
    st = startTimes(i);
    H = nan(numel(singleDoses), numel(doseNumbers));

    for j = 1:numel(singleDoses)
        for k = 1:numel(doseNumbers)
            idx = Result.start_time == st & ...
                  Result.single_dose == singleDoses(j) & ...
                  Result.num_doses == doseNumbers(k);

            if any(idx)
                H(j, k) = Result.log10_rebound_peak_V(idx);
            end
        end
    end

    Heat{i} = H;
end

allHeatValues = [];
for i = 1:numel(Heat)
    allHeatValues = [allHeatValues; Heat{i}(:)];
end
clim = [min(allHeatValues, [], 'omitnan'), max(allHeatValues, [], 'omitnan')];

%% =========================
% 3. Figure
% ==========================
fig = figure('Units', 'centimeters', ...
    'Position', [2, 2, 19.0, 18.0], ...
    'Color', 'w');

tl = tiledlayout(fig, 3, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% =========================
% Panel A: representative viral dynamics
% ==========================
axA = nexttile(tl, 1, [1 2]);
hold(axA, 'on');

plot_representative_curve(axA, Representative.lowDoseEarlyShort, ...
    [0.20 0.35 0.80], 'Low dose, early, short');

plot_representative_curve(axA, Representative.highDoseEarlyLong, ...
    [0.75 0.10 0.15], 'High dose, early, long');

plot_representative_curve(axA, Representative.lowDoseLateShort, ...
    [0.20 0.60 0.35], 'Low dose, late, short');

plot_representative_curve(axA, Representative.highDoseLateLong, ...
    [0.10 0.10 0.10], 'High dose, late, long');

xlabel(axA, 'Time');
ylabel(axA, 'Viral load, log_{10}(V + 1)');
title(axA, 'Representative viral trajectories');
legend(axA, 'Location', 'northeast', 'FontSize', 6);
format_axis(axA);
panel_label(axA, 'a');

%% =========================
% Panel B: global relationship, dose number vs rebound
% ==========================
axB = nexttile(tl, 3);
hold(axB, 'on');

cmapDose = nature_lines(numel(singleDoses));

for i = 1:numel(singleDoses)
    dose = singleDoses(i);

    yMean = nan(size(doseNumbers));
    ySem = nan(size(doseNumbers));

    for k = 1:numel(doseNumbers)
        idx = Result.single_dose == dose & ...
              Result.num_doses == doseNumbers(k);

        vals = Result.log10_rebound_peak_V(idx);
        yMean(k) = mean(vals, 'omitnan');
        ySem(k) = std(vals, 'omitnan') ./ sqrt(sum(~isnan(vals)));
    end

    shaded_error(axB, doseNumbers, yMean, ySem, cmapDose(i, :));
    plot(axB, doseNumbers, yMean, '-', ...
        'Color', cmapDose(i, :), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('10^{%d}', round(log10(dose))));
end

xlabel(axB, 'Number of doses');
ylabel(axB, 'Rebound peak, log_{10}(V + 1)');
title(axB, 'Dose-length response');
legend(axB, 'Location', 'northeast', 'FontSize', 6);
format_axis(axB);
panel_label(axB, 'b');

%% =========================
% Panels C-G: heatmaps for each start time
% ==========================
for i = 1:numel(startTimes)
    ax = nexttile(tl, i + 3);

    imagesc(ax, doseNumbers, log10(singleDoses), Heat{i});
    set(ax, 'YDir', 'normal');

    colormap(ax, nature_colormap());
    caxis(ax, clim);

    xlabel(ax, 'Number of doses');
    ylabel(ax, 'Dose, log_{10}');
    title(ax, sprintf('Start time = %g', startTimes(i)));

    xticks(ax, [1 5 10 15 20]);
    yticks(ax, log10(singleDoses));
    yticklabels(ax, compose('10^{%d}', round(log10(singleDoses))));

    format_axis(ax);
    panel_label(ax, char('c' + i - 1));

    if i == numel(startTimes)
        cb = colorbar(ax);
        cb.TickDirection = 'out';
        cb.Box = 'off';
        cb.Label.String = 'Rebound peak, log_{10}(V + 1)';
        cb.Label.FontSize = 7;
    end
end

%% =========================
% Panel H: rebound probability
% ==========================
axH = nexttile(tl, 9);
hold(axH, 'on');

for i = 1:numel(startTimes)
    st = startTimes(i);

    reboundProb = nan(size(doseNumbers));

    for k = 1:numel(doseNumbers)
        idx = Result.start_time == st & ...
              Result.num_doses == doseNumbers(k);

        reboundProb(k) = mean(Result.rebound_flag(idx), 'omitnan');
    end

    plot(axH, doseNumbers, reboundProb, '-o', ...
        'LineWidth', 1.3, ...
        'MarkerSize', 3.5, ...
        'DisplayName', sprintf('Start %g', st));
end

xlabel(axH, 'Number of doses');
ylabel(axH, 'Rebound probability');
ylim(axH, [-0.02 1.02]);
title(axH, 'Rebound suppression by regimen length');
legend(axH, 'Location', 'northeast', 'FontSize', 6);
format_axis(axH);
panel_label(axH, 'h');

%% =========================
% Export
% ==========================
pdfName = fullfile(figDir, 'Figure_drug_rebound_parameter_scan_nature.pdf');
pngName = fullfile(figDir, 'Figure_drug_rebound_parameter_scan_nature.png');
svgName = fullfile(figDir, 'Figure_drug_rebound_parameter_scan_nature.svg');

exportgraphics(fig, pdfName, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

exportgraphics(fig, svgName, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

exportgraphics(fig, pngName, ...
    'Resolution', 600, ...
    'BackgroundColor', 'white');

fprintf('Nature-style figure saved:\n%s\n%s\n%s\n', ...
    pdfName, svgName, pngName);
end

%% ========================================================================
% Helper functions
% ========================================================================

function plot_representative_curve(ax, sim, color, labelText)
if isempty(sim) || ~isfield(sim, 't') || isempty(sim.t)
    return;
end

plot(ax, sim.t, log10(sim.V + 1), ...
    'Color', color, ...
    'LineWidth', 1.5, ...
    'DisplayName', labelText);
end

function shaded_error(ax, x, y, e, color)
x = x(:)';
y = y(:)';
e = e(:)';

valid = isfinite(x) & isfinite(y) & isfinite(e);
x = x(valid);
y = y(valid);
e = e(valid);

if isempty(x)
    return;
end

fill(ax, [x fliplr(x)], [y - e fliplr(y + e)], color, ...
    'FaceAlpha', 0.12, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
end

function format_axis(ax)
set(ax, ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'LineWidth', 0.75, ...
    'FontSize', 8, ...
    'Layer', 'top');

ax.Title.FontWeight = 'normal';
ax.XColor = [0.1 0.1 0.1];
ax.YColor = [0.1 0.1 0.1];
grid(ax, 'off');
end

function panel_label(ax, labelText)
text(ax, -0.18, 1.12, labelText, ...
    'Units', 'normalized', ...
    'FontWeight', 'bold', ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');
end

function cmap = nature_colormap()
n = 256;

c1 = [0.231 0.298 0.753];
c2 = [1.000 1.000 1.000];
c3 = [0.706 0.016 0.150];

x = linspace(0, 1, n)';
cmap = zeros(n, 3);

for i = 1:n
    if x(i) <= 0.5
        a = x(i) / 0.5;
        cmap(i, :) = (1 - a) * c1 + a * c2;
    else
        a = (x(i) - 0.5) / 0.5;
        cmap(i, :) = (1 - a) * c2 + a * c3;
    end
end
end

function cols = nature_lines(n)
base = [
    0.230 0.299 0.754
    0.706 0.016 0.150
    0.000 0.533 0.337
    0.865 0.407 0.094
    0.494 0.184 0.556
    0.100 0.100 0.100
];

if n <= size(base, 1)
    cols = base(1:n, :);
else
    x = linspace(1, size(base, 1), n);
    cols = interp1(1:size(base, 1), base, x);
end
end
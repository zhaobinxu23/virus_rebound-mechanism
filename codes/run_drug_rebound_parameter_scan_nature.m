%% run_drug_rebound_parameter_scan_nature.m
% Parameter scan for drug dose, treatment start time, and total number of doses.
% Total simulations: 5 x 5 x 20 = 500.
%
% Required existing files:
%   main_consider_drug_dose.m
%   sys_ode_modified_new.m
%   compute_vin_remdesivir_total_fast.m
%   drug_concentration_fast.m

%% run_drug_rebound_parameter_scan_nature.m
clear;
clc;

%% =========================
% 1. Load base model
% ==========================
% 注意：main_consider_drug_dose.m 里面如果有 clear，
% 会清除本脚本前面定义的变量，所以必须先运行它。
run('main_consider_drug_dose.m');

%% =========================
% 2. Re-create output folders after running main script
% ==========================
outDir = 'drug_rebound_scan_results';
figDir = fullfile(outDir, 'figures');
dataDir = fullfile(outDir, 'raw_simulations');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

if ~exist(figDir, 'dir')
    mkdir(figDir);
end

if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

%% =========================
% 3. Initial condition
% ==========================
if ~exist('y0', 'var')
    y0 = [T0; I0; V0; Tc0; x0; cum_death_0];
end

%% =========================
% 4. Simulation time span
% ==========================
if exist('t_span', 'var')
    tspan = t_span;
elseif ~exist('tspan', 'var')
    tspan = [0 4000];
end

%% =========================
% 5. ODE options
% ==========================
if exist('options', 'var')
    ode_options = options;
elseif ~exist('ode_options', 'var')
    ode_options = odeset( ...
        'RelTol', 1e-4, ...
        'AbsTol', 1e-8, ...
        'NonNegative', 1:numel(y0));
end

%% =========================
% 6. Parameter grid
% ==========================
startTimes = [0 100 200 300 400];
singleDoses = [1e5 1e6 1e7 1e8];
doseNumbers = 5:5:20;
doseInterval = 100;

nStart = numel(startTimes);
nDose = numel(singleDoses);
nNum = numel(doseNumbers);
nTotal = nStart * nDose * nNum;

fprintf('Starting drug rebound parameter scan: %d simulations\n', nTotal);

%% =========================
% 3. Preallocate result table
% ==========================
Result = table;
row = 0;

% Store representative curves for figure
Representative = struct;
Representative.lowDoseEarlyShort = [];
Representative.highDoseEarlyLong = [];
Representative.lowDoseLateShort = [];
Representative.highDoseLateLong = [];

%% =========================
% 4. Run 500 simulations
% ==========================
tic;

for iStart = 1:nStart
    for iDose = 1:nDose
        for iNum = 1:nNum

            row = row + 1;

            protocol.start_time = startTimes(iStart);
            protocol.single_dose = singleDoses(iDose);
            protocol.num_doses = doseNumbers(iNum);
            protocol.dose_interval = doseInterval;
            protocol.dose_times = protocol.start_time + ...
                (0:protocol.num_doses - 1) * protocol.dose_interval;
            protocol.dose_amounts = protocol.single_dose * ...
                ones(size(protocol.dose_times));
            protocol.last_dose_time = protocol.dose_times(end);
            protocol.treatment_duration = ...
                (protocol.num_doses - 1) * protocol.dose_interval;

            fprintf('[%03d/%03d] start=%g, dose=%.1e, n=%d\n', ...
                row, nTotal, protocol.start_time, ...
                protocol.single_dose, protocol.num_doses);

            oneRunTic = tic;
            sim = simulate_one_drug_protocol(p, y0, tspan, ode_options, protocol);
            fprintf('    finished in %.2f seconds, success=%d\n', toc(oneRunTic), sim.success);
            % Save raw simulation
            rawFile = fullfile(dataDir, sprintf( ...
                'sim_start_%04g_dose_%1.0e_ndose_%02d.mat', ...
                protocol.start_time, protocol.single_dose, protocol.num_doses));

            save(rawFile, 'sim', 'protocol', '-v7.3');

            % Add row to result table
            Result.row(row, 1) = row;
            Result.start_time(row, 1) = protocol.start_time;
            Result.single_dose(row, 1) = protocol.single_dose;
            Result.log10_single_dose(row, 1) = log10(protocol.single_dose);
            Result.num_doses(row, 1) = protocol.num_doses;
            Result.dose_interval(row, 1) = protocol.dose_interval;
            Result.last_dose_time(row, 1) = protocol.last_dose_time;
            Result.treatment_duration(row, 1) = protocol.treatment_duration;

            Result.peak_V_all(row, 1) = sim.metrics.peak_V_all;
            Result.time_peak_V_all(row, 1) = sim.metrics.time_peak_V_all;
            Result.min_V_after_treatment(row, 1) = sim.metrics.min_V_after_treatment;
            Result.time_min_V_after_treatment(row, 1) = sim.metrics.time_min_V_after_treatment;
            Result.rebound_peak_V(row, 1) = sim.metrics.rebound_peak_V;
            Result.time_rebound_peak(row, 1) = sim.metrics.time_rebound_peak;
            Result.rebound_amplitude(row, 1) = sim.metrics.rebound_amplitude;
            Result.log10_rebound_peak_V(row, 1) = log10(sim.metrics.rebound_peak_V + 1);
            Result.log10_rebound_amplitude(row, 1) = log10(sim.metrics.rebound_amplitude + 1);
            Result.rebound_flag(row, 1) = sim.metrics.rebound_flag;
            Result.auc_V(row, 1) = sim.metrics.auc_V;
            Result.auc_log10V(row, 1) = sim.metrics.auc_log10V;
            Result.final_V(row, 1) = sim.metrics.final_V;
            Result.raw_file{row, 1} = rawFile;

            % Representative examples for plotting
            if protocol.start_time == 0 && protocol.single_dose == 1e3 && protocol.num_doses == 1
                Representative.lowDoseEarlyShort = sim;
            end

            if protocol.start_time == 0 && protocol.single_dose == 1e7 && protocol.num_doses == 20
                Representative.highDoseEarlyLong = sim;
            end

            if protocol.start_time == 400 && protocol.single_dose == 1e3 && protocol.num_doses == 1
                Representative.lowDoseLateShort = sim;
            end

            if protocol.start_time == 400 && protocol.single_dose == 1e7 && protocol.num_doses == 20
                Representative.highDoseLateLong = sim;
            end
        end
    end
end

elapsedTime = toc;

fprintf('Finished %d simulations in %.2f minutes.\n', ...
    nTotal, elapsedTime / 60);

%% =========================
% 5. Save summary data
% ==========================
summaryMatFile = fullfile(outDir, 'drug_rebound_parameter_scan_summary.mat');
summaryCsvFile = fullfile(outDir, 'drug_rebound_parameter_scan_summary.csv');

save(summaryMatFile, ...
    'Result', 'Representative', ...
    'startTimes', 'singleDoses', 'doseNumbers', ...
    'doseInterval', 'tspan', '-v7.3');

writetable(Result, summaryCsvFile);

fprintf('Summary saved:\n%s\n%s\n', summaryMatFile, summaryCsvFile);

%% =========================
% 6. Plot Nature-style figure
% ==========================
plot_drug_rebound_parameter_scan_nature(Result, Representative, figDir);

fprintf('All done.\n');
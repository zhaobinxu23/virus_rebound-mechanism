%% run_figure2_targetcell_immunity.m
% Figure 2: separation of target-cell depletion and adaptive immune priming
% Models:
%   1 WT
%   2 No T depletion (dT/dt = 0, adaptive immunity intact)
%   3 No adaptive immunity, saturating infection beta*T*V/(V+Km)
%   4 No adaptive immunity, mass-action infection beta*T*V

clear; clc; close all;
run('main_consider_drug_dose.m');   % loads p, y0, options

% Consistent parameters for all four models (same as model_1/2/3)
p.k4 = 1e-3;
p.km = 2e6;
p.para(15) = 5;
p.Tc_generation = 5e-5;

tspan1 = [0 4000];
tspan2 = [0 2000];
options = odeset('RelTol',1e-4,'AbsTol',1e-8, ...
                 'NonNegative',1:numel(y0));

models = struct('name',{},'ode',{},'immune',{});
models(1).name = 'WT';
models(1).ode  = @sys_ode_modified_new;
models(1).immune = true;

models(2).name = 'No T depletion';
models(2).ode  = @sys_ode_modified_model_1;
models(2).immune = true;

models(3).name = 'No adaptive (saturating)';
models(3).ode  = @sys_ode_modified_model_2;
models(3).immune = false;

models(4).name = 'No adaptive (mass action)';
models(4).ode  = @sys_ode_modified_model_3;
models(4).immune = false;

Res = struct();
for m = 1:4
    fprintf('Running model %d/%d: %s\n', m, 4, models(m).name);

    p_m = p;
    if models(m).immune
        p_m.para(15) = 5;
        p_m.Tc_generation = 5e-5;
    else
        % Disable adaptive immunity
        p_m.para(15) = 0;
        p_m.Tc_generation = 0;
        p_m.k_adcc_max = 0;
        p_m.k_tc_max = 0;
        p_m.para_new(:) = 0;
        p_m.para_new_1(:) = 0;
        p_m.para(1) = 0;
        p_m.para(2) = 0;
        p_m.para(7) = 0;
        p_m.para(8) = 0;
    end

    % ---------- Primary infection ----------
    [t1, y1] = ode15s(@(t,y) models(m).ode(t,y,p_m), tspan1, y0, options);

    N = p_m.N;
    T1 = max(y1(:,1),0);
    I1 = max(sum(y1(:,2:N+1),2),0);
    V1 = max(y1(:,N+2),0);
    Tc1 = max(y1(:,N+3),0);

    % ---------- Secondary challenge ----------
    % Equalise target cells, clear residual infection, reset V = 10,
    % preserve antibody and CD8+ T-cell memory variables.
    y_sec = y1(end,:)';
    y_sec(1) = y0(1);          % T = T0 for all models
    y_sec(2:N+1) = 0;          % clear infected cells
    y_sec(N+2) = 10;           % re-exposure dose

    [t2, y2] = ode15s(@(t,y) models(m).ode(t,y,p_m), tspan2, y_sec, options);

    T2 = max(y2(:,1),0);
    I2 = max(sum(y2(:,2:N+1),2),0);
    V2 = max(y2(:,N+2),0);
    Tc2 = max(y2(:,N+3),0);

    Res(m).name = models(m).name;
    Res(m).t1 = t1; Res(m).T1 = T1; Res(m).I1 = I1; Res(m).V1 = V1; Res(m).Tc1 = Tc1;
    Res(m).t2 = t2; Res(m).T2 = T2; Res(m).I2 = I2; Res(m).V2 = V2; Res(m).Tc2 = Tc2;
    Res(m).primary_peak_I = max(I1);
    Res(m).primary_peak_V = max(V1);
    Res(m).primary_final_V = V1(end);
    Res(m).secondary_peak_I = max(I2);
    Res(m).secondary_peak_V = max(V2);
end

%% ---------- Figure 2 ----------
figure('Color','w','Units','centimeters','Position',[2 2 22 16]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
colors = lines(4);
styles = {'-','--','-.',':'};

nexttile;
for m = 1:4
    semilogy(Res(m).t1, max(Res(m).I1,1e-6), styles{m}, ...
        'Color',colors(m,:), 'LineWidth',1.3); hold on;
end
xlabel('Time'); ylabel('Total infected cells');
title('a  Primary infection: infected cells');
legend({Res.name}, 'Location','best', 'Box','off'); box off;

nexttile;
for m = 1:4
    semilogy(Res(m).t1, max(Res(m).V1,1e-6), styles{m}, ...
        'Color',colors(m,:), 'LineWidth',1.3); hold on;
end
xlabel('Time'); ylabel('Virus V');
title('b  Primary infection: virus'); box off;

nexttile;
for m = 1:4
    semilogy(Res(m).t2, max(Res(m).V2,1e-6), styles{m}, ...
        'Color',colors(m,:), 'LineWidth',1.3); hold on;
end
xlabel('Time after rechallenge'); ylabel('Virus V');
title('c  Secondary challenge: virus'); box off;

nexttile;
for m = 1:4
    semilogy(Res(m).t2, max(Res(m).I2,1e-6), styles{m}, ...
        'Color',colors(m,:), 'LineWidth',1.3); hold on;
end
xlabel('Time after rechallenge'); ylabel('Total infected cells');
title('d  Secondary challenge: infected cells'); box off;

exportgraphics(gcf,'Figure2_targetcell_vs_immunity.png','Resolution',300);

%% ---------- Summary table ----------
T = table( ...
    {Res.name}', ...
    [Res.primary_peak_I]', [Res.primary_peak_V]', [Res.primary_final_V]', ...
    [Res.secondary_peak_I]', [Res.secondary_peak_V]', ...
    'VariableNames', {'Model','PrimaryPeakI','PrimaryPeakV', ...
    'PrimaryFinalV','SecondaryPeakI','SecondaryPeakV'});
writetable(T,'Figure2_summary.csv');
disp(T);
save('Figure2_results.mat','Res');
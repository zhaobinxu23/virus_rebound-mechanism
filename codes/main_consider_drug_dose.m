% Description: 求解基于感染龄(Infection Age)结构的病毒感染模型
% Solver: ode15s (适合刚性问题)
% =========================================================================

clear; clc; close all;

%% 1. 参数设置  first infection
p.N  = 400;      % 感染龄的分段数量 (bins)
p.da = 1.0;      % 感染龄步长 (Delta a, 单位通常是小时或步长单位)
p.dt_step = 1;   % 仅仅用于定义物理上的流速，通常为1(时间流逝速率)

%%
p.para(1) = 1e-20; % environmental antigen kon
p.para(2) = 0.5;% environmental antigen koff
p.para(3) = 2.2e17+1e13; % replenish constant pi 1 
p.para(4) = 0.5e13;% replenish constant pi 2
p.para(5) = 0.01; % decay constant of BCR IgM
p.para(6) = 0.005;% decay constant of BCR IgG
p.para(7) = 1; % k2  feedback constant of enviromental antigen-antibody complex
p.para(8) = 5e-7;% k2' feedback constant on PC cell regeneration
p.para(9) = 0.1; % decay constant of plasma Cell IgM
p.para(10) = 0.05;%% decay constant of plasma Cell IgG
p.para(11) = 4.4e9;% production constant of IgM
p.para(12) = 1.2e10;% production constant of IgG

p.para(13) = 0.05;% decay constant of IgM
p.para(14) = 0.025;% decay constant of IgG
p.para(15) = 5;% amplification constant of virus antigen
p.para(16) = 0.02;% transformation constant from IgM to IgG memory cell
p.para(17) = 0.1;% maximal production percentage of plasma cell 
p.para(18) = 0.5;% decay constant of complex
p.para(19) = 0; % virus replication constant
p.para(20) = 1e5;



p.para_new(1) = 1e-22; 
p.para_new(2) = 1e-21;
p.para_new(3) = 1e-20; 
p.para_new(4) = 1e-19;
p.para_new(5) = 1e-18;
p.para_new(6) = 1e-17;
p.para_new(7) = 1e-16; 
p.para_new(8) = 1e-15;
p.para_new(9) = 1e-14; 
p.para_new(10) = 1e-13;

p.para_new_1(1) = 1e0; 
p.para_new_1(2) = 1e1;
p.para_new_1(3) = 1e2; 
p.para_new_1(4) = 1e3;
p.para_new_1(5) = 1e4; 
p.para_new_1(6) = 1e5;
p.para_new_1(7) = 1e6; 
p.para_new_1(8) = 1e7;
p.para_new_1(9) = 1e8; 
p.para_new_1(10) = 1e9;

% --- 病毒与易感细胞参数 ---
p.k4 = 1.0e-3;   % 病毒入侵常数
p.km = 2.0e6;    % 半饱和常数
p.k6 = 1.0e8;       % 易感细胞细胞再生率
p.k7 = 0.01;     % 易感细胞细胞自然死亡率
p.c_clear = 0.1; % 胞外病毒清除率
p.Tc_generation = 5e-5;%% 抗原-抗体复合物刺激生成Tc细胞的速率 5e-5

% --- 胞内动力学参数 ---
p.v_start = 1;   % 刚感染时刻胞内病毒量
p.k5 = 0.1;     % 胞内病毒复制速率 (指数增长系数)
p.a_vec = (0:p.N-1)' * p.da; % 感染龄向量 [0, 1, 2, ..., 399]'
p.k_leak = 0; % 病毒渗漏速率
% 预计算胞内病毒量 (假设只跟感染龄有关)，这是一个 N x 1 向量
p.vin_vec = p.v_start .* exp(p.k5 .* p.a_vec);

% Remdesivir 给药设置
p.para(21) = 1e-7;
p.para(22) = 1e-14;
p.para(23) = 0.01;

p.C1_replication_ratio = 0.01;

p.drug_dose_times = [];
p.drug_dose_amounts = [];

p.drug_internal_dt = 1;
p.vin_cap = 1e50;

p.extinction_V_threshold = 1e-2;
p.extinction_I_total_threshold = 1e-2;
p.use_extinction_switch = false;
p.use_extinction_event_restart = true;
%%
% p.drug_on_time  = 2000;  % 给药时间
% p.drug_off_time = 2000;  % 停药时间
% p.drug_k5_ratio = 0.01;  % 给药期间复制速率下降为原来的 1%
% 
% % 防止 exp 数值溢出
% p.vin_exp_cap = 700;

% --- 裂解阈值参数 (Hill函数参数) ---
p.Tc_binding = 1e-5; % forward binding constant between Tc cell and infected cell
p.theta_lysis = 1e12; % 自然裂解阈值
p.theta_adcc  = 1e8;  % ADCC 结合分阈值
p.theta_tc = 1e7; % Tc 裂解的阈值
p.n_hill      = 1;    % Hill 系数 (控制平滑度/陡峭度)
p.k_lysis_max = 1;  % 超过阈值后的最大自然裂解速率
p.k_adcc_max  = 1; % 超过阈值后的最大ADCC裂解速率
p.k_tc_max = 1;% 超过阈值后的最大ADCC裂解速率

% --- Tc 细胞参数 (示例) ---
p.k_kill_tc = 0.1;   % Tc 杀伤系数

%% 2. 初始条件 (Initial Conditions)
% 状态向量结构 y = [T; I_1; ...; I_N; V; Tc; x0]
% 总长度 = 1 + N + 1 + 1 + 1401 = 1804

T0 = 1e10;                % 初始易感细胞
I0 = zeros(p.N, 1);      % 初始感染细胞 (各龄均为0)
V0 = 10;                % 初始胞外病毒
Tc0 = 1e5;               % 初始 Tc 细胞 1e5

%% 关于antibody部分：
x0 = zeros(1401,1);
%%  IgM distribution
mu = 9; % 均值
sigma = 0.8; % 标准差


% prob_A(1) = 0.001;
% prob_A(2:9) = 0;
prob_A(1) = normcdf(5, mu, sigma);
% 计算累积分布概率
for i = 2:5
prob_A(i) = normcdf(5+i-1, mu, sigma) - normcdf(5+i-2, mu, sigma);
end
for i = 6:10
prob_A(i) = prob_A(11-i);
end

for i = 1:10
    for j = 1:10
        AA(10*(i-1)+j) = prob_A(i)*prob_A(j);
    end
end

M_1 = 1e15;% 总的IgM-BCR
M_2 = 4e18;% 总的IgM
G_1 = 1e15;% 总的IgG-BCR
G_2 = 4e19;% 总的IgG
E = 1e18;% E 代表了环境类抗原的浓度，或者是自身抗原类物质的浓度
C_1 = 1e13;% 总的IgM-BCR-E 复合物
C_2 = 4e16;% 总的IgM-E 复合物
C_3 = 1e13;% 总的IgG-BCR-E 复合物
C_4 = 4e17;% 总的IgG-E 复合物
P_M = 5e7;% 总的IgM-plasma 细胞浓度
P_G = 1e8;% 总的IgG-plasma 细胞浓度


for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j) = M_1*prob_A(i)*prob_A(j);% IgM-BCR
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+100) = M_2*prob_A(i)*prob_A(j);% IgM
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+200) = G_1*prob_A(i)*prob_A(j);% IgG-BCR
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+300) = G_2*prob_A(i)*prob_A(j);% IgG
    end
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+400) = C_1*prob_A(i)*prob_A(j); % IgM-BCR-E complex
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+500) = C_2*prob_A(i)*prob_A(j);% IgM-E complex
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+600) = C_3*prob_A(i)*prob_A(j);% IgG-BCR-E complex
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+700) = C_4*prob_A(i)*prob_A(j);% IgG-E complex
    end
end

for i = 801:1200
   
 x0(i) = 0;  % IgM-BCR-V complex; IgM-V complex; IgG-BCR-V complex; IgG-V complex
  
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+1200) = P_M*prob_A(i)*prob_A(j); % plasma IgM cell
    end
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+1300) = P_G*prob_A(i)*prob_A(j); % plasma IgG cell
    end
end


x0(1401) = E;
% x0(1402) = 10;
cum_death_0 = [0; 0; 0]; % 初始化这三个计数器为0


p.AA = AA;

y0 = [T0; I0; V0; Tc0; x0; cum_death_0]; 

%% 3. 运行模拟
t_span = [0, 4000]; % 模拟时间范围

% 使用 ode15s 求解，传入参数结构体 p
non_neg_indices = 1:length(y0);

% 修改 options
options = odeset('RelTol', 1e-4, ...
                 'AbsTol', 1e-8, ...
                 'NonNegative', non_neg_indices); % 强制所有变量非负

% 调用求解器

% [t, y] = ode15s(@(t,y) sys_ode_modified_new(t, y, p), t_span, y0,options);
% % % 
% % % % Total_I = sum(y(:, 2:p.N+1), 2);
% % % 
% % % % plot(t, Total_I, 'r'); 
% % % % hold on
% % 4. Nature 风格绘图（修正布局重叠问题）
% 
% % 设置全局绘图参数
% set(groot, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 9);
% set(groot, 'DefaultLineLineWidth', 1.5);
% set(groot, 'DefaultTextFontName', 'Arial');
% 
% % 创建图形窗口 (单位厘米，18x16适合插入文档)
% figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 16]);
% 
% % --- 关键修改：使用 tiledlayout 替代 subplot ---
% % 'TileSpacing', 'compact' 让子图间距更紧凑
% % 'Padding', 'compact' 减少边缘留白
% tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
% % --------------------------
% % 子图1：易感细胞和总感染细胞 (左上)
% % --------------------------
% nexttile; % 替代 subplot
% Total_I = sum(y(:, 2:p.N+1), 2);
% plot(t, y(:,1), 'b'); hold on;
% plot(t, Total_I, 'r'); hold off;
% title('a) Cell Populations', 'FontWeight', 'bold'); % 加a) b)序号符合论文习惯
% ylabel('Cell Count', 'FontWeight', 'bold');
% legend({'Susceptible (T)', 'Infected (I)'}, 'Location', 'best', 'Box', 'off');
% grid off; box on; xlim([0 max(t)]);
% 
% % --------------------------
% % 子图2：胞外病毒 (右上)
% % --------------------------
% nexttile;
% plot(t, y(:, p.N+2), 'm');
% title('b) Extracellular Virus (V)', 'FontWeight', 'bold');
% ylabel('Virus Titer', 'FontWeight', 'bold');
% set(gca, 'YScale', 'log'); % 对数坐标
% grid off; box on; xlim([0 max(t)]);
% 
% % --------------------------
% % 子图3：Tc 细胞 (左下)
% % --------------------------
% nexttile;
% plot(t, y(:, p.N+3), 'k');
% title('c) Cytotoxic T Cells (Tc)', 'FontWeight', 'bold');
% xlabel('Time (h)', 'FontWeight', 'bold');
% ylabel('Tc Count', 'FontWeight', 'bold');
% grid off; box on; xlim([0 max(t)]);
% 
% % --------------------------
% % 子图4：细胞死亡速率 (右下)
% % --------------------------
% nexttile;
% Cum_Natural = y(:, end-2);
% Cum_ADCC    = y(:, end-1);
% Cum_Tc      = y(:, end);
% % 计算速率
% Rate_Natural = diff(Cum_Natural) ./ diff(t);
% Rate_ADCC    = diff(Cum_ADCC)    ./ diff(t);
% Rate_Tc      = diff(Cum_Tc)      ./ diff(t);
% t_mid        = (t(1:end-1) + t(2:end)) / 2;
% 
% plot(t_mid, Rate_Natural, 'k'); hold on;
% plot(t_mid, Rate_ADCC,    'b');
% plot(t_mid, Rate_Tc,      'r'); hold off;
% title('d) Dynamics of Cell Death', 'FontWeight', 'bold');
% xlabel('Time (h)', 'FontWeight', 'bold');
% ylabel('Rate (cells/h)', 'FontWeight', 'bold');
% legend({'Natural', 'ADCC', 'Tc Killing'}, 'Location', 'northeast', 'Box', 'off');
% grid off; box on; xlim([0 max(t)]);
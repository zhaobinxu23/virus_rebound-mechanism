function dydt = sys_ode_modified_new(t, y, p)
    % Unpack State Vector
    % y(1)       : T (Susceptible)
    % y(2:N+1)   : I (Infected cohorts, vector length N)
    % y(N+2)     : V (Extracellular Virus)
    % y(N+3)     : Tc (Cytotoxic T Cells)
    y = max(0,y);
    % N = p.N;
    % 
    % 
    % T  = y(1);
    % I  = y(2 : N+1); % 这是一个 N x 1 的向量，对应 I(0), I(1)... I(N-1)
    % V  = y(N+2);
    % 
    % 
    % Tc = y(N+3);
    % 
    % 
    % X  = y(N+4:end);


N = p.N;

T = y(1);
I = y(2:N+1);
V = y(N+2);
Tc = y(N+3);
X = y(N+4:N+1404);

if isfield(p, 'use_extinction_switch')
    use_extinction_switch = p.use_extinction_switch;
else
    use_extinction_switch = true;
end

if isfield(p, 'extinction_V_threshold')
    extinction_V_threshold = p.extinction_V_threshold;
else
    extinction_V_threshold = 1e-2;
end

if isfield(p, 'extinction_I_total_threshold')
    extinction_I_total_threshold = p.extinction_I_total_threshold;
else
    extinction_I_total_threshold = 1e-2;
end

I_total_now = sum(I);

is_extinct = use_extinction_switch && ...
             V < extinction_V_threshold && ...
             I_total_now < extinction_I_total_threshold;

if is_extinct
    I(:) = 0;
    V = 0;
end
    
    dydt = zeros(size(y));
    
    %% -------------------------------------------------------------
    % Part 1: 辅助变量计算 (Tc 和 Antibody 需要你根据逻辑填充)
    % -------------------------------------------------------------
    
    % [用户填充区域] 计算抗体结合分数 Phi
    % 假设 Phi 是 vin 的函数乘以当前时间的抗体水平(这里仅做简单示例)
    % 如果你有 Antibody(t) 的逻辑，请在这里计算
    % current_Ab_level = 1.0; % 示例：假设抗体水平为 1
    vin_per_cell = compute_vin_remdesivir_total_fast(t, p.a_vec, p);
        % 构造索引向量 (1到100)
    idx_base = 1:100; 
    
    % 提取需要的状态变量部分
    X_IgM = X(idx_base + 100); % 对应原来的 X(...+100)
    X_IgG = X(idx_base + 300); % 对应原来的 X(...+300)
    
    % 计算权重向量 (p.para_new 是 1x10，需要扩展成 100x1 的形式以对应 mm 循环)
    % 原逻辑 p.para_new(mm) 对每一行的 10 个 nn 是相同的
    weight_vec = repelem(p.para_new(:), 10); % 将 10个参数扩展为 100个，对应 mm=1..10
    
    % 计算 current_Ab_level (点乘后求和)
    % 原公式: sum( p.para_new(mm) * (5*IgM + IgG) )
    current_Ab_level = sum( weight_vec .* (5 * X_IgM + X_IgG) );

    Phi_vec = vin_per_cell * current_Ab_level; 
    
    % --- 同样的逻辑优化 IgM_V_IgG_V ---
    % X(...+900) 和 X(...+1100)
    X_Complex_1 = X(idx_base + 900);
    X_Complex_2 = X(idx_base + 1100);
    
    % 直接求和，不需要中间变量 array
    Complex_total = sum(X_Complex_1 + X_Complex_2);
    
    Tc_generation = Complex_total * p.Tc_generation;
    
    % [用户填充区域] Tc 动态
    % 这里假设 Tc 数量是状态变量 y(end)，不需要额外计算，直接使用 Tc 即可
    
    %% -------------------------------------------------------------
    % Part 2: 计算特定条件下的死亡率 (Hazard Functions)
    % 这里实现了“平滑化”的关键逻辑
    % -------------------------------------------------------------
    
    % A. 自然裂解率 mu_lysis (基于 Hill 方程) - Vectorized
    % 这里的逻辑：当 vin 接近 theta_lysis 时，死亡率飙升
    mu_lysis = p.k_lysis_max .* (vin_per_cell.^p.n_hill) ./ ...
               (vin_per_cell.^p.n_hill + p.theta_lysis^p.n_hill);
           
    % B. ADCC 裂解率 mu_adcc (基于 Hill 方程) - Vectorized
    mu_adcc = p.k_adcc_max .* (Phi_vec.^p.n_hill) ./ ...
              (Phi_vec.^p.n_hill + p.theta_adcc^p.n_hill);
          
    % C. Tc 杀伤率 mu_tc
    % 假设对所有感染细胞杀伤均等，也可以改成跟 antigen 表达有关
    % mu_tc = p.k_kill_tc * Tc * (vin_per_cell/p.Tc_lysis);
    Phi_vec_2 = vin_per_cell * Tc*p.Tc_binding; % 注意：这是 N x 1 向量
    mu_tc = p.k_tc_max.* (Phi_vec_2.^p.n_hill) ./ ...
              (Phi_vec_2.^p.n_hill + p.theta_tc^p.n_hill);
          
    

    
    % --- 汇总死亡率 ---
    % 1. 导致爆裂释放病毒的死亡率 (Lysis + ADCC)
    mu_burst = mu_lysis + mu_adcc;
    
    % 2. 总死亡率 (用于从 I 种群中移除细胞)
    mu_total = mu_burst + mu_tc;
    
    %% -------------------------------------------------------------
    % Part 3: 微分方程组 (RHS)
    % -------------------------------------------------------------
    
    % --- 1. 易感细胞 T ---
    % 新感染产生的速率 (流入 I 的第一个 bin)
    new_infection_rate = p.k4 * T * V / (V + p.km);
    
    dydt(1) = p.k6 - p.k7 * T - new_infection_rate;
    
    % --- 2. 感染细胞 I (PDE 离散化: Upwind Scheme) ---
    % 物理方程: dI/dt + dI/da = - mu * I
    % 离散化:   dI_j/dt = (I_{j-1} - I_j)/da - mu_j * I_j
    
    % 通量 (Flux) 计算: 代表细胞随感染龄增加从一个格子流到下一个格子
    % Flux = I / da (因为 da/dt = 1)
    flux = I ./ p.da;
    
    % 构造 dI (N x 1 向量)
    dI = zeros(N, 1);
    
    % 边界条件 (Grid j=1): 来源是 new_infection_rate
    % 注意：new_infection_rate 是总数/时间，要进入 bin 需要除以 da 变成密度，
    % 或者如果我们直接追踪“bin内的细胞数”而不是密度，则直接用速率平衡。
    % 这里假设 I 代表该 bin 内的细胞**数量**。
    dI(1) = new_infection_rate - flux(1) - mu_total(1) * I(1);
    
    % 内部格点 (Grid j=2 to N): 来源是上一格的流出 flux(j-1)
    dI(2:end) = flux(1:end-1) - flux(2:end) - mu_total(2:end) .* I(2:end);
    
    dydt(2:N+1) = dI;

    %%

 % y = max(0,y);
 F = zeros(1402,1);


for i = 1:10
    for j = 1:10
        F(10*(i-1)+j,1) = p.para(4)*p.AA(10*(i-1)+j)-p.para(1)*X(10*(i-1)+j)*X(1401)+p.para(2)*X(10*(i-1)+j+400)+p.para(7)*X(10*(i-1)+j+400)...
            -5*p.para_new(i)*X(10*(i-1)+j)*V+p.para_new_1(j)*X(10*(i-1)+j+800)...
            +1.05*p.para(15)*(1-p.para(16))*(1-p.para(17)*X(10*(i-1)+j+800)/(X(10*(i-1)+j+800)+X(10*(i-1)+j)))*X(10*(i-1)+j+800)^1.25-p.para(5)*X(10*(i-1)+j);
      
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+100,1) = -p.para(1)*X(10*(i-1)+j+100)*X(1401)+p.para(2)*X(10*(i-1)+j+500)+p.para(11)*X(10*(i-1)+j+1200)...
            -5*p.para_new(i)*X(10*(i-1)+j+100)*V+p.para_new_1(j)*X(10*(i-1)+j+900)...
            -p.para(13)*X(10*(i-1)+j+100);
        
    end
end


for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+200,1) = -p.para(1)*X(10*(i-1)+j+200)*X(1401)+p.para(2)*X(10*(i-1)+j+600)+p.para(7)*X(10*(i-1)+j+600)...
            -p.para_new(i)*X(10*(i-1)+j+200)*V+p.para_new_1(j)*X(10*(i-1)+j+1000)...
            +1.0*p.para(15)*p.para(16)*(1-p.para(17)*X(10*(i-1)+j+800)/(X(10*(i-1)+j+800)+X(10*(i-1)+j)))*X(10*(i-1)+j+800)^1.25...
            +1.0*p.para(15)*(1-p.para(17)*X(10*(i-1)+j+1000)/(X(10*(i-1)+j+1000)+X(10*(i-1)+j+200)))*X(10*(i-1)+j+1000)^1.25-p.para(6)*X(10*(i-1)+j+200);
        
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+300,1) =  -p.para(1)*X(10*(i-1)+j+300)*X(1401)+p.para(2)*X(10*(i-1)+j+700)+p.para(12)*X(10*(i-1)+j+1300)...
            -p.para_new(i)*X(10*(i-1)+j+300)*V+p.para_new_1(j)*X(10*(i-1)+j+1100)...
            -p.para(14)*X(10*(i-1)+j+300);
        

    end
end


%% c1-c4

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+400,1) =  p.para(1)*X(10*(i-1)+j)*X(1401)-p.para(2)*X(10*(i-1)+j+400)-p.para(18)*X(10*(i-1)+j+400);
        F_E_1(i,j) = -p.para(1)*X(10*(i-1)+j)*X(1401)+ p.para(2)*X(10*(i-1)+j+400);

    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+500,1) = p.para(1)*X(10*(i-1)+j+100)*X(1401)-p.para(2)*X(10*(i-1)+j+500)-p.para(18)*X(10*(i-1)+j+500);
        F_E_2(i,j) = -p.para(1)*X(10*(i-1)+j+100)*X(1401)+p.para(2)*X(10*(i-1)+j+500);
        
    end
end


for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+600,1) =  p.para(1)*X(10*(i-1)+j+200)*X(1401)-p.para(2)*X(10*(i-1)+j+600)-p.para(18)*X(10*(i-1)+j+600);
        F_E_3(i,j) = -p.para(1)*X(10*(i-1)+j+200)*X(1401)+p.para(2)*X(10*(i-1)+j+600);

    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+700,1) = p.para(1)*X(10*(i-1)+j+300)*X(1401)-p.para(2)*X(10*(i-1)+j+700)-p.para(18)*X(10*(i-1)+j+700);
        F_E_4(i,j) = -p.para(1)*X(10*(i-1)+j+300)*X(1401)+p.para(2)*X(10*(i-1)+j+700);
    end
end
%% c1'-c4'

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+800,1) = 5*p.para_new(i)*X(10*(i-1)+j)*V-p.para_new_1(j)*X(10*(i-1)+j+800)-p.para(18)*X(10*(i-1)+j+800);
        F_V_1(i,j) = -5*p.para_new(i)*X(10*(i-1)+j)*V+p.para_new_1(j)*X(10*(i-1)+j+800);    
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+900,1) = 5*p.para_new(i)*X(10*(i-1)+j+100)*V-p.para_new_1(j)*X(10*(i-1)+j+900)-p.para(18)*X(10*(i-1)+j+900);
        F_V_2(i,j) = -5*p.para_new(i)*X(10*(i-1)+j+100)*V+p.para_new_1(j)*X(10*(i-1)+j+900);         
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1000,1) = p.para_new(i)*X(10*(i-1)+j+200)*V-p.para_new_1(j)*X(10*(i-1)+j+1000)-p.para(18)*X(10*(i-1)+j+1000);
        F_V_3(i,j) = -p.para_new(i)*X(10*(i-1)+j+200)*V+p.para_new_1(j)*X(10*(i-1)+j+1000);      
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1100,1) = p.para_new(i)*X(10*(i-1)+j+300)*V-p.para_new_1(j)*X(10*(i-1)+j+1100)-p.para(18)*X(10*(i-1)+j+1100);
        F_V_4(i,j) = -p.para_new(i)*X(10*(i-1)+j+300)*V+p.para_new_1(j)*X(10*(i-1)+j+1100);      
    end
end

%% plasma_M
for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1200,1) = p.para(8)*X(10*(i-1)+j+400)+1.0*p.para(15)*(1-p.para(16))*p.para(17)*X(10*(i-1)+j+800)/(X(10*(i-1)+j+800)+X(10*(i-1)+j))*X(10*(i-1)+j+800)^1.25/p.para(20)-p.para(9)*X(10*(i-1)+j+1200);
             
    end
end

%% plasma_G
for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1300,1) = p.para(8)*X(10*(i-1)+j+600)+1.0*p.para(15)*p.para(17)*X(10*(i-1)+j+1000)/(X(10*(i-1)+j+1000)+X(10*(i-1)+j+200))*X(10*(i-1)+j+1000)^1.25/p.para(20)-p.para(10)*X(10*(i-1)+j+1300);
             
    end
end





F(1401,1) = p.para(3)+sum(sum(F_E_1))+sum(sum(F_E_2))+sum(sum(F_E_3))+sum(sum(F_E_4));
F(1402,1) = p.para(19)*V+sum(sum(F_V_1))+sum(sum(F_V_2))+sum(sum(F_V_3))+sum(sum(F_V_4));
    % --- 3. 胞外病毒 V (积分项的实现) ---
    % 公式: Source = Integral[ i(a) * mu_burst(a) * vin(a) da ]
    % 离散化: Sum[ I_j * mu_burst_j * vin_j ]
    
    % 计算每一组此刻释放病毒的速率
    % I .* mu_burst 计算的是“此刻每组有多少个细胞正在爆裂”
    viral_release_per_bin = I .* mu_total .* vin_per_cell; 

    
    % 对所有 bin 求和 -> 得到总释放率
    total_viral_source = sum(viral_release_per_bin);
    Secretion_Flux = sum(I .* vin_per_cell .* p.k_leak);
    
    % dV/dt
    dydt(N+2) = total_viral_source + F(1402,1) + Secretion_Flux - p.k4 * T * V / (V + p.km) - p.c_clear*V;

% --- 4. Tc 细胞 ---
    % [用户填充区域]

    
dTc_dt = sum(-0.1*I .* mu_tc) + Tc_generation; % sum(-0.1*I .* mu_tc) 表示Tc细胞的exhuastion
dydt(N+3) = dTc_dt;
dydt(N+4:end-3) = F(1:1401,1);
%% --- [新增] 追踪细胞裂解累积量 ---
% 计算当前时刻三种机制造成的瞬时死亡总数 (所有感染龄求和)
rate_death_lysis = sum(I .* mu_lysis);
rate_death_adcc  = sum(I .* mu_adcc);
rate_death_tc    = sum(I .* mu_tc);

% 将这些速率记录到新的状态变量的导数中
% 假设原来的 y 长度是 len_old
% 我们在最后追加 3 个变量: 
% Cum_Lysis (自然), Cum_ADCC, Cum_Tc

% 注意：你需要确保 sys_ode 的输出向量维度会自动扩展，或者手动指定索引
% 原来的最后索引是 N + 3 + 1401
idx_end = N + 3 + 1401; 

dydt(idx_end + 1) = rate_death_lysis; % 累积自然裂解数
dydt(idx_end + 2) = rate_death_adcc;  % 累积 ADCC 裂解数
dydt(idx_end + 3) = rate_death_tc;    % 累积 Tc 裂解数

% I_idx = 2:N+1;
% V_idx = N+2;
% 
% below_I = I_raw < I_threshold;
% dI_now = dydt(I_idx);
% 
% dI_now(below_I & dI_now < 0) = 0;
% dydt(I_idx) = dI_now;
% 
% if V_raw < V_threshold && dydt(V_idx) < 0
%     dydt(V_idx) = 0;
% end


% if is_extinct
%     dydt(2:N+1) = 0;
%     dydt(N+2) = 0;
% end

if size(dydt,2) > 1
    dydt = dydt';
end

end
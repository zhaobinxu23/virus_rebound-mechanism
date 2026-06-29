function dist = compute_affinity_distribution_cell(t, y, t_grid)

% Interpolate IgM and IgG trajectories onto a uniform time grid
nT = length(t_grid);

data_IgM = zeros(100, nT);   % soluble IgM: 101:200
data_IgG = zeros(100, nT);   % soluble IgG: 301:400

for i = 1:100
    data_IgM(i,:) = interp1(t, y(:,i+403+100), t_grid, 'linear', 'extrap');
    data_IgG(i,:) = interp1(t, y(:,i+403+300), t_grid, 'linear', 'extrap');
end

% Aggregate by kon
data_k_on_IgM = zeros(10, nT);
data_k_on_IgG = zeros(10, nT);

for i = 1:100
    row_index = fix((i - 1)/10) + 1;
    data_k_on_IgM(row_index,:) = data_k_on_IgM(row_index,:) + data_IgM(i,:);
    data_k_on_IgG(row_index,:) = data_k_on_IgG(row_index,:) + data_IgG(i,:);
end

% Aggregate by koff
data_k_off_IgM = zeros(10, nT);
data_k_off_IgG = zeros(10, nT);

for i = 1:100
    row_index = mod(i - 1, 10) + 1;
    data_k_off_IgM(row_index,:) = data_k_off_IgM(row_index,:) + data_IgM(i,:);
    data_k_off_IgG(row_index,:) = data_k_off_IgG(row_index,:) + data_IgG(i,:);
end

% Aggregate by binding-energy / Kd class
data_kd_IgM = zeros(19, nT);
data_kd_IgG = zeros(19, nT);

for i = 1:100
    row_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;
    data_kd_IgM(row_index,:) = data_kd_IgM(row_index,:) + data_IgM(i,:);
    data_kd_IgG(row_index,:) = data_kd_IgG(row_index,:) + data_IgG(i,:);
end

dist.t_grid = t_grid;
dist.data_IgM = data_IgM;
dist.data_IgG = data_IgG;

dist.k_on_IgM = data_k_on_IgM;
dist.k_on_IgG = data_k_on_IgG;

dist.k_off_IgM = data_k_off_IgM;
dist.k_off_IgG = data_k_off_IgG;

dist.kd_IgM = data_kd_IgM;
dist.kd_IgG = data_kd_IgG;

% binding-energy class labels
dist.kd_values = -31:1:-13;

end
function vin_total = compute_vin_remdesivir_total_fast(t, a_vec, p)

    a_vec = a_vec(:);
    N = length(a_vec);

    vin_total = zeros(N, 1);

    k5    = p.k5;
    kon   = p.para(21);
    koff  = p.para(22);

    if isfield(p, 'k_leak')
        k_leak = p.k_leak;
    else
        k_leak = 0;
    end

    if isfield(p, 'C1_replication_ratio')
        rho_C = p.C1_replication_ratio;
    else
        rho_C = 0.01;
    end

    if isfield(p, 'da')
        da = p.da;
    else
        da = a_vec(2) - a_vec(1);
    end

    if isfield(p, 'vin_cap')
        vin_cap = p.vin_cap;
    else
        vin_cap = 1e300;
    end

    D = drug_concentration_fast(t, p);

    A = [k5 - kon * D - k_leak,      koff;
         kon * D,                   rho_C * k5 - koff - k_leak];

    M = expm(A * da);

    y = [p.v_start; 0];

    for j = 1:N
        if j > 1
            y = M * y;
        end

        y(y < 0) = 0;
        y(y > vin_cap) = vin_cap;

        vin_total(j) = y(1) + y(2);
    end

    vin_total(vin_total > vin_cap) = vin_cap;
end
function D = drug_concentration_fast(t, p)

    if ~isfield(p, 'drug_dose_times') || isempty(p.drug_dose_times)
        D = 0;
        return;
    end

    dose_times = p.drug_dose_times(:);

    if isfield(p, 'drug_dose_amounts')
        dose_amounts = p.drug_dose_amounts(:);
    else
        dose_amounts = ones(size(dose_times));
    end

    kD = p.para(23);

    idx = dose_times <= t;

    if ~any(idx)
        D = 0;
        return;
    end

    dt = t - dose_times(idx);

    D = sum(dose_amounts(idx) .* exp(-kD .* dt));
end
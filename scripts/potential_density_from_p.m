function sigma = potential_density_from_p(S, T, P, Lon, Lat)
    sa = gsw_SA_from_SP(S, P, Lon, Lat);
    ct = gsw_CT_from_t(sa, T, P);
    sigma = gsw_sigma0(sa, ct);
end

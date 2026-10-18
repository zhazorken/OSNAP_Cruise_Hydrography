% Potential temperature from p
function pt = potential_temperature_from_p(S, T, P, Lon, Lat)
    sa = gsw_SA_from_SP(S, P, Lon, Lat);
    pt = gsw_pt_from_t(sa, T, P, 0.0);
end

% maketransect_from2014_2022mean.m - Plot mean hydrographic transects (2014-2022)
%
% DESCRIPTION:
%   This script loads the multi-year mean gridded hydrographic product and
%   creates publication-quality plots of Temperature, Salinity, and velocity
%   sections for both OSNAP West and OSNAP East transects. It also performs
%   transport calculations in salinity space.
%
% REQUIRED INPUT FILE:
%   - meangrid.mat (created by hydromean.m)
%
% REQUIRED DEPENDENCIES:
%   - GSW Oceanographic Toolbox
%   - M_Map toolbox
%   - cmocean colormaps
%   - smooth2a function
%   - Functions: fillNaNsrow, fillNaNsWithDecay, findClosestElement
%


%% Load mean gridded data
load meangrid.mat

%% Plot mean hydrographic sections
figure

%==========================================================================
% OSNAP WEST - Temperature
%==========================================================================
subplot(3,2,1)
contourf(xnew, zvec, wbathmasknew.*Twest, 100, 'linecolor', 'none')
shading flat
cmocean('thermal')
caxis([-1 8])
hold on
% Density contours
[CC,h] = contour(xnew, zvec, wbathmasknew.*sigwest, [27:.1:29], 'k');
clabel(CC, h, 'LabelSpacing', 600);
[CC2,h2] = contour(xnew, zvec, wbathmasknew.*sigwest, [25:.5:27.8], 'r', 'linewidth', 2);
clabel(CC2, h2, 'LabelSpacing', 600);
% Station positions
x = xx_westsort';
plot(x, 30*ones(size(x)), 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% Mooring positions
plot(xx_moorw, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% Bathymetry
plot(xx_westfinesort, wtopotransectfine, 'k', 'linewidth', 2)
set(gca, 'Xdir', 'reverse')
xticklabels([])
ylabel('Depth (m)')
title(strcat(year, ' \Theta (OSNAP WEST)'))
ylim([-1400 50])
xlim([0 90])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty

%==========================================================================
% OSNAP WEST - Salinity
%==========================================================================
subplot(3,2,3)
contourf(xnew, zvec, wbathmasknew.*Swest, [10:.01:35.2], 'linecolor', 'none')
shading flat
cmocean('haline')
caxis([33 35.2])
hold on
% Density contours
[CC,h] = contour(xnew, zvec, wbathmasknew.*sigwest, [27:.1:29], 'k');
clabel(CC, h, 'LabelSpacing', 600);
[CC2,h2] = contour(xnew, zvec, wbathmasknew.*sigwest, [25:.5:27.8], 'r', 'linewidth', 2);
clabel(CC2, h2, 'LabelSpacing', 600);
% Station positions
x = xx_westsort';
plot(x, 30*ones(size(x)), 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% Mooring positions
plot(xx_moorw, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% Bathymetry
plot(xx_westfinesort, wtopotransectfine, 'k', 'linewidth', 2)
set(gca, 'Xdir', 'reverse')
xticklabels([])
ylabel('Depth (m)')
title('S_A')
ylim([-1400 50])
xlim([0 90])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty

%==========================================================================
% OSNAP EAST - Temperature
%==========================================================================
subplot(3,2,2)
contourf(xnew, zvec, ebathmasknew.*Teast, 100, 'linecolor', 'none')
colorbar
shading flat
cmocean('thermal')
yticklabels([])
caxis([-1 8])
hold on
% Density contours
[CC,h] = contour(xnew, zvec, ebathmasknew.*sigeast, [27:.1:29], 'k');
clabel(CC, h, 'LabelSpacing', 600);
[CC2,h2] = contour(xnew, zvec, ebathmasknew.*sigeast, [25:.5:27.8], 'r', 'linewidth', 2);
clabel(CC2, h2, 'LabelSpacing', 600);
% Station positions
x = xx_eastsort';
plot(x, 30*ones(size(x)), 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% Mooring positions
plot(xx_moor, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% Bathymetry
plot(xx_eastfinesort, etopotransectfine, 'k', 'linewidth', 2)
xticklabels([])
title(strcat(year, ' \Theta (OSNAP EAST)'))
ylim([-1400 50])
xlim([0 90])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty

%==========================================================================
% OSNAP EAST - Salinity
%==========================================================================
subplot(3,2,4)
contourf(xnew, zvec, ebathmasknew.*Seast, [10:.01:35.2], 'linecolor', 'none')
colorbar
shading flat
cmocean('haline')
yticklabels([])
caxis([33 35.2])
hold on
% Density contours
[CC,h] = contour(xnew, zvec, ebathmasknew.*sigeast, [27:.1:29], 'k');
clabel(CC, h, 'LabelSpacing', 600);
[CC2,h2] = contour(xnew, zvec, ebathmasknew.*sigeast, [25:.5:27.8], 'r', 'linewidth', 2);
clabel(CC2, h2, 'LabelSpacing', 600);
% Station positions
x = xx_eastsort';
plot(x, 30*ones(size(x)), 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% Mooring positions
plot(xx_moor, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% Bathymetry
plot(xx_eastfinesort, etopotransectfine, 'k', 'linewidth', 2)
xticklabels([])
title('S_A')
ylim([-1400 50])
xlim([0 90])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty

%==========================================================================
% OSNAP WEST - Velocity
%==========================================================================
subplot(3,2,5)
contourf(xneww_v, zvec, wbathmasknew(:,:).*vwest', 100, 'linecolor', 'none')
shading flat
caxis([0 0.6])
hold on
% Density contours
[CC,h] = contour(xnew, zvec, wbathmasknew.*sigwest, [27:.1:29], 'k');
clabel(CC, h, 'LabelSpacing', 600);
[CC2,h2] = contour(xnew, zvec, wbathmasknew.*sigwest, [25:.5:27.8], 'r', 'linewidth', 2);
clabel(CC2, h2, 'LabelSpacing', 600);
colormap(gca, pmkmp(100,'swtth'))
% Station positions
x = xx_westsort';
plot(x, 30*ones(size(x)), 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% Mooring positions
plot(xx_moorw, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% Bathymetry
plot(xx_westfinesort, wtopotransectfine, 'k', 'linewidth', 2)
set(gca, 'Xdir', 'reverse')
xlabel('Dist. (km)')
ylabel('Depth (m)')
title('Velocity [m/s]')
ylim([-1400 50])
xlim([0 90])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty

%==========================================================================
% OSNAP EAST - Velocity
%==========================================================================
subplot(3,2,6)
contourf(xnewe_v, zvec, ebathmasknew(:,:).*veast', 100, 'linecolor', 'none')
colorbar
shading flat
yticklabels([])
caxis([0 .6])
hold on
% Density contours
[CC,h] = contour(xnew, zvec, ebathmasknew.*sigeast, [27:.1:29], 'k');
clabel(CC, h, 'LabelSpacing', 600);
[CC2,h2] = contour(xnew, zvec, ebathmasknew.*sigeast, [25:.5:27.8], 'r', 'linewidth', 2);
clabel(CC2, h2, 'LabelSpacing', 600);
colormap(gca, pmkmp(100,'swtth'))
% Station positions
x = xx_eastsort';
plot(x, 30*ones(size(x)), 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% Mooring positions
plot(xx_moor, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% Bathymetry
plot(xx_eastfinesort, etopotransectfine, 'k', 'linewidth', 2)
xlabel('Dist. (km)')
title('Velocity [m/s]')
ylim([-1400 50])
xlim([0 90])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty

%% Define parameters and add paths
westlon_orig = -45.2654; 
westlat_orig = 59.9854;
eastlon_orig = -43.1189; 
eastlat_orig = 60.1022;

% Add required paths (update to match your directory structure)
addpath ~/Desktop/matlab/gsw_matlab_v3_06_16/
addpath ~/Desktop/matlab/gsw_matlab_v3_06_16/library/
addpath ~/Desktop/SIO_2024/Work/matlab_scripts/osnap_calc
addpath ~/Desktop/matlab
addpath ~/Desktop/matlab/m_map

%% Set grid parameters
zmax = 2000;
zvec = 2:2:2*zmax;
xres = 0.25;

%% Define mooring positions
ii = 1; eastlatmoor(ii) = 59 + 55.367/60; eastlonmoor(ii) = 41 + 26.02/60; eastz(ii) = 1901;
ii = 2; eastlatmoor(ii) = 59 + 57.277/60; eastlonmoor(ii) = 41 + 44.648/60; eastz(ii) = 1829;
ii = 3; eastlatmoor(ii) = 59 + 59.087/60; eastlonmoor(ii) = 42 + 1.563/60; eastz(ii) = 1260;
ii = 4; eastlatmoor(ii) = 60 + 0.302/60; eastlonmoor(ii) = 42 + 12.34/60; eastz(ii) = 384;
ii = 5; eastlatmoor(ii) = 60 + 1.851/60; eastlonmoor(ii) = 42 + 25.708/60; eastz(ii) = 184;
ii = 6; eastlatmoor(ii) = 60 + 2.853/60; eastlonmoor(ii) = 42 + 35.975/60; eastz(ii) = 178;
ii = 7; eastlatmoor(ii) = 60 + 4.208/60; eastlonmoor(ii) = 42 + 49.527/60; eastz(ii) = 170;

for ii = 1:7
    xx_moor(ii) = m_lldist([-eastlonmoor(ii) eastlon_orig], [eastlatmoor(ii) eastlat_orig]);
end

%% Re-grid data onto standardized fine grid for transport calculations
% New grid with finer horizontal resolution
xnewfine = 0:0.1:ceil(max(xnew));  % 100 m horizontal resolution

zmax = 2000;
zvecfine = 2:2:2*zmax;

% Create grids for interpolation
[X, Z] = meshgrid(xnew, zvec);
[Xfine, Zfine] = meshgrid(xnewfine, -zvecfine);

% Interpolate data onto fine grid
sigma_gridded = interp2(X, Z, sigeast, Xfine, Zfine, 'linear');
temp_gridded = interp2(X, Z, Teast, Xfine, Zfine, 'linear');
salt_gridded = interp2(X, Z, Seast, Xfine, Zfine, 'linear');

[X, Zadcp] = meshgrid(xnewe_v, zvec);
[Xfine, Zfine] = meshgrid(xnewfine, -zvecfine);
v_gridded = interp2(X, Zadcp, veast', Xfine, Zfine, 'linear');

%% Create bathymetry mask for fine grid
% Interpolate lat/lon positions onto fine grid
latvec = interp1(xnew, eastlat, xnewfine);
lonvec = interp1(xnew, eastlon, xnewfine);

% Extract topography along transect
for ii = 1:numel(latvec)
    [latind_interp(ii), latvecinterp(ii), lonind_interp(ii), lonvecinterp(ii)] = ...
        findClosestElement(latgridbox, latvec(ii), longridbox, lonvec(ii));
    wtopotransect(ii) = topobox(lonind_interp(ii), latind_interp(ii));
end

% Initialize the bathymetry mask with NaNs
ebathmask = nan(numel(zvecfine), numel(xnewfine));

% Create mask (NaN = water, 1 = below seafloor)
for ii = 1:size(ebathmask, 2)
    MM = ebathmask(:, ii);
    MM(-zvecfine' <= wtopotransect(ii)) = NaN;  % Water column
    MM(-zvecfine' > wtopotransect(ii)) = 1;     % Below seafloor
    ebathmasknew(:, ii) = MM;
end

%% Fill NaN values and smooth
sigma_gridded_f = fillNaNsrow(sigma_gridded);
temp_gridded_f = fillNaNsrow(temp_gridded);
salt_gridded_f = fillNaNsrow(salt_gridded);
v_gridded_f = fillNaNsWithDecay(v_gridded, 60, 100);

sigma_gridded_s = fillNaNsrow(smooth2a(sigma_gridded_f, 5, 5)) .* ebathmasknew;
temp_gridded_s = fillNaNsrow(smooth2a(temp_gridded_f, 5, 5)) .* ebathmasknew;
salt_gridded_s = fillNaNsrow(smooth2a(salt_gridded_f, 5, 5)) .* ebathmasknew;
v_gridded_s = fillNaNsrow(smooth2a(v_gridded_f, 5, 5)) .* ebathmasknew;

%% Plot gridded and smoothed data
figure

% Temperature
subplot(4,1,1)
contourf(xnewfine, -zvecfine, temp_gridded_s, 100, 'linecolor', 'none')
colorbar
shading flat
cmocean('thermal')
caxis([-1.5 4])
hold on
[CC,h] = contour(xnewfine, -zvecfine, sigma_gridded_s, [27:.1:29], 'k');
clabel(CC, h, 'LabelSpacing', 600);
[CC2,h2] = contour(xnewfine, -zvecfine, sigma_gridded_s, [25:.5:27.8], 'r', 'linewidth', 2);
clabel(CC2, h2, 'LabelSpacing', 600);
x = xx_eastsort';
plot(x, 10*ones(size(x)), 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
plot(xx_moor, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
xticklabels([])
ylabel('Depth (m)')
title('Stat. 105-116, \Theta')
ylim([-200 30])
xlim([0 50])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty

% Salinity
subplot(4,1,2)
contourf(xnewfine, -zvecfine, salt_gridded_s, [10:.01:35.2], 'linecolor', 'none')
colorbar
shading flat
cmocean('haline')
caxis([31 35.2])
hold on
[CC,h] = contour(xnewfine, -zvecfine, sigma_gridded_s, [27:.1:29], 'k');
clabel(CC, h, 'LabelSpacing', 600);
[CC2,h2] = contour(xnewfine, -zvecfine, sigma_gridded_s, [25:.5:27.8], 'r', 'linewidth', 2);
clabel(CC2, h2, 'LabelSpacing', 600);
x = xx_eastsort';
plot(x, 10*ones(size(x)), 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
plot(xx_moor, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
xticklabels([])
ylabel('Depth (m)')
title('S_A')
ylim([-200 30])
xlim([0 50])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty

% Velocity
subplot(4,1,3)
contourf(xnewfine, -zvecfine, v_gridded_s, 100, 'linecolor', 'none')
colorbar
shading flat
caxis([-.6 .6])
cmocean('balance')
hold on
[CC,h] = contour(xnewfine, -zvecfine, sigma_gridded_s, [27:.1:29], 'k');
clabel(CC, h, 'LabelSpacing', 600);
[CC2,h2] = contour(xnewfine, -zvecfine, sigma_gridded_s, [25:.5:27.8], 'r', 'linewidth', 2);
clabel(CC2, h2, 'LabelSpacing', 600);
x = xx_eastsort';
plot(x, 10*ones(size(x)), 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
plot(xx_moor, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
ylabel('Depth (m)')
title('Velocity (m/s)')
ylim([-200 30])
xlim([0 50])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty

%% Calculate freshwater transport in salinity space
clear integrand_interp filtered_integrand_interp

% Run external scripts for freshwater flux calculations
Sflux;      % Main freshwater flux calculation
Sfluxalt;   % Alternative calculation for south transect

% Calculate total freshwater transport
dSbins = S_bins(2) - S_bins(1);
dx = (xnewfine(2) - xnewfine(1)) * 1000;  % Convert km to m
integral = sum(sum(integrand_interp(:,:))) * dSbins * dx;

% Plot freshwater transport
subplot(4,1,4)
filtered_integrand_interp = integrand_interp;
filtered_integrand_interp(integrand_interp > 5) = 5;  % Cap for visualization
pcolor(xnewfine, S_bins, filtered_integrand_interp)
colorbar
shading interp
colormap(pmkmp(40,'swtth'))
ylabel('Abs. Salinity (g/kg)')
xlabel('x (km)')
title('Freshwater Transport (m^2/s/(g/kg))')
ylim([29.5 35])
xlim([0 50])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty
caxis([0 5])
set(gca, 'YDir', 'reverse');

%% Plot freshwater flux profiles
figure;

% Freshwater flux as function of salinity
subplot(2,1,1)
plot(sum(integrand_interp(:,1:min(numel(xnewfine),500)),2)*dx, S_bins, 'b-', 'LineWidth', 2);
ylabel('Salinity (g/kg)');
xlabel('v(S) * ((S_{ref} - S) / S_{ref}) * (dz / dS)');
title('FWF(S)');
grid on;
set(gca, 'YDir', 'reverse');

% Cumulative freshwater flux
subplot(2,1,2)
cumulative_integral = cumsum((S_bins(2)-S_bins(1)) * sum(integrand_interp(:,1:min(numel(xnewfine),500))*dx, 2));
plot(cumulative_integral, S_bins, 'r-', 'LineWidth', 2);
ylabel('Salinity (g/kg)');
xlabel('Cumulative Integral');
title('Cumulative Integral FWF(S)');
grid on;

%% Calculate volume transport
transport;  % Run external transport calculation script

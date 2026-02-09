% maketransect.m - Create hydrographic transect from CTD and ADCP data
%
% DESCRIPTION:
%   This script processes CTD (Temperature and Salinity) and ADCP (velocity)
%   data to create gridded hydrographic transects. The workflow includes:
%   1. Loading and preprocessing raw CTD/ADCP data
%   2. Computing derived variables (potential density, absolute salinity, etc.)
%   3. Sorting stations by distance along transect
%   4. Gridding data onto regular horizontal and vertical grids
%   5. Filling missing values using specialized interpolation
%   6. Creating bathymetry masks
%   7. Generating diagnostic plots
%   8. Computing freshwater and volume transports
%
% REQUIRED DEPENDENCIES:
%   - GSW Oceanographic Toolbox (gsw_matlab_v3_06_16)
%   - M_Map toolbox
%   - Custom functions: fillNaNsrow, fillNaNsWithDecay, findClosestElement
%   - cmocean colormaps
%   - smooth2a function
%
% DATA REQUIREMENTS:
%   - CTD data: tempmat, salmat, depthmat, lonmat, latmat
%   - ADCP data: u2, v2, depthadcp2
%   - Bathymetry data: topobox, latgridbox, longridbox
%
% OUTPUTS:
%   - Gridded Temperature, Salinity, density, and velocity fields
%   - Diagnostic plots of transect sections
%   - Transport calculations
%


%% Define transect endpoints
westlon_orig = -45.2654; 
westlat_orig = 59.9854;
eastlon_orig = -43.1189; 
eastlat_orig = 60.1022;

%% Add required paths
% NOTE: Update these paths to match your local directory structure
addpath ~/Desktop/matlab/gsw_matlab_v3_06_16/
addpath ~/Desktop/matlab/gsw_matlab_v3_06_16/library/
addpath ~/Desktop/SIO_2024/Work/matlab_scripts/osnap_calc
addpath ~/Desktop/matlab
addpath ~/Desktop/matlab/m_map

%% Compute derived oceanographic variables
% Potential density from in-situ measurements
sigma = potential_density_from_p(salmat, tempmat, depthmat, lonmat, latmat);

% Potential temperature
pt = potential_temperature_from_p(salmat, tempmat, depthmat, lonmat, latmat);

% Absolute salinity using TEOS-10
sa = gsw_SA_from_SP(salmat, depthmat, lonmat, latmat);

%% Set vertical grid parameters
zmax = 2000;  % Maximum depth (m)
zvec = 2:2:2*zmax;  % Vertical grid: 2m resolution from 2m to 4000m

%% Set horizontal resolution
xres = 0.25;  % Horizontal resolution (km)

%% Define mooring positions along east transect
% Mooring locations specified in degrees and decimal minutes
ii = 1; eastlatmoor(ii) = 59 + 55.367/60; eastlonmoor(ii) = 41 + 26.02/60; eastz(ii) = 1901;
ii = 2; eastlatmoor(ii) = 59 + 57.277/60; eastlonmoor(ii) = 41 + 44.648/60; eastz(ii) = 1829;
ii = 3; eastlatmoor(ii) = 59 + 59.087/60; eastlonmoor(ii) = 42 + 1.563/60; eastz(ii) = 1260;
ii = 4; eastlatmoor(ii) = 60 + 0.302/60; eastlonmoor(ii) = 42 + 12.34/60; eastz(ii) = 384;
ii = 5; eastlatmoor(ii) = 60 + 1.851/60; eastlonmoor(ii) = 42 + 25.708/60; eastz(ii) = 184;
ii = 6; eastlatmoor(ii) = 60 + 2.853/60; eastlonmoor(ii) = 42 + 35.975/60; eastz(ii) = 178;
ii = 7; eastlatmoor(ii) = 60 + 4.208/60; eastlonmoor(ii) = 42 + 49.527/60; eastz(ii) = 170;

% Calculate mooring distances from transect origin
for ii = 1:7
    xx_moor(ii) = m_lldist([-eastlonmoor(ii) eastlon_orig], [eastlatmoor(ii) eastlat_orig]);
end

%% Select CTD stations for transect
% Choose one of the following transect options:

% Option 1: Mid transect (stations 85-98)
eastil = 85; 
eastir = 98;
eastlon_orig = -43.085; 
eastlat_orig = 60.1384;

% Option 2: North/trough transect (stations 105-116)
% eastil = 105;
% eastir = 116;
% eastlon_orig = -43.0975;
% eastlat_orig = 60.2343;

% Option 3: South transect (stations 65-75, 81-83)
% indseceast = [65:75 81:83];
% eastlon_orig = -43.11;
% eastlat_orig = 60.0967;
% eastlon = [lonmat(65:75) lonmat(81:83)]';
% eastlat = [latmat(65:75) latmat(81:83)]';

%% Extract station coordinates
indseceast = [eastil:eastir];
eastlon = [lonmat(eastil:eastir)]';
eastlat = [latmat(eastil:eastir)]';

%% Calculate station distances from transect origin
clear xx_east

for ii = 1:size(eastlon, 1)
    xx_east(ii) = m_lldist([eastlon(ii) eastlon_orig], [eastlat(ii) eastlat_orig]);
end

% Sort stations by distance
[xx_eastsort, Ie] = sort(xx_east);
indseceastsort = indseceast(Ie);
xnew = xx_eastsort;

%% Plot raw CTD and velocity data
zindmax = 1000;  % Maximum depth for plotting
figure

% Temperature plot
subplot(3,1,1)
contourf(xnew, -zvec, tempmat(:,indseceastsort), 100, 'linecolor', 'none')
colorbar
shading flat
cmocean('thermal')
caxis([-1.5 4])
hold on
[CC,h] = contour(xnew, -zvec, sigma(:,indseceastsort), [27:.1:29], 'k');
clabel(CC, h, 'LabelSpacing', 600);
[CC2,h2] = contour(xnew, -zvec, sigma(:,indseceastsort), [25:.5:27.8], 'r', 'linewidth', 2);
clabel(CC2, h2, 'LabelSpacing', 600);
x = xx_eastsort';
plot(x, 10*ones(size(x)), 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
xticklabels([])
ylabel('Depth (m)')
title('Stat. 105-116, \Theta')
ylim([-200 30])
xlim([0 50])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty
plot(xx_moor, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

% Salinity plot
subplot(3,1,2)
contourf(xnew, -zvec, salmat(:,indseceastsort), [10:.01:35.2], 'linecolor', 'none')
colorbar
shading flat
cmocean('haline')
caxis([31 35.2])
hold on
[CC,h] = contour(xnew, -zvec, sigma(:,indseceastsort), [27:.1:29], 'k');
clabel(CC, h, 'LabelSpacing', 600);
[CC2,h2] = contour(xnew, -zvec, sigma(:,indseceastsort), [25:.5:27.8], 'r', 'linewidth', 2);
clabel(CC2, h2, 'LabelSpacing', 600);
x = xx_eastsort';
plot(x, 10*ones(size(x)), 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
xticklabels([])
ylabel('Depth (m)')
title('S_A')
ylim([-200 30])
xlim([0 50])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty
plot(xx_moor, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

%% Calculate velocity component perpendicular to transect
% Transect angle relative to north
angle_deg = -4.5383;
angle_rad = deg2rad(angle_deg);

% Unit vector along transect
unit_vector_x = cos(angle_rad);
unit_vector_y = sin(angle_rad);
unit_vector = [unit_vector_x, unit_vector_y];

% Perpendicular unit vector (normal to transect)
perp_unit_vector = [unit_vector_y -unit_vector_x];

% Extract velocity components for selected stations
U2 = u2(:, indseceastsort);
V2 = v2(:, indseceastsort);

clear projection
% Project velocity onto perpendicular direction
for i = 1:size(U2, 1)
    for j = 1:size(U2, 2)
        current_vector = [U2(i, j), V2(i, j)];
        projection(i, j) = dot(current_vector, perp_unit_vector);
    end
end

vel = projection;

% Velocity plot
subplot(3,1,3)
contourf(xnew, -depthadcp2(:,1), projection, 100, 'linecolor', 'none')
colorbar
shading flat
caxis([-.6 .6])
cmocean('balance')
hold on
[CC,h] = contour(xnew, -zvec, sigma(:,indseceastsort), [27:.1:29], 'k');
clabel(CC, h, 'LabelSpacing', 600);
[CC2,h2] = contour(xnew, -zvec, sigma(:,indseceastsort), [25:.5:27.8], 'r', 'linewidth', 2);
clabel(CC2, h2, 'LabelSpacing', 600);
x = xx_eastsort';
plot(x, 10*ones(size(x)), 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
ylabel('Depth (m)')
title('Velocity (m/s)')
ylim([-200 30])
xlim([0 50])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty
plot(xx_moor, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

%% Grid data onto fine resolution grid
% New grid with finer horizontal resolution
xnewfine = 0:0.1:ceil(max(xnew));  % 100 m horizontal resolution
zvecfine = zvec;  % Maintain 2 m vertical resolution

% Create grids for interpolation
[X, Z] = meshgrid(xnew, -zvec);
[Xfine, Zfine] = meshgrid(xnewfine, -zvecfine);

% Interpolate CTD data onto fine grid
sigma_gridded = interp2(X, Z, sigma(:,indseceastsort), Xfine, Zfine, 'linear');
temp_gridded = interp2(X, Z, tempmat(:,indseceastsort), Xfine, Zfine, 'linear');
salt_gridded = interp2(X, Z, salmat(:,indseceastsort), Xfine, Zfine, 'linear');

% Interpolate velocity data onto fine grid
[X, Zadcp] = meshgrid(xnew, -depthadcp2(:,1));
[Xfine, Zfine] = meshgrid(xnewfine, -zvecfine);
v_gridded = interp2(X, Zadcp, vel, Xfine, Zfine, 'linear');

%% Create bathymetry mask
% Interpolate lat/lon positions onto fine grid
latvec = interp1(xnew, eastlat, xnewfine);
lonvec = interp1(xnew, eastlon, xnewfine);

% Extract topography along transect
for ii = 1:numel(latvec)
    [latind_interp(ii), latvecinterp(ii), lonind_interp(ii), lonvecinterp(ii)] = ...
        findClosestElement(latgridbox, latvec(ii), longridbox, lonvec(ii));
    wtopotransect(ii) = topobox(lonind_interp(ii), latind_interp(ii));
end

% Create mask (1 = water, NaN = below seafloor)
ebathmask = nan(numel(zvecfine), numel(xnewfine));
clear ebathmasknew;
for ii = 1:size(ebathmask, 2)
    MM = ebathmask(:, ii);
    MM(-zvec > wtopotransect(ii)) = ones;  % Mark water column
    ebathmasknew(:, ii) = MM;
end

%% Fill NaN values using specialized interpolation
% Fill using nearest neighbor row interpolation
sigma_gridded_f = fillNaNsrow(sigma_gridded);
temp_gridded_f = fillNaNsrow(temp_gridded);
salt_gridded_f = fillNaNsrow(salt_gridded);

% Fill velocity using decay-based interpolation
v_gridded_f = fillNaNsWithDecay(v_gridded, 60, 100);

% Apply smoothing and mask below seafloor
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
xticklabels([])
ylabel('Depth (m)')
title('Stat. 105-116, \Theta')
ylim([-200 30])
xlim([0 50])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty
plot(xx_moor, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

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
xticklabels([])
ylabel('Depth (m)')
title('S_A')
ylim([-200 30])
xlim([0 50])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty
plot(xx_moor, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

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
ylabel('Depth (m)')
title('Velocity (m/s)')
ylim([-200 30])
xlim([0 50])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, .3, 0.96])
makepretty
plot(xx_moor, 30, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

%% Calculate freshwater transport in salinity space
clear integrand_interp filtered_integrand_interp

% Run external scripts for freshwater flux calculations
Sflux;       % Main freshwater flux calculation
Sfluxalt;    % Alternative calculation for south transect

% Calculate total freshwater transport
dSbins = S_bins(2) - S_bins(1);
dx = (xnewfine(2) - xnewfine(1)) * 1000;  % Convert km to m
integral = sum(sum(integrand_interp(:,:))) * dSbins * dx;

% Plot freshwater transport in salinity space
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

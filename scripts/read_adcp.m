% read_adcp.m - Load and process shipboard ADCP data for OSNAP transects
%
% DESCRIPTION:
%   This script loads shipboard ADCP (Acoustic Doppler Current Profiler) data
%   for OSNAP East and West transects, calculates velocity magnitudes, and
%   prepares smoothed velocity profiles for use in transect analysis.
%
% USAGE:
%   Set the 'year' variable to select data year (2014, 2016, 2018, 2020, 2022)
%   Ensure 'eastlon_orig', 'eastlat_orig', 'westlon_orig', 'westlat_orig' 
%   are defined before running this script.
%
% REQUIRED INPUTS (must be defined before running):
%   - year: Data year (e.g., 2014, 2016, 2018, 2020, 2022)
%   - eastlon_orig, eastlat_orig: OSNAP East transect origin coordinates
%   - westlon_orig, westlat_orig: OSNAP West transect origin coordinates
%
% OUTPUTS:
%   - newxe_adcp: Horizontal coordinates for OSNAP East ADCP (5 km grid)
%   - vunderwaylo: Smoothed velocity magnitude for OSNAP East (upper 100m)
%   - newxw_adcp: Horizontal coordinates for OSNAP West ADCP (5 km grid)
%   - vwunderwaylo: Smoothed velocity magnitude for OSNAP West (upper 100m)
%
% DATA STRUCTURE:
%   ADCP data are organized by year with the following naming conventions:
%   2014: sect1os_os_dt.mat (East), sect4os_os_dt.mat (West)
%   2016: sect1os_nb_dt.mat (East), sect4os_nb_dt.mat (West)
%   2018: sect1os_nb_dt.mat (East), sect3os_nb_dt.mat (West)
%   2020: sect1_os150nb_dt.mat (East), sect4_on_os150nb_dt.mat (West)
%   2022: ar693_CFsec2_os150_dt_v2.mat (East), ar693_LSsec1_os150_dt_v2.mat (West)
%
% Author: [Your name]
% Date: [Date]

%% ========================================================================
%  USER CONFIGURATION - SET YEAR HERE
%  ========================================================================
year = 2016;  % Options: 2014, 2016, 2018, 2020, 2022

%% ========================================================================
%  CONFIGURE DATA PATHS AND FILENAMES
%  ========================================================================
% Base directory for ADCP data
% NOTE: Update this path to match your local directory structure
base_path = '/Users/zhazorken/Desktop/SIO_2024/Work/cruise_data/';

% Define year-specific filenames for East and West sections
switch year
    case 2014
        east_file = 'sect1os_os_dt.mat';
        west_file = 'sect4os_os_dt.mat';
        
    case 2016
        east_file = 'sect1os_nb_dt.mat';
        west_file = 'sect4os_nb_dt.mat';
        
    case 2018
        east_file = 'sect1os_nb_dt.mat';
        west_file = 'sect3os_nb_dt.mat';
        
    case 2020
        east_file = 'sect1_os150nb_dt.mat';
        west_file = 'sect4_on_os150nb_dt.mat';
        
    case 2022
        east_file = 'ar693_CFsec2_os150_dt_v2.mat';
        west_file = 'ar693_LSsec1_os150_dt_v2.mat';
        
    otherwise
        error('Invalid year. Choose from: 2014, 2016, 2018, 2020, 2022');
end

% Construct full file paths
east_path = fullfile(base_path, num2str(year), 'vmadcp', 'sections', east_file);
west_path = fullfile(base_path, num2str(year), 'vmadcp', 'sections', west_file);

fprintf('Processing ADCP data for year %d\n', year);
fprintf('  East section: %s\n', east_file);
fprintf('  West section: %s\n', west_file);

%% ========================================================================
%  PROCESS OSNAP EAST SECTION
%  ========================================================================
fprintf('\nLoading OSNAP East ADCP data...\n');
load(east_path, 'vm_data');
vm_stat = vm_data;

% Calculate along-transect distances
vmlon = vm_data.lon;
vmlat = vm_data.lat;

clear xx_eastadcp
for ii = 1:size(vmlon, 2)
    xx_eastadcp(ii) = m_lldist([vm_data.lon(ii) eastlon_orig], ...
                                [vm_data.lat(ii) eastlat_orig]);
end
xx_eastadcp = sort(xx_eastadcp);

% Calculate velocity magnitude from u and v components
% Flip data for proper orientation and compute magnitude with sign from v
u_mag = sqrt((fliplr(vm_stat.u)).^2 + (fliplr(vm_stat.v)).^2) .* ...
        (fliplr(vm_stat.v)) ./ abs(fliplr(vm_stat.v));

% For detided data (if available)
if isfield(vm_stat, 'u_dt') && isfield(vm_stat, 'v_dt')
    udt_mag = sqrt((fliplr(vm_stat.u_dt)).^2 + (fliplr(vm_stat.v_dt)).^2) .* ...
              (fliplr(vm_stat.v_dt)) ./ abs(fliplr(vm_stat.v_dt));
    use_detided = true;
    fprintf('  Using detided velocity data\n');
else
    udt_mag = u_mag;
    use_detided = false;
    fprintf('  Detided data not available, using raw velocities\n');
end

% Average over upper 100m (assuming 4m vertical resolution, 25 bins)
% NOTE: Adjust if your vertical resolution differs
upper_100m_bins = 1:25;
udt_mean_upper = mean(udt_mag(upper_100m_bins, :), 'omitnan');

% Interpolate to 5 km horizontal grid
newxe_adcp = xx_eastadcp(1):5:xx_eastadcp(end-1);
udt_magnew = interp1(xx_eastadcp, udt_mean_upper, newxe_adcp);

% Smooth the velocity profile
vunderwaylo = smooth(udt_magnew);

fprintf('  OSNAP East processed:\n');
fprintf('    Distance range: %.1f - %.1f km\n', xx_eastadcp(1), xx_eastadcp(end));
fprintf('    Grid points: %d (5 km spacing)\n', length(newxe_adcp));

%% ========================================================================
%  PROCESS OSNAP WEST SECTION
%  ========================================================================
fprintf('\nLoading OSNAP West ADCP data...\n');
load(west_path, 'vm_data');
vm_statw = vm_data;

% Calculate along-transect distances
vmlon = vm_data.lon;
vmlat = vm_data.lat;

clear xx_westadcp
for ii = 1:size(vmlon, 2)
    xx_westadcp(ii) = m_lldist([vm_data.lon(ii) westlon_orig], ...
                                [vm_data.lat(ii) westlat_orig]);
end

% Calculate velocity magnitude from u and v components
u_mag = sqrt((fliplr(vm_statw.u)).^2 + (fliplr(vm_statw.v)).^2) .* ...
        (fliplr(vm_statw.v)) ./ abs(fliplr(vm_statw.v));

% For detided data (if available)
if isfield(vm_statw, 'u_dt') && isfield(vm_statw, 'v_dt')
    udt_mag = sqrt((fliplr(vm_statw.u_dt)).^2 + (fliplr(vm_statw.v_dt)).^2) .* ...
              (fliplr(vm_statw.v_dt)) ./ abs(fliplr(vm_statw.v_dt));
    use_detided_west = true;
    fprintf('  Using detided velocity data\n');
else
    udt_mag = u_mag;
    use_detided_west = false;
    fprintf('  Detided data not available, using raw velocities\n');
end

% Average over upper 100m (assuming 4m vertical resolution, 25 bins)
upper_100m_bins = 1:25;
u_mean_upper = mean(udt_mag(upper_100m_bins, :), 'omitnan');

% Optional offset adjustment (set to 0 by default)
OFFSET = 0;

% Interpolate to 5 km horizontal grid
newxw_adcp = xx_westadcp(1):5:xx_westadcp(end-1);
u_magnew = interp1(xx_westadcp + OFFSET, u_mean_upper, newxw_adcp);

% Smooth the velocity profile
vwunderwaylo = smooth(u_magnew);

fprintf('  OSNAP West processed:\n');
fprintf('    Distance range: %.1f - %.1f km\n', xx_westadcp(1), xx_westadcp(end));
fprintf('    Grid points: %d (5 km spacing)\n', length(newxw_adcp));

%% ========================================================================
%  SUMMARY
%  ========================================================================
fprintf('\n========================================\n');
fprintf('ADCP Processing Complete for Year %d\n', year);
fprintf('========================================\n');
fprintf('Output variables:\n');
fprintf('  OSNAP East:\n');
fprintf('    newxe_adcp     : Horizontal coordinates [km]\n');
fprintf('    vunderwaylo    : Smoothed velocity magnitude [m/s]\n');
fprintf('  OSNAP West:\n');
fprintf('    newxw_adcp     : Horizontal coordinates [km]\n');
fprintf('    vwunderwaylo   : Smoothed velocity magnitude [m/s]\n');
fprintf('========================================\n');

%% ========================================================================
%  OPTIONAL: DIAGNOSTIC PLOTS
%  ========================================================================
% Uncomment the following section to generate diagnostic plots

% figure('Name', sprintf('ADCP Data %d', year), 'Position', [100 100 1200 600]);
% 
% % OSNAP East velocity profile
% subplot(2,1,1)
% plot(newxe_adcp, vunderwaylo, 'b-', 'LineWidth', 2);
% grid on;
% xlabel('Distance along transect (km)');
% ylabel('Velocity magnitude (m/s)');
% title(sprintf('OSNAP East - Upper 100m Mean Velocity (%d)', year));
% 
% % OSNAP West velocity profile
% subplot(2,1,2)
% plot(newxw_adcp, vwunderwaylo, 'r-', 'LineWidth', 2);
% grid on;
% xlabel('Distance along transect (km)');
% ylabel('Velocity magnitude (m/s)');
% title(sprintf('OSNAP West - Upper 100m Mean Velocity (%d)', year));

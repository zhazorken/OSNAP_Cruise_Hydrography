% hydromean.m - Calculate multi-year mean hydrographic fields (2014-2022)
%
% DESCRIPTION:
%   This script computes the mean Temperature (T), Salinity (S), density (sigma),
%   and velocity (v) fields from gridded hydrographic data spanning 2014-2022.
%   The script loads individual year data files, averages them, and saves the
%   result as a mean gridded product.
%
% REQUIRED INPUT FILES:
%   - 2014grid.mat
%   - 2016grid.mat
%   - 2018grid.mat
%   - 2020grid.mat
%   - 2022grid.mat
%
% OUTPUT FILE:
%   - meangrid.mat - Contains mean fields and grid parameters
%
% GRID PARAMETERS:
%   - Vertical: 2m resolution, starts at 1500m depth
%   - Horizontal: 0.25 km resolution, extends to 90 km
%


%% Set grid parameters
zbegind = 1500;  % Starting depth index (to focus on upper ocean)
xend = 90 * 4;   % Horizontal extent (90 km at 0.25 km resolution)

% Define file list for multi-year average
filelist = {'2014grid.mat', '2016grid.mat', '2018grid.mat', ...
            '2020grid.mat', '2022grid.mat'};

%% Load and extract data from each year
for jj = 1:size(filelist, 2)
    load(filelist{jj})
    fprintf('Loading file %d of %d: %s\n', jj, length(filelist), filelist{jj});
    
    % Extract bathymetry masks (subset to depth range and horizontal extent)
    wbathmasknewmat{jj} = wbathmasknew(end-zbegind:end, :);
    ebathmasknewmat{jj} = ebathmasknew(end-zbegind:end, :);
    
    % Extract Temperature fields
    Twestmat{jj} = Twest(end-zbegind:end, 1:xend);
    Teastmat{jj} = Teast(end-zbegind:end, 1:xend);
    
    % Extract Salinity fields
    Swestmat{jj} = Swest(end-zbegind:end, 1:xend);
    Seastmat{jj} = Seast(end-zbegind:end, 1:xend);
    
    % Extract density (sigma) fields
    sigwestmat{jj} = sigwest(end-zbegind:end, 1:xend);
    sigeastmat{jj} = sigeast(end-zbegind:end, 1:xend);
    
    % Extract velocity fields (note: dimensions are transposed)
    vwestmat{jj} = vwest(1:xend, end-zbegind:end);
    veastmat{jj} = veast(1:xend, end-zbegind:end);
end

%% Compute multi-year means
% Average individual 2D T/S/sigma/v fields over all years
Twest = (Twestmat{1} + Twestmat{2} + Twestmat{3} + Twestmat{4} + Twestmat{5}) / 5;
Teast = (Teastmat{1} + Teastmat{2} + Teastmat{3} + Teastmat{4} + Teastmat{5}) / 5;
Swest = (Swestmat{1} + Swestmat{2} + Swestmat{3} + Swestmat{4} + Swestmat{5}) / 5;
Seast = (Seastmat{1} + Seastmat{2} + Seastmat{3} + Seastmat{4} + Seastmat{5}) / 5;
sigwest = (sigwestmat{1} + sigwestmat{2} + sigwestmat{3} + sigwestmat{4} + sigwestmat{5}) / 5;
sigeast = (sigeastmat{1} + sigeastmat{2} + sigeastmat{3} + sigeastmat{4} + sigeastmat{5}) / 5;
vwest = (vwestmat{1} + vwestmat{2} + vwestmat{3} + vwestmat{4} + vwestmat{5}) / 5;
veast = (veastmat{1} + veastmat{2} + veastmat{3} + veastmat{4} + veastmat{5}) / 5;

%% Subset bathymetry masks to match data extent
wbathmasknew = wbathmasknew(end-zbegind:end, 1:xend);
ebathmasknew = ebathmasknew(end-zbegind:end, 1:xend);

%% Define grid parameters for output
xres = 0.25;  % Horizontal resolution in km

% Depth vector (subset to focus range)
zvec = zvec(end-zbegind:end);

% Horizontal coordinate vectors (in km)
xnew = 0.25:0.25:xend/4;
xneww = 0.25:0.25:xend/4;
xnewe = 0.25:0.25:xend/4;
xneww_v = 0.25:0.25:xend/4;
xnewe_v = 0.25:0.25:xend/4;

% Year label for output
year = 'mean';

%% Save mean gridded product
save('meangrid')

fprintf('\nMean grid saved to meangrid.mat\n');
fprintf('Grid dimensions:\n');
fprintf('  Depth: %d levels\n', length(zvec));
fprintf('  Horizontal: %d points (%.2f km)\n', length(xnew), max(xnew));

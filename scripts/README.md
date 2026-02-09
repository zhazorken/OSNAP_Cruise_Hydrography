# OSNAP Hydrographic Gridded Product Scripts (2014-2022)

## Overview
This repository contains MATLAB scripts for processing CTD (Conductivity-Temperature-Depth) and ADCP (Acoustic Doppler Current Profiler) data to create gridded hydrographic products along ocean transects. The scripts were developed for the OSNAP (Overturning in the Subpolar North Atlantic Program) dataset covering 2014-2022.

## Citation
If you use these scripts, please cite:
Zhao, Straneo, Holte, Myers, Torres et al. 2026 (JPO)

## Repository Contents

### Main Processing Scripts
- **maketransect.m** - Main script for processing individual cruise data and creating gridded transect products
- **maketransect_from2014_2022mean.m** - Script for plotting and analyzing the multi-year mean gridded product
- **hydromean.m** - Computes multi-year mean fields from individual year data files

### Utility Functions
- **fillNaNs.m** - Fill NaN values using nearest neighbor interpolation
- **fillNaNsrow.m** - Fill NaN values in rows using linear interpolation
- **fillNaNsWithDecay.m** - Fill NaN values using exponentially-decayed interpolation
- **fillRemainingNaNs.m** - Fill remaining NaN values by averaging adjacent cells
- **findClosestElement.m** - Find closest grid elements to specified lat/lon coordinates

## Data Requirements

### Input Data
The scripts require the following data types:

1. **CTD Data** (available via CCHDO repository)
   - Temperature (in-situ and potential)
   - Salinity (practical and absolute)
   - Pressure/Depth
   - Geographic coordinates (latitude, longitude)
   - Access: https://cchdo.ucsd.edu/cruise/33VB20220819

2. **ADCP Data** (shipboard velocity measurements)
   - U and V velocity components
   - Depth levels
   - Location: "vmadcp" folder in the data repository

3. **Bathymetry Data**
   - Gridded topography/bathymetry
   - Latitude and longitude grids

### Data Format
The scripts expect data to be pre-loaded into MATLAB workspace with the following variable names:
- `tempmat` - Temperature matrix [depth × station]
- `salmat` - Salinity matrix [depth × station]
- `depthmat` - Depth matrix [depth × station]
- `lonmat` - Longitude vector [1 × station]
- `latmat` - Latitude vector [1 × station]
- `u2`, `v2` - Velocity components [depth × station]
- `depthadcp2` - ADCP depth levels

## Installation and Setup

### Required MATLAB Toolboxes and Functions
1. **GSW Oceanographic Toolbox** (v3.06.16 or later)
   - Download from: http://www.teos-10.org/software.htm
   - Used for: Absolute salinity, potential temperature, and density calculations

2. **M_Map** (v1.4 or later)
   - Download from: https://www.eoas.ubc.ca/~rich/map.html
   - Used for: Geographic distance calculations

3. **cmocean** colormaps
   - Download from: https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean
   - Used for: Oceanographic color schemes

4. **Additional Functions**
   - `smooth2a.m` - 2D smoothing function: https://www.mathworks.com/matlabcentral/fileexchange/23287-smooth2a
   - `makepretty.m` - Figure formatting function: Placeholder for improving figure formatting (feel free to replace with your own)
   - `pmkmp.m` - Perceptually uniform colormaps: https://www.mathworks.com/matlabcentral/fileexchange/28982-perceptually-improved-colormaps

### Path Configuration
Before running the scripts, update the `addpath` commands in `maketransect.m` and `maketransect_from2014_2022mean.m` to match your local directory structure:

```matlab
addpath /your/path/to/gsw_matlab_v3_06_16/
addpath /your/path/to/gsw_matlab_v3_06_16/library/
addpath /your/path/to/m_map/
addpath /your/path/to/custom_functions/
```

## Usage

### Workflow Overview
1. Load raw CTD and ADCP data into MATLAB workspace
2. Run `maketransect.m` to process individual cruise data
3. Repeat for each year (2014, 2016, 2018, 2020, 2022)
4. Run `hydromean.m` to compute multi-year mean
5. Run `maketransect_from2014_2022mean.m` to visualize mean product

### Processing Individual Cruises
```matlab
% 1. Load your CTD and ADCP data
load('your_ctd_data.mat')

% 2. Configure transect parameters in maketransect.m
% Edit lines 52-76 to select desired CTD stations

% 3. Run the script
maketransect

% 4. Save gridded output
save('2022grid.mat')  % Use appropriate year
```

### Creating Multi-Year Mean
```matlab
% 1. Ensure all individual year files exist:
%    - 2014grid.mat
%    - 2016grid.mat
%    - 2018grid.mat
%    - 2020grid.mat
%    - 2022grid.mat

% 2. Run the mean calculation
hydromean

% 3. Output: meangrid.mat
```

### Visualizing Mean Product
```matlab
% Load and plot the mean gridded product
maketransect_from2014_2022mean
```

## Grid Specifications

### Spatial Resolution
- **Horizontal**: 0.25 km (250 m) standard grid, 0.1 km (100 m) fine grid
- **Vertical**: 2 m resolution
- **Depth range**: Surface to 4000 m (subset to 1500 m for some analyses)

### Coordinate System
- **Along-transect distance**: Calculated using great circle distances (m_lldist)
- **Origin points**: 
  - OSNAP East: 60.1022°N, 43.1189°W
  - OSNAP West: 59.9854°N, 45.2654°W

## Output Variables

### Gridded Fields
- `temp_gridded_s` - Smoothed temperature [°C]
- `salt_gridded_s` - Smoothed absolute salinity [g/kg]
- `sigma_gridded_s` - Smoothed potential density [kg/m³]
- `v_gridded_s` - Smoothed cross-transect velocity [m/s]

### Grid Vectors
- `xnewfine` - Horizontal coordinates [km]
- `zvecfine` - Depth coordinates [m]

### Masks
- `ebathmasknew` - Bathymetry mask (1 = water, NaN = below seafloor)

## Key Processing Steps

### 1. Station Selection and Sorting
Stations are selected based on geographic criteria and sorted by distance along the transect.

### 2. Velocity Projection
ADCP velocities are projected onto the direction perpendicular to the transect using:
```matlab
angle_deg = -4.5383;  % Transect angle
perp_velocity = u*sin(angle) - v*cos(angle)
```

### 3. Grid Interpolation
Data are interpolated from irregular station spacing onto regular grids using MATLAB's `interp2` function with linear interpolation.

### 4. NaN Filling Strategy
Multiple passes are used to fill missing values:
1. `fillNaNsrow` - Nearest neighbor interpolation along rows
2. `fillNaNsWithDecay` - Exponentially-decayed interpolation (for velocity)
3. Smoothing using `smooth2a` (5×5 window)
4. Application of bathymetry mask

### 5. Transport Calculations
Freshwater and volume transports are calculated in salinity space using external scripts:
- `Sflux.m` - Freshwater flux calculation
- `transport.m` - Volume transport calculation

## Customization Options

### Transect Selection
Sometimes, we want to make our own transects or mix-match from other OSNAP Hydrography datasets. An example of this to make your own transect with manually selected stations is shown with three pre-configured transect options are available in `maketransect.m`:
- **Mid transect**: Stations 85-98
- **North/trough transect**: Stations 105-116 (default)
- **South transect**: Stations 65-75, 81-83

Edit lines 52-76 to select or define custom transects.

### Interpolation Parameters
Decay-based interpolation parameters can be adjusted:
```matlab
v_gridded_f = fillNaNsWithDecay(v_gridded, 60, 100);
%                                         ^    ^
%                                   row decay  column decay
```

### Smoothing Window
2D smoothing window size (default 5×5):
```matlab
sigma_gridded_s = fillNaNsrow(smooth2a(sigma_gridded_f, 5, 5))
%                                                       ^  ^
```

## Troubleshooting

### Common Issues

1. **Missing dependencies**
   - Error: `Undefined function or variable 'gsw_SA_from_SP'`
   - Solution: Install GSW Oceanographic Toolbox and add to path

2. **Path errors**
   - Error: `Unable to read file...`
   - Solution: Update `addpath` commands to match your directory structure

3. **Memory issues with large datasets**
   - Solution: Process data in chunks or reduce grid resolution

4. **Missing data**
   - Some NaN values may remain after interpolation if isolated from valid data
   - Use `fillRemainingNaNs.m` for additional filling

## Version History
- v1.0 (2024) - Initial release for publication

## Contact and Support
For questions or issues, please contact:
kenzhao@unc.edu

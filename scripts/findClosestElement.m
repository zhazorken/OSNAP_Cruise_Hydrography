function [closestLatIndex, closestLatValue, closestLonIndex, closestLonValue] = ...
    findClosestElement(latgrid, latval, longrid, lonval)
% findClosestElement - Find closest grid elements to given lat/lon values
%
% SYNTAX:
%   [closestLatIndex, closestLatValue, closestLonIndex, closestLonValue] = ...
%       findClosestElement(latgrid, latval, longrid, lonval)
%
% DESCRIPTION:
%   Finds the indices and values of grid elements closest to specified
%   latitude and longitude values. Useful for matching observations to
%   gridded datasets.
%
% INPUTS:
%   latgrid - Vector of latitude grid points
%   latval  - Target latitude value
%   longrid - Vector of longitude grid points
%   lonval  - Target longitude value
%
% OUTPUTS:
%   closestLatIndex - Index of closest latitude in latgrid
%   closestLatValue - Value of closest latitude
%   closestLonIndex - Index of closest longitude in longrid
%   closestLonValue - Value of closest longitude
%
% EXAMPLE:
%   [iLat, lat, iLon, lon] = findClosestElement(latGrid, 60.5, lonGrid, -43.2);

    % Find the closest latitude element
    [~, closestLatIndex] = min(abs(latgrid - latval));
    closestLatValue = latgrid(closestLatIndex);
    
    % Find the closest longitude element
    [~, closestLonIndex] = min(abs(longrid - lonval));
    closestLonValue = longrid(closestLonIndex);
end

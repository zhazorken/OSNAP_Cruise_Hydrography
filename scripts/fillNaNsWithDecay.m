function filledMatrix = fillNaNsWithDecay(matrix, rowDecay, colDecay)
% fillNaNsWithDecay - Fill NaN values using exponentially decayed nearest values
%
% SYNTAX:
%   filledMatrix = fillNaNsWithDecay(matrix, rowDecay, colDecay)
%
% DESCRIPTION:
%   Fills NaN values using an exponential decay function based on distance
%   to nearest non-NaN values in both row and column directions. The decay
%   factors control how quickly influence decreases with distance.
%
% INPUTS:
%   matrix   - Input matrix with NaN values to fill
%   rowDecay - Decay length scale for row direction (higher = slower decay)
%   colDecay - Decay length scale for column direction (higher = slower decay)
%
% OUTPUTS:
%   filledMatrix - Matrix with NaN values filled using decayed interpolation
%
% EXAMPLE:
%   A = fillNaNsWithDecay(myMatrix, 60, 100);
%   % Uses decay lengths of 60 for rows and 100 for columns

    filledMatrix = matrix;
    [rows, cols] = size(matrix);
    
    % Iterate over each element
    for i = 1:rows
        for j = 1:cols
            if isnan(matrix(i, j))
                % Find the nearest non-NaN value in the row and column
                [nearestRowVal, nearestRowDist] = findNearestNonNaN(matrix(i, :), j);
                [nearestColVal, nearestColDist] = findNearestNonNaN(matrix(:, j), i);
                
                % Calculate decay factors based on distance
                % Decay = exp(-distance / decay_length)
                rowDecayFactor = exp(-nearestRowDist / rowDecay);
                colDecayFactor = exp(-nearestColDist / colDecay);
                
                % Use the decayed values to fill in the NaN
                if isnan(nearestRowVal)
                    % Only column value available
                    filledMatrix(i, j) = nearestColVal * colDecayFactor;
                elseif isnan(nearestColVal)
                    % Only row value available
                    filledMatrix(i, j) = nearestRowVal * rowDecayFactor;
                else
                    % Average the decayed values if both are available
                    filledMatrix(i, j) = (nearestRowVal * rowDecayFactor + ...
                                         nearestColVal * colDecayFactor) / 2;
                end
            end
        end
    end
end

function [nearestVal, nearestDist] = findNearestNonNaN(vec, idx)
% findNearestNonNaN - Find nearest non-NaN value and its distance in a vector
%
% INPUTS:
%   vec - Vector to search (row or column)
%   idx - Current index position
%
% OUTPUTS:
%   nearestVal  - Value of nearest non-NaN element
%   nearestDist - Distance to nearest non-NaN element

    % Find the non-NaN indices
    nonNaNIndices = find(~isnan(vec));
    
    if isempty(nonNaNIndices)
        nearestVal = NaN;
        nearestDist = Inf;
        return;
    end
    
    % Calculate distances to all non-NaN values
    distances = abs(nonNaNIndices - idx);
    [nearestDist, minIdx] = min(distances);
    nearestVal = vec(nonNaNIndices(minIdx));
end

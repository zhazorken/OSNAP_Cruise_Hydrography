function filledMatrix = fillNaNs(matrix)
% fillNaNs - Fill NaN values in a matrix using nearest neighbor interpolation
%
% SYNTAX:
%   filledMatrix = fillNaNs(matrix)
%
% DESCRIPTION:
%   Fills NaN values in each row by replacing them with the closest non-NaN
%   value (left or right). If equidistant, chooses the left value.
%
% INPUTS:
%   matrix - Input matrix with potential NaN values
%
% OUTPUTS:
%   filledMatrix - Matrix with NaN values filled
%
% EXAMPLE:
%   A = [1 NaN 3 NaN 5];
%   B = fillNaNs(A);  % Returns [1 1 3 3 5]

    filledMatrix = matrix;
    [rows, cols] = size(matrix);
    
    % Iterate over each row
    for i = 1:rows
        row = matrix(i, :);
        
        % Find the indices of non-NaN values
        nonNaNIndices = find(~isnan(row));
        
        % Iterate over each NaN value in the row
        nanIndices = find(isnan(row));
        for j = 1:length(nanIndices)
            nanIndex = nanIndices(j);
            
            % Find the closest non-NaN value to the left
            leftNonNaN = findClosestNonNaN(row, nonNaNIndices, nanIndex, true);
            
            % Find the closest non-NaN value to the right
            rightNonNaN = findClosestNonNaN(row, nonNaNIndices, nanIndex, false);
            
            % Choose the closest non-NaN value
            if isnan(leftNonNaN)
                closestNonNaN = rightNonNaN;
            elseif isnan(rightNonNaN)
                closestNonNaN = leftNonNaN;
            else
                % If both exist, choose the closer one
                if abs(nanIndex - nonNaNIndices(leftNonNaN)) < abs(nanIndex - nonNaNIndices(rightNonNaN))
                    closestNonNaN = leftNonNaN;
                else
                    closestNonNaN = rightNonNaN;
                end
            end
            
            % Replace the NaN value with the closest non-NaN value
            if ~isnan(closestNonNaN)
                filledMatrix(i, nanIndex) = row(nonNaNIndices(closestNonNaN));
            end
        end
    end
end

function closestNonNaNIndex = findClosestNonNaN(row, nonNaNIndices, nanIndex, searchLeft)
% findClosestNonNaN - Helper function to find closest non-NaN value index
%
% INPUTS:
%   row            - Current row being processed
%   nonNaNIndices  - Indices of non-NaN values in the row
%   nanIndex       - Index of the NaN value being filled
%   searchLeft     - If true, search left; if false, search right
%
% OUTPUTS:
%   closestNonNaNIndex - Index into nonNaNIndices array of closest value

    closestNonNaNIndex = NaN;
    if searchLeft
        % Search from right to left for the first non-NaN to the left
        for k = length(nonNaNIndices):-1:1
            if nonNaNIndices(k) < nanIndex
                closestNonNaNIndex = k;
                break;
            end
        end
    else
        % Search from left to right for the first non-NaN to the right
        for k = 1:length(nonNaNIndices)
            if nonNaNIndices(k) > nanIndex
                closestNonNaNIndex = k;
                break;
            end
        end
    end
end

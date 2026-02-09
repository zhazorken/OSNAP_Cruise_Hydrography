function filledMatrix = fillRemainingNaNs(matrix)
% fillRemainingNaNs - Fill remaining NaN values using adjacent cell averaging
%
% SYNTAX:
%   filledMatrix = fillRemainingNaNs(matrix)
%
% DESCRIPTION:
%   Fills NaN values by averaging the non-NaN values in the 4 adjacent
%   cells (left, right, above, below). Useful as a final pass after other
%   interpolation methods.
%
% INPUTS:
%   matrix - Input matrix with potential NaN values
%
% OUTPUTS:
%   filledMatrix - Matrix with NaN values filled using local averaging
%
% NOTE:
%   NaN values with no non-NaN neighbors remain as NaN

    filledMatrix = matrix;
    [rows, cols] = size(matrix);
    
    for i = 1:rows
        for j = 1:cols
            if isnan(filledMatrix(i, j))
                % Collect surrounding non-NaN values
                surroundingVals = [];
                
                % Check left neighbor
                if j > 1 && ~isnan(filledMatrix(i, j-1))
                    surroundingVals = [surroundingVals filledMatrix(i, j-1)];
                end
                
                % Check right neighbor
                if j < cols && ~isnan(filledMatrix(i, j+1))
                    surroundingVals = [surroundingVals filledMatrix(i, j+1)];
                end
                
                % Check above neighbor
                if i > 1 && ~isnan(filledMatrix(i-1, j))
                    surroundingVals = [surroundingVals filledMatrix(i-1, j)];
                end
                
                % Check below neighbor
                if i < rows && ~isnan(filledMatrix(i+1, j))
                    surroundingVals = [surroundingVals filledMatrix(i+1, j)];
                end
                
                % Fill with average of surrounding values if any exist
                if ~isempty(surroundingVals)
                    filledMatrix(i, j) = mean(surroundingVals);
                end
            end
        end
    end
end

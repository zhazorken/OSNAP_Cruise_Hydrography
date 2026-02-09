function filledMatrix = fillNaNsrow(matrix)
% fillNaNsrow - Fill NaN values in each row using nearest neighbor interpolation
%
% SYNTAX:
%   filledMatrix = fillNaNsrow(matrix)
%
% DESCRIPTION:
%   Fills NaN values in each row using MATLAB's interp1 with 'nearest'
%   method and extrapolation. Requires at least 2 non-NaN values per row.
%
% INPUTS:
%   matrix - Input matrix with potential NaN values
%
% OUTPUTS:
%   filledMatrix - Matrix with NaN values filled via interpolation
%
% NOTE:
%   Rows with fewer than 2 non-NaN values remain unchanged

    filledMatrix = matrix;
    [rows, cols] = size(matrix);
    
    % Iterate over each row
    for i = 1:rows
        row = matrix(i, :);
        
        % Find the indices of non-NaN values
        nonNaNIndices = find(~isnan(row));
        
        % If the row contains NaN values, interpolate to fill them
        if ~isempty(nonNaNIndices) && length(nonNaNIndices) < cols
            
            % Perform interpolation if there are at least two non-NaN values
            if numel(nonNaNIndices) >= 2
                % Use nearest neighbor interpolation with extrapolation
                filledMatrix(i, :) = interp1(nonNaNIndices, row(nonNaNIndices), ...
                                             1:cols, 'nearest', 'extrap');
            end
        end
    end
end

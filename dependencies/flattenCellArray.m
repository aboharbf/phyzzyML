function flattenedArray = flattenCellArray(cellArray)
    flattenedArray = {};
    for i = 1:numel(cellArray)
        if iscell(cellArray{i})
            % Recursively flatten nested cell arrays
            flattenedSubArray = flattenCellArray(cellArray{i});
            flattenedArray = [flattenedArray, flattenedSubArray];
        else
            % Append non-cell elements to the flattened array
            flattenedArray = [flattenedArray, cellArray{i}];
        end
    end
end
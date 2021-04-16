function maxStetchLength = findStretches(dataMat)

% Pad the sides to make sure each stretch begins and ends.
dataMat = [zeros(size(dataMat,1), 1), dataMat, zeros(size(dataMat,1), 1)];

% diffVec = [ones(size(vector,1), 1), diff(vector, 1, 2), -1 * ones(size(vector,1), 1)];
diffVec = diff(dataMat, 1, 2);
maxStetchLength = zeros(size(diffVec,1), 1);

for unit_i = 1:length(maxStetchLength)
  % Look for starts and stops
  rowData = diffVec(unit_i, :);
  starts = find(rowData == 1);
  stops = find(rowData == -1);
  lengths = stops - starts;
  
  if ~isempty(lengths)
    maxStetchLength(unit_i) = max(lengths);
  end
  
end

end
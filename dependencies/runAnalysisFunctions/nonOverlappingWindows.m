function elements2Keep = nonOverlappingWindows(startStopInds, priorityRank)
% A function which accepts a set of inputs and outputs, and determines if
% there is overlap. 
% Inputs
% - startStopInds: a n*2 matrix with start and stop indicies (positive
% ints)
% - priorityRank: an n*1 vector showing the current primacy of the arrays.
% In instances of overlap, the lower ranking element is removed.
% Outputs:
% - element2keep = an n*1 vector of which elements to keep to avoid any
% overlaps while respecting rank.

% Procedure takes place from top to bottom (checking and removing for the
% highest rank object, then the second, etc).

% convert priority ranks into a sorting
[rankSorted, newInd] = sort(priorityRank, 'descend');

startStopIndsSorted = startStopInds(newInd, :);

% Create an array
tmpArray = false(length(rankSorted), max(startStopIndsSorted(:)));

for row_i = 1:size(startStopIndsSorted,1)
  tmpArray(row_i, startStopIndsSorted(row_i,1):startStopIndsSorted(row_i,2)) = 1;
end

tmpArrayOrig = tmpArray;

elements2KeepSorted = true(size(rankSorted));
for ii = 1:length(elements2KeepSorted)
  if ~elements2KeepSorted(ii)
    continue
  end
  
  % Grab the row, and see if any in the other rows match it
  rowInd = false(size(elements2KeepSorted));
  rowInd(ii) = true;
  
  priorityRow = tmpArray(rowInd, :);
  otherRows = tmpArray(~rowInd, priorityRow);
  otherRowsWithElements = any(otherRows,2);
  otherRowNum = find(~rowInd);
  
  otherRowNum2Clear = otherRowNum(otherRowsWithElements)';
  % Empty out those rows with something in them
  for row_i = otherRowNum2Clear
    tmpArray(row_i,:) = false(size(tmpArray(row_i,:)));
  end
end

elements2KeepSorted = any(tmpArray,2);

% Plotting
if 0
  figure();
  
  subplot(1,5,1:2);
  imagesc(tmpArrayOrig);
  title('All Stretches');
  
  subplot(1,5,3);
  imagesc(priorityRank);
  title('R^2 Sorted');
  
  subplot(1,5,4:5);
  imagesc(tmpArray);
  title(sprintf('All Stretches Left after exclusive (%d)', sum(elements2KeepSorted)))
end

% Unsort
elements2Keep(newInd) = elements2KeepSorted;

end
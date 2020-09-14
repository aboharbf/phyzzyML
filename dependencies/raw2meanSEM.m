function [meanStack, semStack, nStack] = raw2meanSEM(rawDataCellArray, grandMeanSwitch)
% Converts a cell array of raw values into 2 arrays, one with a mean for
% each cell, and the other with the SEM for each cell.


[meanStack, semStack] = deal(zeros(length(rawDataCellArray), size(rawDataCellArray{1},2)));
nStack = zeros(length(rawDataCellArray),1);

for cell_i = 1:length(rawDataCellArray)
  dataCell = rawDataCellArray{cell_i};
  meanStack(cell_i,:) = nanmean(dataCell);
  nStack(cell_i) = size(dataCell,1);
  semStack(cell_i,:) = nanstd(dataCell)/sqrt(nStack(cell_i));
end

if grandMeanSwitch
  grandStack = vertcat(rawDataCellArray{:});
  nStack = [nStack; sum(nStack)];
  meanStack = [meanStack; nanmean(grandStack)];
  semStack = [semStack; nanstd(grandStack)/sqrt(nStack(end))];
end 

end
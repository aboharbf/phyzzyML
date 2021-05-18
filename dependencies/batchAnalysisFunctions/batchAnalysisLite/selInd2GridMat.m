function [selIndGridMat, selIndGridMatLabels] = selInd2GridMat(selTable)
% a function which converts a selTable, with columns labeled selInd, and
% turns it into a matrix.
% Input:
% - selTable, with columns ending in '_selInd'
% Output:
% - selIndGridMat: an ML * AP * D Matrix spanning all possible positions, represented in the
% selTable. the 3rd dimension is the number of '_selInd' variables in the
% table.
% - labelArray: the labels corresponding to the 3rd Dimension of the selIndGridMat.

% Recover the events represented
labelArray = selTable.Properties.VariableNames(contains(selTable.Properties.VariableNames, '_selInd'));

% Make Grid hole arrays - identify all the possible gridHoles.
gridHoles = unique(selTable.gridHole);
gridHoleA = unique(string(extractBefore(gridHoles, 2)));
gridHoleB = string(unique(double(string(extractAfter(gridHoles, 1)))));
% [gridRuns, gridChannels, gridUnits]  = deal(zeros(length(gridHoleA), length(gridHoleB)));

selIndGridMat = zeros(length(gridHoleA), length(gridHoleB), length(labelArray) + 1, 2);

selIndexMat = selTable{:, labelArray};

% Add a final item for All Units.
labelArray = [labelArray, 'All Units'];
selIndexMat = [selIndexMat, true(size(selIndexMat,1),1)];

% Extract the gridHole elements
gridHoleMLVec = cellfun(@(x) x(1), selTable{:, 'gridHole'});
gridHoleAPVec = cellfun(@(x) x(2:end), selTable{:, 'gridHole'}, 'UniformOutput', false);

% Change the elements to indices
[~, gridHoleMLind] = ismember(gridHoleMLVec, gridHoleA);
[~, gridHoleAPind] = ismember(gridHoleAPVec, gridHoleB);

% Get a unit count per grid hole
for grid_i = 1:length(gridHoles)
  % Identify spot to add
  aInd = find(strcmp(gridHoles{grid_i}(1), gridHoleA));
  bInd = find(strcmp(gridHoles{grid_i}(2:end), gridHoleB));
  
  % Cycle through units
  for unit_i = 1:2
    
    if unit_i == 1
      unitIndex = contains(selTable.unitType, 'MUA');
    else
      unitIndex = ~contains(selTable.unitType, 'MUA');
    end
    
    % Cycle through Events
    for event_i = 1:length(labelArray)
      selIndGridMat(aInd, bInd, event_i, unit_i) = sum(gridHoleMLind == aInd & gridHoleAPind == bInd & selIndexMat(:,event_i) & unitIndex);
    end
    
  end
end

selIndGridMatLabels = {gridHoleA, gridHoleB, labelArray, {'MUA', 'U&US'}};

function selectivityPerEyeEvent(selTable, paradigm, params)
% a function which plots selectivity for fixed events in the paradigm, as
% defined in the events variable below. Also if desired, include
% directional saccade selectivity (probably not what we want).

includeSaccade = false;
UnitTypes = params.UnitTypes;
UnitTypePlot = params.UnitTypePlot;
alpha = 0.05;                                         % Alpha that the runs were done at.
outDir = fullfile(params.outputDir, 'selectivityPerEvent');

columnsOfInterest = {'saccDir_selInd', 'subSel_saccades_selInd', 'subSel_pre-saccades_selInd'};
eventNames = {'Directional Saccade', 'Post Saccade', 'Pre Saccade'};

% Generate bar plots showing total counts in each category
for unitType_i = 1:length(UnitTypes)
  % Filter the table per unit
  unitInd = contains(selTable.unitType, UnitTypes{unitType_i});
  selTableParadigmUnit = selTable(unitInd, :);
  
  % Plot 1 - fixation/saccade
  unitCount = sum(unitInd);
  chanceUnitCount = round(unitCount * alpha);

  selectivityIndex = selTableParadigmUnit{:, columnsOfInterest};
  selectivityCounts = sum(selectivityIndex,1);
  
  % Plot
  figTitle = sprintf('%s activity selective for non-Stimulus Events during %s', UnitTypePlot{unitType_i}, paradigm);
  createBarPlotWithChanceLine(eventNames, sum(selectivityIndex), alpha*2, unitCount, figTitle, []);
  saveFigure(outDir, figTitle, [], params.figStruct, [])
  
end

end
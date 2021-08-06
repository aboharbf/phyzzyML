function selectivityPerEyeEvent(selTable, paradigm, params)
% a function which plots selectivity for fixed events in the paradigm, as
% defined in the events variable below. Also if desired, include
% directional saccade selectivity (probably not what we want).

UnitTypes = params.UnitTypes;
UnitTypePlot = params.UnitTypePlot;
alpha = 0.05;                                         % Alpha that the runs were done at.
outDir = fullfile(params.outputDir, 'selectivityPerEvent');

figureTitles = {'%s direction selective activity during task events during %s', '%s activity selective for non-Stimulus Events during %s', '%s Activity Overlap between directional and non-directional during %s'};
columnsOfInterest = {{'saccDir_fixOnset_selInd', 'saccDir_stimOnset_selInd', 'saccDir_stim_selInd', 'saccDir_all_selInd'}, ...
  {'subSel_saccades_selInd', 'subSel_pre-saccades_selInd'}, {'subSel_pre-saccades_selInd', 'saccDir_all_selInd'}};
eventNames = {{'Fixation Dot Onset', 'Stimulus Onset', 'Stimulus', 'All'}, {'Pre Saccade', 'Post Saccade'}, {'Non-Directional', 'Directional'}};

% Generate bar plots showing total counts in each category
for unitType_i = 1:length(UnitTypes)
  % Filter the table per unit
  unitInd = contains(selTable.unitType, UnitTypes{unitType_i});
  selTableParadigmUnit = selTable(unitInd, :);
  unitCount = sum(unitInd);
  chanceUnitCount = round(unitCount * alpha);
  
  for plotType_i = 1:length(figureTitles)
    % Grab columns of interest
    selectivityIndex = selTableParadigmUnit{:, columnsOfInterest{plotType_i}};
    selectivityCounts = sum(selectivityIndex,1);
    
    % Plot
    figTitle = sprintf(figureTitles{plotType_i}, UnitTypePlot{unitType_i}, paradigm);
    createBarPlotWithChanceLine(eventNames{plotType_i}, selectivityCounts, alpha*2, unitCount, figTitle, []);
    saveFigure(outDir, figTitle, [], params.figStruct, [])
    
    % Venn Diagram
    if plotType_i == 1 || 3
      if plotType_i == 1
        sigInd = selectivityIndex(:, [1 3 4]);
        eventLabels = eventNames{plotType_i}([1 3 4]);
      else
        sigInd = selectivityIndex;
        eventLabels = eventNames{plotType_i};
      end
      
      figTitle = sprintf(figureTitles{plotType_i}, UnitTypePlot{unitType_i}, paradigm);
      vennXExpanded(sigInd, figTitle, eventLabels)
      saveFigure(outDir, sprintf('VD - %s', figTitle), [], params.figStruct, [])
      
    end
    
  end  
end

end
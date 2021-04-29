function selectivityPerEventBarGraphs(selTable, paradigm, params)
% a function which plots selectivity for fixed events in the paradigm, as
% defined in the events variable below. Also if desired, include
% directional saccade selectivity (probably not what we want).

includeSaccade = false;
UnitTypes = params.UnitTypes;
UnitTypePlot = params.UnitTypePlot;
alpha = 0.05;                                         % Alpha that the runs were done at.

% Generate bar plots showing total counts in each category
for unitType_i = 1:length(UnitTypes)
  % Filter the table per unit
  unitInd = contains(selTable.unitType, UnitTypes{unitType_i});
  selTableParadigmUnit = selTable(unitInd, :);
  
  if ~isempty(selTableParadigmUnit)
    % Plot 1 - fixation/saccade
    unitCount = sum(unitInd);
    chanceUnitCount = round(sum(unitCount) * alpha);
    saccSelCount = sum(selTableParadigmUnit.saccSel ~= 0);
    fixSelCount = sum(selTableParadigmUnit.BaseVFix_PVal < alpha);
    
    % Reward processing - take the rewardAbsent vec, and use the other reward
    % paradigm to fill in parts where that paradigm wasn't used.
    rewardVec = selTableParadigmUnit.subSel_rewardAbsent;
    rewardAbsMissing = isnan(rewardVec);
    rewardVec(rewardAbsMissing) = selTableParadigmUnit.subSel_reward(rewardAbsMissing);
    rewardVec = sum(rewardVec ~= 0);
    
    % Collect Data, labels.
    if includeSaccade
      events = {'Directional Saccade', 'Fixation', 'Reward'};
      dataMat = [saccSelCount; fixSelCount; rewardVec];
    else
      events = {'Fixation', 'Reward'};
      dataMat = [fixSelCount; rewardVec];
    end
    
    % Plot
    figTitle = sprintf('%s activity selective for non-Stimulus Events during %s', UnitTypePlot{unitType_i}, paradigm);
    createBarPlotWithChanceLine(events, dataMat, alpha, unitCount, figTitle, [])
    saveFigure(params.outputDir, figTitle, [], params.figStruct, [])
    
  end
end

end
function selectivityPerEventBarGraphs(selTable, params)
% a function which plots selectivity for fixed events in the paradigm, as
% defined in the events variable below. 

alpha = 0.05;                                         % Alpha that the runs were done at.
outDir = fullfile(params.outputDir, 'selectivityPerEvent');

% Generate bar plots showing total counts in each category
% Plot 1 - Task Events
figTitle = sprintf('%s activity selective for Fix and Reward during %s', params.unitTag, params.paradigmTag);
plotTitle = sprintf('Task Event Unit Selectivity');

eventTableTag = {'fixation', 'stimOnset', 'stimOffset', 'rewardCombo'};
eventNames = {'Fixation', 'Stimulus Onset', 'Stimulus Offset', 'Reward'};

unitCount = size(selTable, 1);
eventCounts = sum(selTable{:, strcat('subSel_', eventTableTag, '_selInd')});

% Plot
createBarPlotWithChanceLine(eventNames, eventCounts, alpha*2, unitCount, figTitle, []);
title(plotTitle);
saveFigure(outDir, figTitle, [], params.figStruct, [])

% Plot 2 - stimEvents
figTitle = sprintf('%s activity selective for Fix and Reward during %s', params.unitTag, params.paradigmTag);
plotTitle = sprintf('Stimulus Event Unit Selectivity');

eventTableTag = {'headTurn_right', 'headTurn_left', 'eyeContact', 'bodyTurn'};
eventNames = {'Head turn, right', 'Head turn, left', 'Eye contact', 'Body turn'};

if strcmp(params.paradigmTag, 'headTurnCon')
  eventTableTag = eventTableTag(1:2);
  eventNames = eventNames(1:2);
end

eventCounts = sum(selTable{:, strcat('subSel_', eventTableTag, '_selInd')});

% Plot
createBarPlotWithChanceLine(eventNames, eventCounts, alpha*2, unitCount, plotTitle, []);
saveFigure(outDir, figTitle, [], params.figStruct, [])

end
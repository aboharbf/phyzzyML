function selectivityPerEventBarGraphs(selTable, params)
% a function which plots selectivity for fixed events in the paradigm, as
% defined in the events variable below. Also if desired, include
% directional saccade selectivity (probably not what we want).

alpha = 0.05;                                         % Alpha that the runs were done at.
outDir = fullfile(params.outputDir, 'selectivityPerEvent');

% Generate bar plots showing total counts in each category
% Plot 1 - fixation/saccade
unitCount = size(selTable, 1);
fixSelCount = sum(selTable.subSel_fixation_selInd);

% Reward processing - take the rewardAbsent vec, and use the other reward
% paradigm to fill in parts where that paradigm wasn't used.
rewardVec = selTable.subSel_rewardCombo_selInd;
rewardVec = sum(rewardVec);

% Collect Data, labels.
events = {'Fixation', 'Reward'};
dataMat = [fixSelCount; rewardVec];

% Plot
figTitle = sprintf('%s activity selective for Fix and Reward during %s', params.unitTag, params.paradigmTag);
createBarPlotWithChanceLine(events, dataMat, alpha*2, unitCount, figTitle, []);
saveFigure(outDir, figTitle, [], params.figStruct, [])

end
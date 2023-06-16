function selectivityPerEventBarGraphsParComp(spikeBank, params)
% a function which plots selectivity for fixed events in the paradigm, as
% defined in the events variable below. 

alpha = 0.05;                                         % Alpha that the runs were done at.
outDir = fullfile(params.outputDir, 'selectivityPerEvent');
paradigmList = spikeBank.paradigmName;
uniquePar = unique(paradigmList);

% Comparing paradigms across events below - Filter each table to allow for
% concatonation
eventNames = {'fixation', 'stimOnset', 'stimOffset', 'rewardCombo', 'headTurn_right', 'headTurn_left'};
eventNames = strcat('subSel_', eventNames, '_selInd');
selTable = spikeBank.selTable;
for tab_i = 1:length(selTable)
  selTable{tab_i} = selTable{tab_i}(contains(selTable{tab_i}.unitType, digitsPattern), eventNames);
end
eventPlotNames = {'Fixation', 'Stimulus Onset', 'Stimulus Offset', 'Reward', 'Head turn, right', 'Head turn, left'};

eventData = nan(length(eventPlotNames), 2);
for par_i = 1:2
  parInd = strcmp(paradigmList, uniquePar(par_i));
  selTableAll = vertcat(selTable{parInd});
  eventData(:,par_i) = sum(selTableAll{:,:});
  unitCount(par_i) = size(selTableAll,1);
end

% Generate bar plots showing total counts in each category
% Plot 1 - Task Events
plotTitle = sprintf('Task and Stimulus Event Unit Selectivity');

legend2Add = strcat({'Natural', 'Animated'}, ' (', string(unitCount), ')');

% Plot
createBarPlotWithChanceLine(eventPlotNames, eventData, alpha*2, unitCount, plotTitle, legend2Add);
title(plotTitle);
figH = gcf;
figH.Position = [0.2957    0.3257    0.5332    0.4507];

% Space out long label names. 
eventNames = {'Fixation', 'StimulusXOnset', 'StimulusXOffset', 'Reward', 'Head turn,Xright', 'Head turn,Xleft'};
eventNames = cellfun(@(x) strrep(x,'X','\newline'), eventNames,'UniformOutput',false);
a = gca;
a.XTickLabel = eventNames;

saveFigure(outDir, plotTitle, [], params.figStruct, [])

end
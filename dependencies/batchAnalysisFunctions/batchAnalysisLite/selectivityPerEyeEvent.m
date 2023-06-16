function selectivityPerEyeEvent(selTable, params)
% a function which plots selectivity for fixed events in the paradigm, as
% defined in the events variable below. Also if desired, include
% directional saccade selectivity (probably not what we want).

alpha = 0.05;                                         % Alpha that the runs were done at.
outDir = fullfile(params.outputDir, 'selectivityPerEvent');

figureTitles = {'%s traces selective for (pre-)saccade period during %s', '%s Activity Overlap between directional and non-directional during %s',...
  '%s trace selective for (pre-)saccade (nonStim) during %s'};

columnsOfInterest = {{'subSel_saccades_selInd', 'subSel_pre_saccades_selInd'}, {'subSel_pre_saccades_selInd', 'saccDir_all_selInd'},...
  {'subSel_pre_saccadesNonStim_selInd', 'subSel_saccadesNonStim_selInd'}};

eventNames = {{'Pre Saccade', 'Post Saccade'}, {'Non-Directional', 'Directional'}...
  {'Pre Saccade', 'Post Saccade'}};

% Generate bar plots showing total counts in each category
unitCount = size(selTable, 1);

for plotType_i = 3%1:length(figureTitles)
  % Grab columns of interest
  selectivityIndex = selTable{:, columnsOfInterest{plotType_i}};
  selectivityCounts = sum(selectivityIndex,1);
  
  % Plot
  figTitle = sprintf(figureTitles{plotType_i}, params.unitTag, params.paradigmTag);
  createBarPlotWithChanceLine(eventNames{plotType_i}, selectivityCounts, alpha*2, unitCount, figTitle, []);
  saveFigure(outDir, figTitle, [], params.figStruct, [])
  
  % Venn Diagram
  if plotType_i == 2 || 4
%     if plotType_i == 2
%       sigInd = selectivityIndex(:, [1 3 4]);
%       eventLabels = eventNames{plotType_i}([1 3 4]);
%     else
      sigInd = selectivityIndex;
      eventLabels = eventNames{plotType_i};
%     end
    
    figTitle = sprintf(figureTitles{plotType_i}, params.unitTag, params.paradigmTag);
    vennXExpanded(sigInd, figTitle, eventLabels)
    saveFigure(outDir, sprintf('VD - %s', figTitle), [], params.figStruct, [])
    
  end
  
end
end


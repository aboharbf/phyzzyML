function selectivityPerAnyAndHeadTurn(selTable, params)
% a function which plots selectivity for fixed events in the params.paradigmTag, as
% defined in the events variable below. Also if desired, include
% directional saccade selectivity (probably not what we want).

alpha = 0.05;                                         % Alpha that the runs were done at.
outDir = fullfile(params.outputDir, 'selectivityPerSubEvent');

% Generate bar plots showing total counts in each category
tableVars = selTable.Properties.VariableNames';
subEventSelVars = tableVars(contains(tableVars, 'subSel_') & contains(tableVars, 'selInd') & ~contains(tableVars, 'reward') & ~contains(tableVars, 'blink') & ~contains(tableVars, 'all'));
subSelInfo = selTable{:,subEventSelVars};
subEventPlotNames = strrep(extractBetween(subEventSelVars, 'subSel_', '_selInd'), '_', ' ');

% Generate counts of selectivity in other regards (broadCat, socVNonSoc)
% for creating Venn Diagrams
epochSelInd = contains(tableVars, 'epochSel') & contains(tableVars, '_selInd') & ~contains(tableVars, 'vBase') & contains(tableVars, 'any');
epochSelIndVars = tableVars(epochSelInd);

  % Count per event
  countPerSubEvent = sum(subSelInfo);
  unitCount = size(subSelInfo, 1);
  atLeastOne = sum(any(subSelInfo,2));

  % Plot bar graphs for individual events
  figTitle = sprintf('%s activity selective for Stimulus sub-events (%d/%d Unique)', params.unitTag, atLeastOne, unitCount);
  figH = createBarPlotWithChanceLine(subEventPlotNames, countPerSubEvent, 0, unitCount, figTitle, []);
  
  % Move the legend
  legendHand = findobj(figH.Children, 'Type', 'Legend');
  legendHand.Location = 'northeastoutside';
  
  % Save the figure
  saveFigure(fullfile(outDir,'Bargraphs'), sprintf('BG %s %s', params.paradigmTag, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])
  
  % Remove individual selectivity columns, keep the combo for Epoch
  % cross over selectivity
  event2KeepIndex = ~(contains(subEventSelVars, 'all') | ~contains(subEventSelVars, 'headTurn'))';
  subEventPlotNamesVenn = subEventPlotNames(event2KeepIndex)';
  subSelInfo = subSelInfo(:,event2KeepIndex);
  epochSelIndPlotVars = extractBetween(epochSelIndVars, 'epochSel_', '_');
  
  
  for epoch_sel_i = 1:length(epochSelIndVars)
    epochSelInd = selTable.(epochSelIndVars{epoch_sel_i});
    
    subEventUnitEpoch = [subSelInfo, epochSelInd];
    
    colNames = [subEventPlotNamesVenn epochSelIndPlotVars(epoch_sel_i)];
    figTitle = sprintf('%s Activity with Head Turning and %s selectivity (%d/%d unique)', params.unitTag, epochSelIndPlotVars{epoch_sel_i}, sum(any(subEventUnitEpoch,2)), length(any(subEventUnitEpoch,2)));
    resortInd = [3 1 2];
    vennXExpanded(subEventUnitEpoch(:,resortInd), figTitle, colNames(resortInd))
    
    % Save the figure
    saveFigure(fullfile(outDir,'VennDiagrams'), sprintf('VD %s %s', params.paradigmTag, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])

    % Binary
    vennXExpanded([any(subSelInfo,2), epochSelInd], figTitle, ['headTurn Any', epochSelIndPlotVars(epoch_sel_i)])
    saveFigure(fullfile(outDir,'VennDiagrams'), sprintf('VD HT Any %s %s', params.paradigmTag, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])
    
  end
  
end


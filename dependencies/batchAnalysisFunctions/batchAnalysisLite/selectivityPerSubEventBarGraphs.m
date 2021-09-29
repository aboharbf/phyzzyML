function selectivityPerSubEventBarGraphs(selTable, params)
% a function which plots selectivity for fixed events in the params.paradigmTag, as
% defined in the events variable below. Also if desired, include
% directional saccade selectivity (probably not what we want).

outDir = fullfile(params.outputDir, 'selectivityPerSubEvent');
eventColors = [[0.9 0 0]; [0.6 0 0]; [0 0 0.8]; [0 0.6 0]; [0.6 0.6 0.6]];

% Generate bar plots showing total counts in each category
tableVars = selTable.Properties.VariableNames';
selIndLog = contains(tableVars, 'subSel_') & contains(tableVars, 'selInd') & ~contains(tableVars, 'reward') & ~contains(tableVars, 'blink');
subEventSelVars = tableVars(selIndLog);

subEvInd = contains(subEventSelVars, 'Turn') | contains(subEventSelVars, 'Contact');
subEvents = {'SubEvents' , 'eyeEvents'};

% Generate counts of selectivity in other regards (broadCat, socVNonSoc)
% for creating Venn Diagrams
epochSelInd = contains(tableVars, 'epochSel') & contains(tableVars, '_selInd') & ~contains(tableVars, 'vBase') & contains(tableVars, 'any');
epochSelIndVars = tableVars(epochSelInd);
eventArrays = {subEventSelVars(subEvInd), subEventSelVars(~subEvInd)};

for subSelType_i = 1:length(eventArrays)
  
  % Grab the events
  subSelInfo = selTable{:,eventArrays{subSelType_i}};
  
  % Names for Plotting
  subEventPlotNames = strrep(extractBetween(eventArrays{subSelType_i}, 'subSel_', '_selInd'), '_', ' ');
  
  % Count per event
  countPerSubEvent = sum(subSelInfo);
  unitCount = size(subSelInfo, 1);
  atLeastOne = sum(any(subSelInfo,2));
  barParams.labels = subEventPlotNames;
  barParams.colors = eventColors;
  
  % Plot bar graphs for individual events
  figTitle = sprintf('%s selective for Stim sub-events %s (%d/%d Unique)', params.unitTag, subEvents{subSelType_i},  atLeastOne, unitCount);
  figH = createBarPlotWithChanceLine({'Counts'}, countPerSubEvent', 0, unitCount, figTitle, subEventPlotNames, barParams);
  xlabel('Stimulus Sub-Event');
  
  % Move the legend
  legendHand = findobj(figH.Children, 'Type', 'Legend');
  legendHand.Location = 'northeastoutside';
  
  % Save the figure
  saveFigure(outDir, sprintf('BG %s %s %s', params.paradigmTag, strrep(figTitle, '/', ' of '), subEvents{subSelType_i}), [], params.figStruct, [])
  
  % Remove individual selectivity columns, keep the combo for Epoch
  % cross over selectivity
  %   event2KeepIndex = ~(contains(subEventSelVars, 'all') | ~contains(subEventSelVars, 'headTurn'))';
  if subSelType_i == 1
    event2KeepIndex = (contains(subEventPlotNames, 'headTurn') & ~contains(subEventPlotNames, 'all'));
  else
    event2KeepIndex = (contains(subEventPlotNames, 'saccade') & contains(subEventPlotNames, 'NonStim'));
  end
  subEventPlotNamesVenn = subEventPlotNames(event2KeepIndex)';
  subSelInfo = subSelInfo(:,event2KeepIndex);
  epochSelIndPlotVars = extractBetween(epochSelIndVars, 'epochSel_', '_');
  
  for epoch_sel_i = 1:length(epochSelIndVars)
    epochSelInd = selTable.(epochSelIndVars{epoch_sel_i});
    
    subEventUnitEpoch = [subSelInfo, epochSelInd];
    
    colNames = [subEventPlotNamesVenn epochSelIndPlotVars(epoch_sel_i)];
    figTitle = sprintf('%s Activity with Head Turning and %s selectivity %s (%d/%d unique)', params.unitTag,  subEvents{subSelType_i}, epochSelIndPlotVars{epoch_sel_i}, sum(any(subEventUnitEpoch,2)), length(any(subEventUnitEpoch,2)));
    resortInd = [3 1 2];
    vennXExpanded(subEventUnitEpoch(:,resortInd), figTitle, colNames(resortInd))
    
    % Save the figure
    saveFigure(outDir, sprintf('VD %s %s', params.paradigmTag, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])
    
  end
  
end

end



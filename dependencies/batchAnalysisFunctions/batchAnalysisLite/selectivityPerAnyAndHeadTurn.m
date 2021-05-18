function selectivityPerAnyAndHeadTurn(selTable, paradigm, params)
% a function which plots selectivity for fixed events in the paradigm, as
% defined in the events variable below. Also if desired, include
% directional saccade selectivity (probably not what we want).

UnitTypes = params.UnitTypes;
% UnitTypePlot = params.UnitTypePlot;
unitTypePlot = {'MUA', 'U&US'};
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

for unitType_i = 1:length(UnitTypes)
  
  if unitType_i == 1
    unitInd = contains(selTable.unitType, 'MUA');
    unitTag = 'MUA';
  else
    unitInd = ~contains(selTable.unitType, 'MUA');
    unitTag = 'U&US';
  end

  % Count per event
  subEventUnit = subSelInfo(unitInd,:);
  countPerSubEvent = sum(subEventUnit);
  unitCount = sum(unitInd);
  atLeastOne = sum(any(subEventUnit,2));

  % Plot bar graphs for individual events
  figTitle = sprintf('%s activity selective for Stimulus sub-events (%d/%d Unique)', unitTypePlot{unitType_i}, atLeastOne, unitCount);
  figH = createBarPlotWithChanceLine(subEventPlotNames, countPerSubEvent, 0, unitCount, figTitle, []);
  
  % Move the legend
  legendHand = findobj(figH.Children, 'Type', 'Legend');
  legendHand.Location = 'northeastoutside';
  
  % Save the figure
  saveFigure(fullfile(outDir,'Bargraphs'), sprintf('BG %s %s', paradigm, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])
  
  % Remove individual selectivity columns, keep the combo for Epoch
  % cross over selectivity
  event2KeepIndex = ~(contains(subEventSelVars, 'all') | ~contains(subEventSelVars, 'headTurn'))';
  subEventPlotNamesVenn = subEventPlotNames(event2KeepIndex)';
  subEventUnit = subEventUnit(:,event2KeepIndex);
  epochSelIndPlotVars = extractBetween(epochSelIndVars, 'epochSel_', '_');
  
  
  for epoch_sel_i = 1:length(epochSelIndVars)
    epochSelInd = selTable.(epochSelIndVars{epoch_sel_i})(unitInd);
    
    subEventUnitEpoch = [subEventUnit, epochSelInd];
    
    colNames = [subEventPlotNamesVenn epochSelIndPlotVars(epoch_sel_i)];
    figTitle = sprintf('%s Activity with Head Turning and %s selectivity (%d/%d unique)', unitTag, epochSelIndPlotVars{epoch_sel_i}, sum(any(subEventUnitEpoch,2)), length(any(subEventUnitEpoch,2)));
    resortInd = [3 1 2];
    vennXExpanded(subEventUnitEpoch(:,resortInd), figTitle, colNames(resortInd))
    
    % Save the figure
    saveFigure(fullfile(outDir,'VennDiagrams'), sprintf('VD %s %s', paradigm, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])

    % Binary
    vennXExpanded([any(subEventUnit,2), epochSelInd], figTitle, ['headTurn Any', epochSelIndPlotVars(epoch_sel_i)])
    saveFigure(fullfile(outDir,'VennDiagrams'), sprintf('VD HT Any %s %s', paradigm, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])
    
  end
  
end

end
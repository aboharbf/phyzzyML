function selectivityPerSubEventBarGraphs(selTable, paradigm, params)
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
selIndLog = contains(tableVars, 'subSel_') & contains(tableVars, 'selInd') & ~contains(tableVars, 'reward') & ~contains(tableVars, 'blink');
subEventSelVars = tableVars(selIndLog);
% subSelInfo = selTable{:,subEventSelVars};

% comboEventTag = {'subSel_headTurn_right_', 'subSel_headTurn_left_' , 'subSel_headTurn_all_', 'subSel_bodyTurn_', 'subSel_eyeContact_'};
% comboEventNewName = {'subSel_headTurn_right_selInd', 'subSel_headTurn_left_selInd', 'subSel_headTurn_all_selInd', 'subSel_bodyTurn_selInd', 'subSel_eyeContact_selInd'};

% for ii = 1:length(comboEventTag)
%   % Identify labels with the tags
%   mergeInd = contains(subEventSelVars, comboEventTag{ii}) & contains(subEventSelVars, 'selInd');
%
%   % Merge those columns
%   mergeCol = any(subSelInfo(:, mergeInd),2);
%
%   % Remove the two original columns from the data
%   subSelInfo = subSelInfo(:, ~mergeInd);
%   subEventSelVars = subEventSelVars(~mergeInd);
%
%   % Add new columns
%   subSelInfo = [subSelInfo, mergeCol];
%   subEventSelVars = [subEventSelVars; comboEventNewName{ii}];
%
% end

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
    figTitle = sprintf('%s selective for Stim sub-events %s (%d/%d Unique)', unitTypePlot{unitType_i}, subEvents{subSelType_i},  atLeastOne, unitCount);
    figH = createBarPlotWithChanceLine(subEventPlotNames, countPerSubEvent, 0, unitCount, figTitle, []);
    xlabel('Stimulus Sub-Event');
    
    % Move the legend
    legendHand = findobj(figH.Children, 'Type', 'Legend');
    legendHand.Location = 'northeastoutside';
    
    % Save the figure
    saveFigure(outDir, sprintf('BG %s %s %s', paradigm, strrep(figTitle, '/', ' of '), subEvents{subSelType_i}), [], params.figStruct, [])
    
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
      figTitle = sprintf('%s Activity with Head Turning and %s selectivity %s (%d/%d unique)', unitTag,  subEvents{subSelType_i}, epochSelIndPlotVars{epoch_sel_i}, sum(any(subEventUnitEpoch,2)), length(any(subEventUnitEpoch,2)));
      resortInd = [3 1 2];
      vennXExpanded(subEventUnitEpoch(:,resortInd), figTitle, colNames(resortInd))
      
      % Save the figure
      saveFigure(outDir, sprintf('VD %s %s', paradigm, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])
      
    end
  end
end

end

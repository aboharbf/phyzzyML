function   jointTuningPlot(selTable, paradigm, params)
% A function which creates visuals demonstrating the selectivity of each
% unit, MUA, and unsorted trace, alongside venn diagrams making the right
% comparisons.

outDir = fullfile(params.outputDir, 'jointTuningPlot');

columnsInTable = {'subSel_fixation_selInd', 'subSel_rewardCombo_selInd', 'subSel_saccades_selInd', ...
    'subSel_headTurn_all_selInd', 'saccDir_all_selInd'};
tableLabels = {'Fix Dot', 'Reward', 'Saccades', ...
    'headTurning', 'Saccade, Directional'};

% Combine across the appropriate ones (mostly the saccade ones)
comboparams.comboEvents = {'subSel_headTurn_all_selInd', };
comboparams.comboSubEvents = {{'subSel_headTurn_all_saccNone_selInd', 'subSel_headTurn_all_sacc_selInd'}};
selTable = expandSelTableComboEvents(selTable, comboparams);

% The image to be produced will be a stack of images, where every pixel is
% either colored in or black. The colors will vary by monkey and unit type.
% Resort the selTable accordingly.
selTableVars = selTable(:, columnsInTable);
dateSubjData = selTable{:, 'dateSubj'};
unitData = selTable{:, 'unitType'};

% Set which combinations will be displayed in Venn Diagrams
vennCombos = {[1 2 3], [3 4 5], [1 4 2]}; 

% Mo, US, Units, MUA
MoUSInd = find(contains(dateSubjData, 'Mo') & contains(unitData, 'US'));
MoUInd = find(contains(dateSubjData, 'Mo') & contains(unitData, digitsPattern));
MoMUAInd = find(contains(dateSubjData, 'Mo') & contains(unitData, 'MUA'));

SamUSInd = find(contains(dateSubjData, 'Sam') & contains(unitData, 'US'));
SamUInd = find(contains(dateSubjData, 'Sam') & contains(unitData, digitsPattern));
SamMUAInd = find(contains(dateSubjData, 'Sam') & contains(unitData, 'MUA'));

resortInd = [MoUSInd; MoUInd; MoMUAInd; SamUSInd; SamUInd; SamMUAInd];

monkey = {'Sam', 'Mo'};
unitTypes = {'MUA', digitsPattern'};
unitTypesPlot = {'MUA', 'Units'};
setInd = 1;

for mo_i = 1:length(monkey)
  for unit_i = 1:length(unitTypes)
    
    % Identify units, pull data
    unitInd = contains(dateSubjData, monkey{mo_i}) & contains(unitData, unitTypes{unit_i});
    sigIndArray = selTableVars{unitInd, :};
    
    % Plot the Image version
    figTitle = sprintf('Joint Selectivity Array - %s - %d %s', monkey{mo_i}, size(sigIndArray,1), unitTypesPlot{unit_i});
    figH = figure('Name', figTitle, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0.0047 0.0370 0.4708 0.9444]);
    imagesc(sigIndArray);
    figH.Children.FontSize = 15;
    title(figTitle)
    
    % Labels
    varCounts = sum(sigIndArray);
    xLabels = strcat(tableLabels', ' (' , string(varCounts'), ')');
    xticks(1:length(tableLabels))
    xticklabels(xLabels);
    xtickangle(45);
    xlabel('Event')
    ylabel('Unit Number')
    
    saveFigure(outDir, sprintf('%d.JS_Array %s - %s', setInd, paradigm, figTitle), [], params.figStruct, [])
    
    % Generate the possible pairing, and create a venn diagram for each
    for venn_i = 1:length(vennCombos)
      
      eventInd = vennCombos{venn_i};
      sigInd = sigIndArray(:, eventInd);
      colNames = tableLabels(eventInd);
      figTitle = '';
      vennXExpanded(sigInd, figTitle, colNames)
      saveFigure(outDir, sprintf('%d.VD %s - %s', setInd, paradigm, strjoin(colNames, '_')), [], params.figStruct, [])
      
    end
    
    setInd = setInd + 1;
    
  end
end


end
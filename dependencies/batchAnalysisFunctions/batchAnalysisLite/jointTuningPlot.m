function   jointTuningPlot(selTable, paradigm, params)
% A function which creates visuals demonstrating the selectivity of each
% unit, MUA, and unsorted trace, alongside venn diagrams making the right
% comparisons.

outDir = fullfile(params.outputDir, 'jointTuningPlot');

columnsInTable = {'subSel_fixation_selInd', 'subSel_rewardCombo_selInd', 'subSel_saccades_selInd', ...
    'subSel_headTurn_all_selInd', 'saccDir_all_selInd', 'epochSel_categories_any_selInd'};
tableLabels = {'Fix Dot', 'Reward', 'Saccades', ...
    'headTurning', 'Saccade, Directional', 'Category Selective'};

% Combine across the appropriate ones (mostly the saccade ones)
% comboparams.comboEvents = {'subSel_headTurn_all_selInd'};
% comboparams.comboSubEvents = {{'subSel_headTurn_all_selInd', 'subSel_headTurn_all_sacc_selInd'}};
% selTable = expandSelTableComboEvents(selTable, comboparams);

% The image to be produced will be a stack of images, where every pixel is
% either colored in or black. The colors will vary by monkey and unit type.
% Resort the selTable accordingly.
selTableVars = selTable{:, columnsInTable};

% Set which combinations will be displayed in Venn Diagrams
vennCombos = {[1 2 3], [3 4 5], [1 4 2]};

% Plot the Image version
figTitle = sprintf('Joint Selectivity Array - %s - %d %s', params.monkeyTag, size(selTable, 1), params.unitTag);
figH = figure('Name', figTitle, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0.0047 0.0370 0.4708 0.9444]);
imagesc(selTableVars);
figH.Children.FontSize = 15;
title(figTitle)

% Labels
varCounts = sum(selTableVars);
xLabels = strcat(tableLabels', ' (' , string(varCounts'), ')');
xticks(1:length(tableLabels))
xticklabels(xLabels); xtickangle(45); xlabel('Event'); ylabel('Unit Number')

saveFigure(outDir, sprintf('JS_Array %s - %s', paradigm, figTitle), [], params.figStruct, [])

% Generate the possible pairing, and create a venn diagram for each
for venn_i = 1:length(vennCombos)
  
  eventInd = vennCombos{venn_i};
  sigInd = selTableVars(:, eventInd);
  colNames = tableLabels(eventInd);
  figTitle = sprintf('VD %s - %s', paradigm, strjoin(colNames, '_'));
  vennXExpanded(sigInd, figTitle, colNames)
  saveFigure(outDir, figTitle, [], params.figStruct, [])

end

end
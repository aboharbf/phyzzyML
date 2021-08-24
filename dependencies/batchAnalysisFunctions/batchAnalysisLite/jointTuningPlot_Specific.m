function   jointTuningPlot_Specific(selTable, paradigm, params)
% A function which creates visuals demonstrating the selectivity of each
% unit, MUA, and unsorted trace, alongside venn diagrams making the right
% comparisons.

outDir = fullfile(params.outputDir, 'jointTuningPlot');

columnsInTable = {{'subSel_rewardCombo_selInd', 'epochSel_categories_any_selInd', {'subSel_saccadesNonStim_selInd' 'subSel_pre-saccadesNonStim_selInd'}} ...
  {'subSel_rewardCombo_selInd', 'epochSel_categories_any_selInd', {'subSel_saccades_selInd' 'subSel_pre-saccades_selInd'}}};
tableLabels = {{'Reward', 'Category Selective', 'Pre + Post Saccade (NonStim)'} {'Reward', 'Category Selective', 'Pre + Post Saccade'}};
figTitleString = {'Joint Selectivity Array (NonStim)', 'Joint Selectivity Array'};

% Combine across the appropriate ones (mostly the saccade ones)
% comboparams.comboEvents = {'subSel_headTurn_all_selInd'};
% comboparams.comboSubEvents = {{'subSel_headTurn_all_selInd', 'subSel_headTurn_all_sacc_selInd'}};
% selTable = expandSelTableComboEvents(selTable, comboparams);

% The image to be produced will be a stack of images, where every pixel is
% either colored in or black. The colors will vary by monkey and unit type.
% Resort the selTable accordingly.

% Set which combinations will be displayed in Venn Diagrams

for ii = 1:length(figTitleString)
  % Columns to contribute
  col2Merge = columnsInTable{ii};
  colData = nan(size(selTable,1),0);
  for col_i = 1:length(col2Merge)
    % Collect info
    colData = [colData, any(selTable{:, col2Merge{col_i}}, 2)];
  end
  
  
  % Plot the Image version
  figTitle = sprintf('%s - %s - %d %s', figTitleString{ii}, params.monkeyTag, size(selTable, 1), params.unitTag);
  figH = figure('Name', figTitle, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0.0047 0.0370 0.4708 0.9444]);
  imagesc(colData);
  figH.Children.FontSize = 15;
  title(figTitle)
  
  % Labels
  varCounts = sum(colData);
  xLabels = strcat(tableLabels{ii}', ' (' , string(varCounts'), ')');
  xticks(1:length(tableLabels{ii}))
  xticklabels(xLabels); xtickangle(45); xlabel('Event'); ylabel('Unit Number')
  
  % Save
  saveFigure(outDir, sprintf('JS_Array %s - %s', paradigm, figTitle), [], params.figStruct, [])
  
  % Generate the possible pairing, and create a venn diagram for each
  colNames = tableLabels{ii};
  figTitle = sprintf('%s %s - %s', figTitleString{ii}, paradigm, strjoin(colNames, ' '));
  vennXExpanded(colData, figTitle, colNames)
  
  % Save
  saveFigure(outDir, figTitle, [], params.figStruct, [])
  
end

end
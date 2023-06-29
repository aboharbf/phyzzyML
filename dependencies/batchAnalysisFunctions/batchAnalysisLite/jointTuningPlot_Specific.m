function   jointTuningPlot_Specific(selTable, paradigm, params)
% A function which creates visuals demonstrating the selectivity of each
% unit, MUA, and unsorted trace, alongside venn diagrams making the right
% comparisons.

outDir = fullfile(params.outputDir, 'jointTuningPlot');

selIdxPresent = selTable.Properties.VariableNames(contains(selTable.Properties.VariableNames, 'selInd'));

columnsInTable = {...
  %{'subSel_rewardCombo_selInd', 'epochSel_categories_any_selInd', {'subSel_saccadesNonStim_selInd' 'subSel_pre_saccadesNonStim_selInd'}} ...
  %{'subSel_rewardCombo_selInd', 'epochSel_categories_any_selInd',                 {'subSel_saccades_selInd' 'subSel_pre_saccades_selInd'}}...
  %{'subSel_rewardCombo_selInd', 'subSel_headTurn_all_selInd',                     {'subSel_saccadesNonStim_selInd' 'subSel_pre_saccadesNonStim_selInd'}} ...
  {'subSel_rewardCombo_selInd', 'subSel_headTurn_all_selInd',                     {'subSel_saccades_selInd' 'subSel_pre_saccades_selInd'}}...
  %{'subSel_rewardCombo_selInd', 'subSel_stimSubEvent_all_selInd',                 {'subSel_saccadesNonStim_selInd' 'subSel_pre_saccadesNonStim_selInd'}} ...
  {'subSel_rewardCombo_selInd', 'subSel_stimSubEvent_all_selInd',                 {'subSel_saccades_selInd' 'subSel_pre_saccades_selInd'}}...
  %{'epochSel_categories_any_selInd' 'subSel_headTurn_all_selInd'}...
  %{'epochSel_categories_any_selInd' 'subSel_headTurn_all_selInd',                 {'subSel_saccadesNonStim_selInd' 'subSel_pre_saccadesNonStim_selInd'}}...
  %{'epochSel_categories_any_selInd' 'subSel_stimSubEvent_all_selInd'}...
  %{'epochSel_categories_any_selInd' 'subSel_pre_saccadesNonStim_selInd'}...
  %{'epochSel_categories_any_selInd' 'subSel_headTurn_right_selInd',                'subSel_headTurn_left_selInd'}...
};
    
vars2Plot = flattenCellArray(columnsInTable);
missingVar = setdiff(vars2Plot, selIdxPresent);

assert(isempty(missingVar), sprintf('Missing some variables from the selTable to make all the currently defined jointPlots \n%s', strjoin(missingVar, '   ')))

tableLabels = {...
                %{'Reward', 'Category Selective', 'Pre/Post Saccade (NonStim)'}...
                %{'Reward', 'Category Selective', 'Pre/Post Saccade'},...
                %{'Reward', 'Head Turn Selective', 'Pre/Post Saccade (NonStim)'}...
                {'Reward', 'Head Turn Selective', 'Pre/Post Saccade'}...
                %{'Reward', 'Stimulus Event', 'Pre/Post Saccade (NonStim)'}...
                {'Reward', 'Stimulus Event', 'Pre/Post Saccade'}...
                %{'Category Selective', 'Head Turn Selective'}...
                %{'Category', 'Head Turn', 'Pre/Post Saccade'}...
                %{'Category Selective', 'Stimulus Event Selective'}...
                %{'Category Selective', 'Pre-Saccade Selectivity'}...
                %{'Category Selective', 'headTurn right', 'headTurn left'}...
                };
              
figTitleString = {...
                  %'Joint Selectivity Array A (NonStim)', 'Joint Selectivity Array A', ...
                  %'Joint Selectivity Array B (NonStim)',...
                  'Joint Selectivity B Array', ...
                  %'Joint Selectivity Array C (NonStim)',...
                  'Joint Selectivity C Array', ...
                  %'Joint Selectivity D Array (HT)',    'Joint Selectivity F Array (HT)', ...
                  %'Joint Selectivity D Array (SubEv)', 'Joint Selectivity D Array (Event + PreSacc)',...
                  %'Joint Selectivity E Array (Cat + HT)'
                  };

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
  if 0
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
  end

  % Generate the possible pairing, and create a venn diagram for each
  colNames = tableLabels{ii};
  
  if 0
    createBarPlotWithChanceLine(colNames, sum(colData), 0.1, size(colData,1), 'bar', []);
    figH = gcf;
    figH.Children(2).XTickLabelRotation = 0;
    
  end
  
  figTitle = sprintf('%s %s - %s', figTitleString{ii}, paradigm, strjoin(colNames, ' '));
  colNames = strrep(colNames, ' (NonStim)', '');
  vennXExpanded(colData, figTitle, colNames)
  
  % Save
  saveFigure(outDir, strrep(figTitle, '/', ' and '), [], params.figStruct, [])

end

end
function selectivityPerEpochBarGraphsAndVennDiagram(selTable, paradigm, params)
% a function which grabs counts of determined selectivities from the
% selTable passed in, and generates bar plots for the different entries on
% a per epoch basis (stim early, late, and reward). columns are generated
% during runAnalysis > epochStats, and it is an epoch based ANOVA (as
% opposed to the sliding scale, which is more granular).

% Inputs:
% - params, containing data which is extracted below.
% - paradigm, a string denoting the paradigm (must match a struct in
% params.

selCheck = params.(params.paradigmTag).selCheck;
colNamePoss = params.colNamePoss;
colNamePlotAll = params.colNamePlotAll;
alpha = 0.05;                             % This is the alpha from the previous processRunBatch. Change it accordingly.
tableVars = selTable.Properties.VariableNames';
outDir = fullfile(params.outputDir, 'selectivityPerEpoch');

uniqueLabelsSorted = {{'socialInteraction', 'agents'}, {'chasing','fighting', 'mounting', 'grooming', 'goalDirected', 'idle', 'objects', 'scene'}, {'socialInteraction', 'goalDirected', 'idle', 'objects', 'scene'}};

% Generate bar plots showing total counts in each category

unitCount = size(selTable,1);

% Sanity Check on atLeastOne being the same across 2 checks
atleastOneVec = nan(unitCount,length(selCheck));

for sel_i = 1:length(selCheck)
  % Identify which columns were part of this analysis
  comparisonIndex = contains(tableVars, ['epochSel_' selCheck{sel_i}]) & ~contains(tableVars, '_any_') & ~contains(tableVars, '_vBase_');
  selColumns = tableVars(comparisonIndex & contains(tableVars, 'selInd'))';
  diffColumns = tableVars(comparisonIndex & contains(tableVars, 'diff'))';
  stimPrefColumns = tableVars(comparisonIndex & contains(tableVars, '_prefStim'))';
  
  % Assuming each only has 3 possible comparisons (onset, pres, reward).
  % Create the groups + their intersections
  colNames = extractAfter(selColumns, sprintf('epochSel_%s_', selCheck{sel_i}));
  colNames = extractBefore(colNames, '_selInd');
  
  [~, A] = intersect(colNames, colNamePoss);
  A = sort(A);
  colNamesPlot = colNames(A);
  selColumns = selColumns(A);
  sigInd = selTable{:, selColumns};
  
  if ~isempty(diffColumns)
    diffColumns = diffColumns(A);
    possMat = selTable{:, diffColumns} > 0;
  end
  
  if ~isempty(stimPrefColumns)
    
    stimPrefColumnData = selTable{:, stimPrefColumns};
    
    % Find unique labels and create a 3rd dimension for each.
    %       uniqueLabels = unique(stimPrefColumnData(:));
    uniqueLabels = uniqueLabelsSorted{sel_i};
    sigIndTmp = zeros([size(stimPrefColumnData), length(uniqueLabels)]);
    for lab_i = 1:length(uniqueLabels)
      sigIndTmp(:, :, lab_i) = strcmp(stimPrefColumnData, uniqueLabels(lab_i));
    end
    
    % Make the dataMat
    objPreMat = squeeze(sum(sigIndTmp,1));
    objPreflegendLabels = uniqueLabels;
    %       atleastOneVec(:,sel_i) = logical(sum(stimPrefColumnData, 2:3));
    
  else
    if ~isempty(diffColumns)
      dataMat = [sum(sigInd & possMat)', sum(sigInd & ~possMat)'];
      legendLabels = {'Increased firing', 'decreased firing'};
      alphaPlot = alpha;
    else
      dataMat = sum(sigInd)';
      legendLabels = [];
      alphaPlot = alpha * 2;
    end
  end
  
  atleastOneVec(:,sel_i) = any(sigInd, 2);
  atleastOne = sum(atleastOneVec(:,sel_i));
  
  % Plot the Bar graphs
  if isempty(stimPrefColumns)
    figTitle = sprintf('%s selectivity during %s test (%d/%d Unique)', params.unitTag, selCheck{sel_i}, atleastOne, size(atleastOneVec,1));
    createBarPlotWithChanceLine(colNamesPlot, dataMat, alphaPlot, unitCount, figTitle, legendLabels);
    saveFigure(fullfile(outDir, 'BarGraph'), sprintf('BG %s %s', paradigm, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])
  end
  
  % If there is object preference data, plot here
  if ~isempty(stimPrefColumns)
    figTitle = sprintf('%s selectivity for Objects during %s test (%d/%d Unique)', params.unitTag, selCheck{sel_i}, atleastOne, size(atleastOneVec,1));
    createBarPlotWithChanceLine(colNamesPlot, objPreMat, 0, unitCount, figTitle, objPreflegendLabels);
    saveFigure(fullfile(outDir, 'objPref'), sprintf('1b.BG %s %s', paradigm, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])
  end
  
  % Plot the Venn diagram
  figTitle = sprintf('%s selectivity for %s (%d/%d Unique)', params.unitTag, selCheck{sel_i}, atleastOne, size(atleastOneVec,1));
  if size(sigInd, 2) == 4
    sigInd = sigInd(:,2:end);
  end
  vennXExpanded(sigInd, figTitle, colNames)
  saveFigure(fullfile(outDir, 'VennDiagram'), sprintf('%s %s', paradigm, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])
  
end


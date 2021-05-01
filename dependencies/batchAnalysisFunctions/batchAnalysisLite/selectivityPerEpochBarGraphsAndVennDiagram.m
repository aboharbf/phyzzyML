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

selCheck = params.(paradigm).selCheck;
UnitTypes = params.UnitTypes;
UnitTypePlot = params.UnitTypePlot;
colNamePoss = params.colNamePoss;
colNamePlotAll = params.colNamePlotAll;
alpha = 0.05;                             % This is the alpha from the previous processRunBatch. Change it accordingly.
tableVars = selTable.Properties.VariableNames';
outDir = fullfile(params.outputDir, 'selectivityPerEpoch');

if ~exist(outDir, 'dir')
  mkdir(outDir);
end

% Generate bar plots showing total counts in each category
for unitType_i = 1:length(UnitTypes)
  % Filter the table per unit
  if unitType_i == 1
    unitInd = contains(selTable.unitType, 'MUA');
  else
    unitInd = ~contains(selTable.unitType, 'MUA');
  end
  
  selTableParadigmUnit = selTable(unitInd, :);
  
  unitCount = sum(unitInd);
  
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
    colNamesPlot = colNamePlotAll(A);
    selColumns = selColumns(A);
    sigInd = selTableParadigmUnit{:, selColumns};
    
    if ~isempty(diffColumns)
      diffColumns = diffColumns(A);
      possMat = selTableParadigmUnit{:, diffColumns} > 0;
    end
    
    if ~isempty(stimPrefColumns)
      
      stimPrefColumnData = selTableParadigmUnit{:, stimPrefColumns};
      
      % Find unique labels and create a 3rd dimension for each.
      uniqueLabels = unique(stimPrefColumnData(:));
      sigIndTmp = zeros([size(stimPrefColumnData), length(uniqueLabels)]);
      for lab_i = 1:length(uniqueLabels)
        sigIndTmp(:, :, lab_i) = strcmp(stimPrefColumnData, uniqueLabels(lab_i));
      end
      
      % Get rid of 'None'.
      stimPrefColumnData = sigIndTmp(:, :, ~strcmp(uniqueLabels, 'None'));
      uniqueLabels = uniqueLabels(~strcmp(uniqueLabels, 'None'));
      
      % Make the dataMat
      objPreMat = squeeze(sum(stimPrefColumnData,1));
      objPreflegendLabels = uniqueLabels;
%       atleastOneVec(:,sel_i) = logical(sum(stimPrefColumnData, 2:3));
%       atleastOne = sum(atleastOneVec(:,sel_i));
      
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
      atleastOneVec(:,sel_i) = any(sigInd, 2);
      atleastOne = sum(atleastOneVec(:,sel_i));
      
    end
    
    % Plot the Bar graphs
    figTitle = sprintf('Bar Graph - %s activity selective for %s during %s (%d Unique)', UnitTypePlot{unitType_i}, selCheck{sel_i}, paradigm, atleastOne);
    createBarPlotWithChanceLine(colNamesPlot, dataMat, alphaPlot, unitCount, figTitle, legendLabels)
    saveFigure(outDir, figTitle, [], params.figStruct, [])
    
    % If there is object preference data, plot here
    if ~isempty(stimPrefColumns)
      figTitle = sprintf('Bar Graph - %s activity selective for Objects in %s during %s (%d Unique)', UnitTypePlot{unitType_i}, selCheck{sel_i}, paradigm, atleastOne);
      createBarPlotWithChanceLine(colNamesPlot, objPreMat, 0, unitCount, figTitle, objPreflegendLabels)
      saveFigure(outDir, figTitle, [], params.figStruct, [])
    end
    
    % Plot the Venn diagram    
    atleastOne = sum(any(sigInd, 2));
    figTitle = sprintf('Venn Diagram - %s activity selective for %s during %s (%d Unique)', UnitTypePlot{unitType_i}, selCheck{sel_i}, paradigm, atleastOne);
    vennXExpanded(sigInd, figTitle, colNames)
    saveFigure(outDir, figTitle, [], params.figStruct, [])
    
  end
  
end
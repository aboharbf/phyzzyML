function selectivityPerEpochBarGraphs(selTable, paradigm, params)
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
alpha = 0.01;                             % This is the alpha from the previous processRunBatch. Change it accordingly.

% Generate bar plots showing total counts in each category
for unitType_i = 1:length(UnitTypes)
  % Filter the table per unit
  unitInd = contains(selTable.unitType, UnitTypes{unitType_i});
  selTableParadigmUnit = selTable(unitInd, :);
  
  if ~isempty(selTableParadigmUnit)
    % Plots 1 - fixation/saccade
    unitCount = sum(unitInd);
    chanceUnitCount = round(sum(unitCount)*(alpha/2));
    
    for sel_i = 1:length(selCheck)
      % Identify which columns were part of this analysis
      columns2Combine = contains(selTableParadigmUnit.Properties.VariableNames, [selCheck{sel_i} 'Sel_']);
      
      if any(columns2Combine)
        % Assuming each only has 3 possible comparisons (onset, pres,
        % reward at most).
        
        % Create the groups + their intersections
        selTableParadigmVarNames = selTableParadigmUnit.Properties.VariableNames(columns2Combine);
        colNames = extractAfter(selTableParadigmVarNames, [selCheck{sel_i} 'Sel_']);
        %         colNames = colNames(~(strcmp(colNames, 'stimWhole') | strcmp(colNames, 'any'))); % Don't include these.
        [~, A] = intersect(colNames, colNamePoss);
        A = sort(A);
        colNamesPlot = colNamePlotAll(A);
        selTableParadigmVarNames = selTableParadigmVarNames(A);
        sigInd = selTableParadigmUnit{:, selTableParadigmVarNames};
        
        if iscell(sigInd)
          
          % Find unique labels and create a 3rd dimension for each.
          uniqueLabels = unique(sigInd(:));
          sigIndTmp = zeros([size(sigInd), length(uniqueLabels)]);
          for lab_i = 1:length(uniqueLabels)
            sigIndTmp(:, :, lab_i) = strcmp(sigInd, uniqueLabels(lab_i));
          end
          
          % Get rid of 'None'.
          sigInd = sigIndTmp(:, :, ~strcmp(uniqueLabels, 'None'));
          uniqueLabels = uniqueLabels(~strcmp(uniqueLabels, 'None'));
          
          % Make the dataMat
          dataMat = squeeze(sum(sigInd,1));
          legendLabels = uniqueLabels;
          
        else
          % Get rid of the NaNs.
          sigInd(isnan(sigInd)) = 0;
%           sigInd = logical(sigInd);
          
          dataMat = [sum(sigInd > 0)',sum(sigInd < 0)'];
          atleastOne = sum(logical(sum(sigInd ~= 0, 2)));
          
          legendLabels = {'Increased firing', 'decreased firing'};
          
        end
                
        % Plot
        figTitle = sprintf('%s activity selective for %s during %s (%d Unique)', UnitTypePlot{unitType_i}, selCheck{sel_i}, paradigm, atleastOne);
        createBarPlotWithChanceLine(colNamesPlot, dataMat, alpha, unitCount, figTitle, legendLabels)
        saveFigure(params.outputDir, figTitle, [], params.figStruct, [])
        
      end
    end
  end
  
end
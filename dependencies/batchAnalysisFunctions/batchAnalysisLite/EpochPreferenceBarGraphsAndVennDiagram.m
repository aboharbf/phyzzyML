function EpochPreferenceBarGraphsAndVennDiagram(selTable, paradigm, params)
% a function which grabs counts of determined selectivities from the
% selTable passed in, and generates bar plots describing the preferred
% epoch activity.
% during runAnalysis > epochStats, and it is an epoch based ANOVA (as
% opposed to the sliding scale, which is more granular).

% Inputs:
% - params, containing data which is extracted below.
% - paradigm, a string denoting the paradigm (must match a struct in
% params.

epochList = {'Fix', 'stimOnset', 'stimPres', 'reward'};

varNames = selTable.Properties.VariableNames';
epochPrefVars = varNames(contains(varNames, 'epochPref_'));
BaseVsVars = varNames(contains(varNames, 'baseV_'));
unitTypePlot = {'MUA', 'U&US'};

% Generate bar plots showing total counts in each category
for unitType_i = 1:2
  
  if unitType_i == 1
    unitInd = contains(selTable.unitType, 'MUA');
  else
    unitInd = ~contains(selTable.unitType, 'MUA');
  end
  
  % Filter the table per unit
  selTableParadigmUnit = selTable(unitInd, :);
  identifierVector = strcat(selTableParadigmUnit.dateSubj,selTableParadigmUnit.runNum, selTableParadigmUnit.channel);
  unitCount = sum(unitInd);
    
  % Extract pVals and means for each comparison  
  pValMat = selTableParadigmUnit{:, BaseVsVars(contains(BaseVsVars, '_pVal'))};
  sigMat = pValMat <= 0.05;
  sigBaselineCount = sum(logical(sum(sigMat, 2)));
  
%   sigMat = selTableParadigmUnit{:, BaseVsVars(contains(BaseVsVars, '_selInd'))};
%   sigBaselineCount = sum(any(sigMat,2));
  sigPerEpochCount = sum(sigMat);
  
  % Plot 1
  alpha = 0.05;
  figTitle = sprintf('%s - %s modulated during each Epoch, alpha = %s (%d of %d Unique)', paradigm, unitTypePlot{unitType_i}, num2str(alpha), sigBaselineCount, unitCount);
  createBarPlotWithChanceLine(epochList, sigPerEpochCount, alpha, unitCount, figTitle, [])
  saveFigure(params.outputDir, figTitle, [], params.figStruct, [])

  % Venn Diagram
  sigMatVenn = [sigMat(:,1), [sigMat(:,2) | sigMat(:,3)], sigMat(:,4)];
  colNames = {'Fix', 'stimAny', 'Reward'};
  figTitle = sprintf('%s - %s modulated during each Epoch, alpha = %s - Intersection ', paradigm, unitTypePlot{unitType_i}, num2str(alpha));
  vennXExpanded(sigMatVenn, figTitle, colNames)
  saveFigure(params.outputDir, figTitle, [], params.figStruct, [])
  
  % Plot 2 - preference index
  prefIndex = selTableParadigmUnit.(epochPrefVars{2}); % highest epoch - mean of the rest / sum;
  epochListSig = selTableParadigmUnit.(epochPrefVars{1});
  
  countPerEpoch = nan(size(epochList));
  figTitle = sprintf('%s - %s Preference Index per Epoch (0.3 = 2x response)', paradigm, unitTypePlot{unitType_i});
  figH = figure('name', figTitle, 'units', 'normalized', 'position', [0.35 0.5 0.6 0.4]);
  sgtitle(figTitle);
  for ep_i = 1:length(epochList)
    
    % find the count for this particular epoch
    indTmp = strcmp(epochListSig, epochList{ep_i});
    indSig = indTmp & sigMat(:,ep_i);
    indSigAllInd = find(indSig);
    
    % Plotting
    subplot(1, length(epochList), ep_i)
    a = histogram(prefIndex(indTmp));
    legendEntryA = sprintf('All (%d)', sum(indTmp));
    hold on
    b = histogram(prefIndex(indSig));
    b.BinEdges = a.BinEdges;
    legendEntryB = sprintf('Sig > baseline (%d)', sum(indSig));

    % Add Labels to plot
    legend(legendEntryA, legendEntryB)
    countPerEpoch(ep_i) = sum(indTmp);
    title(sprintf('Indicies for %s (%d total)', epochList{ep_i}, sum(indTmp)))
    
    % strongest preference report
    [A, B] = sort(prefIndex(indSig), 'descend');
    for b_i = 1:5
      fprintf('Max Pref %s = %s @ %s \n', epochList{ep_i}, num2str(A(b_i)),  identifierVector{indSigAllInd(B(b_i))})
    end
  end
  saveFigure(params.outputDir, figTitle, [], params.figStruct, [])

end

end
function epochPreferenceBarGraphsAndVennDiagram(selTable, params)
% a function which grabs counts of determined selectivities from the
% selTable passed in, and generates bar plots describing the preferred
% epoch activity.
% during runAnalysis > epochStats, and it is an epoch based ANOVA (as
% opposed to the sliding scale, which is more granular).

% Inputs:
% - params, containing data which is extracted below.
% - params.paradigmTag, a string denoting the params.paradigmTag (must match a struct in
% params.

outDir = fullfile(params.outputDir, 'epochPreference');

epochList = {'Fix', 'stimEarly', 'stimLate', 'reward'};
varNames = selTable.Properties.VariableNames';
epochPrefVars = varNames(contains(varNames, 'epochPref_'));
BaseVsVars = varNames(contains(varNames, 'baseV_'));

% Generate bar plots showing total counts in each category

% Filter the table per unit
identifierVector = strcat(selTable.dateSubj, selTable.runNum, selTable.channel);
unitCount = size(selTable,1);

% Extract pVals and means for each comparison
sigMat = selTable{:, BaseVsVars(contains(BaseVsVars, '_selInd'))};
diffMat = selTable{:, BaseVsVars(contains(BaseVsVars, '_diff'))};
possMat = diffMat > 0;

sigBaselineCount = sum(any(sigMat,2));
sigPerEpochCountPos = sum(sigMat & possMat);
sigPerEpochCountNeg = sum(sigMat & ~possMat);

dataMat = [sigPerEpochCountPos; sigPerEpochCountNeg];
legendLabels = {'Increased firing', 'decreased firing'};

% Plot 1
alpha = 0.05;
figTitle = sprintf('%s modulated during each Epoch (%d/%d Unique, a = %s)', params.unitTag, sigBaselineCount, unitCount, num2str(alpha));
createBarPlotWithChanceLine(epochList, dataMat, alpha, unitCount, figTitle, legendLabels);
modDir = fullfile(outDir, 'Modulation');
saveFigure(modDir, sprintf('BG %s - %s', params.paradigmTag, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])

% Venn Diagram
sigMatVenn = [sigMat(:,1), sigMat(:,2) | sigMat(:,3), sigMat(:,4)];
colNames = {'Fix', 'stimAny', 'Reward'};
figTitle = sprintf('%s modulated during each Epoch (%d/%d Unique, a = %s)', params.unitTag, sigBaselineCount, unitCount, num2str(alpha));
vennXExpanded(sigMatVenn, figTitle, colNames)
saveFigure(modDir, sprintf('VD %s - %s', params.paradigmTag, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])

%% Plot 2 - preference index
prefIndex = selTable.(epochPrefVars{2}); % highest epoch - mean of the rest / sum;
epochListSig = selTable.(epochPrefVars{1});

countPerEpoch = nan(size(epochList));
figTitle = sprintf('%s Preference Index per Epoch (0.3 = 2x response)', params.unitTag);
figure('name', figTitle, 'NumberTitle', 'off', 'units', 'normalized', 'position', [0.35 0.5 0.6 0.4]);
sgtitle(figTitle);

% Find the vBase line selInd
tableVars = selTable.Properties.VariableNames';
tableVarsBase = tableVars(contains(tableVars, '_vBase_selInd'));

for ep_i = 1:length(epochList)
  
  % find the count for this particular epoch
  if 1%ep_i == 1
    indTmp = strcmp(epochListSig, epochList{ep_i});
  else
    indTmp = selTable.(tableVarsBase{ep_i-1});
  end
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

saveFigure(fullfile(outDir, 'TargetBaseline'), sprintf('BG %s - %s', params.paradigmTag, strrep(figTitle, '/', ' of ')), [], params.figStruct, [])

end
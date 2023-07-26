function [analysisOutFilename] = runBatchAnalysisLite(spikePathBank, batchAnalysisParams)
%% Overwrite switches with what is currently in file, keeping the spikePathLoadParams, which holds most important stuff.
spikePathLoadParamsKeep = batchAnalysisParams.spikePathLoadParams;
batchAnalysisParams = load(fullfile(batchAnalysisParams.outputDir, 'batchAnalysisParams.mat'));
batchAnalysisParams.spikePathLoadParams = spikePathLoadParamsKeep;

% Load basic inputs
plotSwitch = batchAnalysisParams.plotSwitch;
calcSwitch = batchAnalysisParams.calcSwitch;
figStruct = batchAnalysisParams.figStruct;

%% Pre-Analysis processing

if ~contains('stimPresCount',spikePathBank.Properties.VariableNames)
  spikePathBank = stimulusStatistics(spikePathBank, batchAnalysisParams.stimStructParams);
  spikePathFile = batchAnalysisParams.spikePathLoadParams.batchAnalysisOutput;
  save(spikePathFile, 'spikePathBank', '-append')
end

% stimulusEventCounts = stimulusPresence(spikePathBank);

if ~contains('selTable', spikePathBank.Properties.VariableNames)
  spikePathBank = processAppendSelTable(spikePathBank, batchAnalysisParams);
  spikePathFile = batchAnalysisParams.spikePathLoadParams.batchAnalysisOutput;
  save(spikePathFile, 'spikePathBank', '-append')
end

% spikePathBank.selTable = removeTableVariable(spikePathBank.selTable, {'spikeEyeCorr', 'spikeEyeCorrNull'});

%% Counting
% crossParadigmCheck(spikePathBank, batchAnalysisParams)

if plotSwitch.spikeEyeCorr
  spikeEyeCorrHist(spikePathBank, batchAnalysisParams)
end

if plotSwitch.spikeEyeCorr
  lookingPatternsAnalysis(spikePathBank, batchAnalysisParams)
end

if plotSwitch.rampingAnalysis
  rampingAnalysis(spikePathBank, batchAnalysisParams)
end

if plotSwitch.rewardAnalysis
  rewardAnalysis(spikePathBank, batchAnalysisParams)
end

if plotSwitch.latencyAnalysis
  latencyAnalysis(spikePathBank, batchAnalysisParams)
end

if plotSwitch.saccadeDirAnalysis
  saccadeDirAnalysis(spikePathBank, batchAnalysisParams)
  FixStimOnsetAnalysis(spikePathBank, batchAnalysisParams)
  StimOnSetSaccadeAnalysis(spikePathBank, batchAnalysisParams)
  stimOffsetRewardAnalysis(spikePathBank, batchAnalysisParams)
  stimSubEventAnalysis(spikePathBank, batchAnalysisParams)
  socSaccAnalysis(spikePathBank, batchAnalysisParams)
end

if plotSwitch.selCount
  selCount(spikePathBank, batchAnalysisParams)
end

if plotSwitch.selectivityCurve
  selectivityCurveSelNew(spikePathBank, batchAnalysisParams)
  %selectivityCurveSelNewParComp(spikePathBank, batchAnalysisParams) % Comparing animated and non-animated
end

if plotSwitch.eyeCorr
  allRunEyeCorrelogram(spikePathBank, batchAnalysisParams, figStruct)
end

%% dimensionality reduction

if plotSwitch.dimRed
  % Process data into a more readible format for different dimensionality
  % reduction functions
  dimRedPrep(spikePathBank, batchAnalysisParams)
  
  % Perform a PCA and visualize results
  pcaAndPlot(batchAnalysisParams)
  pcaAndPlotSingleTrial(batchAnalysisParams)
  
  % 
  popDistPlot(batchAnalysisParams)
  
  % Perform a dPCA and visualize results
%   dpcaAndPlot(batchAnalysisParams)
  
  % Perform PSID
%   psidAndPlot(batchAnalysisParams)
  
end

% Prepare things to run in DataHigh
if calcSwitch.dataHigh
  dataHighPrep(spikePathBank, batchAnalysisParams)
end

%% Analyses

if plotSwitch.saccadeAnalysis
  saccadePSTH(spikePathBank, batchAnalysisParams)
end

if plotSwitch.neuralDecodingTB
  batchAnalysisParams.NDTParams.spikePathLoadParams = batchAnalysisParams.spikePathLoadParams;
  batchAnalysisParams.NDTParams.stimParamFile = batchAnalysisParams.stimParamsFilename;
  NeuralDecodingTB(spikePathBank, batchAnalysisParams.NDTParams);
end

if ~any(strcmp(batchAnalysisParams.spikePathLoadParams.files, 'batchAnalyzedData.mat'))
  [spikePathBank, batchAnalysisParams] = pullPSTH(spikePathBank, batchAnalysisParams);
  spikePathFile = batchAnalysisParams.spikePathLoadParams.batchAnalysisOutput;
  save(spikePathFile, 'spikePathBank', 'batchAnalysisParams', '-append')
end

% Combine PSTH across all runs for a particular stimulus.
if plotSwitch.meanPSTH
  meanPSTHLite(spikePathBank, batchAnalysisParams, figStruct);
end

if plotSwitch.subEventPSTH %&& ~exist('meanPSTHStruct','var')
  subEventPSTHStruct = subEventPSTH(spikePathBank, batchAnalysisParams, figStruct);
end

if plotSwitch.rewardPSTH
  rewardPSTHallStack(spikePathBank, batchAnalysisParams)
end

end
function [analysisOutFilename] = runBatchAnalysisLite(spikePathBank, batchAnalysisParams)
%% Overwrite switches with what is currently in file, keeping the
%spikePathLoadParams, which holds most important stuff.
spikePathLoadParamsKeep = batchAnalysisParams.spikePathLoadParams;
batchAnalysisParams = load(fullfile(batchAnalysisParams.outputDir, 'batchAnalysisParams.mat'));
batchAnalysisParams.spikePathLoadParams = spikePathLoadParamsKeep;

% Load basic inputs
runList = spikePathBank.Properties.RowNames;
plotSwitch = batchAnalysisParams.plotSwitch;
calcSwitch = batchAnalysisParams.calcSwitch;
figStruct = batchAnalysisParams.figStruct;
meanPSTHParams = batchAnalysisParams.meanPSTHParams;

%% Pre-Analysis processing

if ~strcmp('stimPresCount',spikePathBank.Properties.VariableNames)
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

% Prepare things to run in DataHigh
if calcSwitch.dataHigh
  dataHighPrep(spikePathBank, batchAnalysisParams)
end

%% Counting
% crossParadigmCheck(spikePathBank, batchAnalysisParams)

if plotSwitch.selCount
  selCount(spikePathBank, batchAnalysisParams)
end

if plotSwitch.selectivityCurve
  selectivityCurveSelNew(spikePathBank, batchAnalysisParams)
end

if plotSwitch.selectivityCounts
%   selectivityCurveCount(spikePathBank, batchAnalysisParams)
  selectivityCurveCountBinned(spikePathBank, batchAnalysisParams)
end

%% Analyses

if plotSwitch.saccadeAnalysis
  saccadePSTH(spikePathBank, batchAnalysisParams)
end

if plotSwitch.neuralDecodingTB
  batchAnalysisParams.NDTParams.spikePathLoadParams = batchAnalysisParams.spikePathLoadParams;
  NeuralDecodingTBLite(spikePathBank, batchAnalysisParams.NDTParams);
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
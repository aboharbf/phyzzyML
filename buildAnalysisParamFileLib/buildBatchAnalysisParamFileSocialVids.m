function [batchAnalysisParamFilename] = buildBatchAnalysisParamFileSocialVids( varargin )
% Generate a file which specifies the parameters for batch analysis. 

[~, machine] = system('hostname');
machine = machine(~isspace(machine));

switch machine
  case 'Alienware_FA'
    analysisDirectory = slashSwap('D:\DataAnalysis\Jan2020');
    outputDir = [analysisDirectory '/batchAnalysis'];
    stimParamsFilename = slashSwap('D:\Onedrive\Lab\ESIN_Ephys_Files\Analysis\phyzzy\stimParamFileLib\StimParamFileSocialVids_Full.mat');   %#ok
  case 'HomeDesktop'
    analysisDirectory = slashSwap('E:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\AnalysisSample');
    outputDir = [analysisDirectory '/batchAnalysis'];
    stimParamsFilename = slashSwap('E:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzy\stimParamFileLib\StimParamFileSocialVids_Full.mat');   %#ok
end

analysisLabel = 'Basic';
preprocessedDataFilenameStem = 'preprocessedData.mat';
analysisParamFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'

saveFig = 1;                
closeFig = 0;               %#ok
exportFig = 0;              %#ok
saveFigData = 0;            %#ok
verbosity = 'INFO';         %other options, 'DEBUG', 'VERBOSE';

%% Switches
calcSwitch.excludeRepeats = 0;
plotSwitch.stimPresCount = 1;         % Figures showing presentation counts across all runs, in development.
plotSwitch.meanPSTH = 1;              % figure showing mean PSTH across all units, MUA, and Unsorted.
plotSwitch.frameFiringRates = 0;      % Figures showing raw, max, mean rates per object depending on viewing during frame.
plotSwitch.novelty = 0;       % Seeing whether 10th percentile values in previous analyses line up with 'novel' stimuli
plotSwitch.slidingWindowANOVA = 0;

%% Parameters
preprocessParams.spikeDataFileName = 'spikeDataBank'; %File ending in .mat, not included to allow for slicing (e.g. 'spikeDataBank_1.mat'...)
%3 variables below are not used, as defining the variables in a load
%command w/ a cell array (or other structures) doesn't work.
preprocessParams.preprocessedVars = {'spikesByEvent','eventIDs','eventCategories','preAlign','postAlign'}; %Variables extracted from preprocessedData.mat
preprocessParams.analyzedVars = {'analysisParamFilename','dateSubject', 'runNum', 'groupLabelsByImage','psthByImage','attendedObjData'}; %Variables extracted from analyzedData.mat
preprocessParams.analysisParamVars = {'psthParams'}; %Variables extracted from analysisParam.mat

cellCountParams.excludePhase2 = 1; % a switch which can be used to remove data from the same neuron collected in subsequent runs. Good for getting accurate counts.
cellCountParams.batchRunxls = fullfile(analysisDirectory,'BatchRunResults.xlsx');                         %Batch analysis xlsx produced by processRunBatch.
cellCountParams.recordingLogxls = 'D:\Onedrive\Lab\ESIN_Ephys_Files\Data\RecordingsMoUpdated.xlsx';  %Used to exclude phase 2 to give accurate unit counts.

meanPSTHParams.outputDir = fullfile(outputDir,'meanPSTH');
meanPSTHParams.stimParamsFilename = stimParamsFilename;
meanPSTHParams.plotHist = 0;
meanPSTHParams.plotTopStim = 0;                 % Only plot stimuli which have been present on at least a certain number of runs.
meanPSTHParams.topStimPresThreshold = 10;       % At least this many stim presentations to be plotted when plotTopStim is on.
meanPSTHParams.broadLabel = 0;                  % Transitions individual stimuli to broad catagory (e.g. chasing).
meanPSTHParams.normalize = 1;                   % Normalizes PSTH values to the recording's fixation period. 1 = Z score.
meanPSTHParams.maxStimOnly = 1;                 % The max value and max index taken from the PSTH is only in the area of the stimulus presentation.
meanPSTHParams.plotLabels = {'chasing','fighting','mounting','grooming','following',...
    'objects','goalDirected','idle','scramble','scene','socialInteraction','animControl','animSocialInteraction','agents','headTurning'}; %If broadLabel is on, all stimuli will have their labels changed to one of the labels in this array.
meanPSTHParams.sortPresCountSort = 1;           % Sorts images based on counts.
meanPSTHParams.fixAlign = 1;                    % For cross catagory comparison lines, shift everything to the mean of the fix period.
meanPSTHParams.topPSTHRunExtract = 3;           % meanPSTH will return a structure of run indices of the top PSTHes by activity (influenced by Z-scoring). This number dictates how many of the top are returned.
meanPSTHParams.type = 'normal';         % options are 'normal', 'baselineSub', 'meanWhite'
meanPSTHParams.psthPre = 800;           % if e.g. +200, then start psth 200ms before trial onset; 
meanPSTHParams.psthImDur = 2800;        % only need to set this for variable length stim runs; else, comes from log file
meanPSTHParams.psthPost = 500;
meanPSTHParams.latency = 0;
meanPSTHParams.movingWin = [200 5];
meanPSTHParams.smoothingWidth = 25;     % psth smoothing width, in ms
meanPSTHParams.Zscore = 1;              % 0: raw PSTH, 1: pre-trial baseline subtraction Z Scored PSTH
meanPSTHParams.errorType = 2;           % chronux convention: 1 is poisfStimson, 2 is trialwise bootstrap, 3 is across trial std for binned spikes, bootstrap for spike times 
meanPSTHParams.errorRangeZ = 1;         % how many standard errors to show
meanPSTHParams.bootstrapSamples = 100;
meanPSTHParams.sortStim = 1;
meanPSTHParams.sortOrder = {'socialInteraction';'goalDirected';'idle';'objects';'scene';'scramble';'allStim'};
meanPSTHParams.psthColormapFilename = 'cocode2.mat'; % a file with one variable, a colormap called 'map'
load(meanPSTHParams.psthColormapFilename);
meanPSTHParams.colormap = map;
meanPSTHParams.tmpFileName = 'tmpStructPrcSigChange.mat';
meanPSTHParams.plotSingleUnitTests = 1; % Avoids running completed plot code.
meanPSTHParams.stimInclude = 2;         % 0 = everything, 1 = Only Animations, 2 = Exclude Animations. 
meanPSTHParams.removeRewardEpoch = 1;   % Removes the reward period activity when generating plots.
meanPSTHParams.plotMeanLine = 0;        % For 'All Chasing' plots, include a additional axis as a line plot.
meanPSTHParams.includeMeanTrace = 0;    % For 'All Chasing' plots, include the mean of all traces at the bottom of the PSTH.
meanPSTHParams.traceCountLabel = 0;     % labels on the catagory specific plots include 'n = X' to highlight trace value.

meanPSTHParams.plotSize = [.7 .7];           % Turns on the 'exportFig' feature of saveFigure, which generates .pngs.
meanPSTHParams.exportFig = 0;           % Turns on the 'exportFig' feature of saveFigure, which generates .pngs.

frameFiringParams.stimParamsFilename = stimParamsFilename;
frameFiringParams.outputDir = fullfile(outputDir,'frameFiring');
frameFiringParams.broadLabelPool = {'chasing','fighting','mounting','grooming','holding','following','observing',...
    'foraging','sitting','objects','goalDirected','idle','scramble','scene','animControl','animSocialInteraction'}; %If broadLabel is on, all stimuli will have their labels changed to one of the labels in this array.
frameFiringParams.broadLabels = 1;
frameFiringParams.useRates = 0;     % collected data per frame can be in rates or spike counts. 1 = Rates, 0 = spikeCounts.
frameFiringParams.delay = 70;       % Hypothetical delay between frame and the activity it causes. deduced from Mean PSTHs.
frameFiringParams.plotRuns = 0;     % Plot Histograms of values across individual runs. 

slidingTestParams.plotTest = 0;   % Plot individual cell pVectors. Saves these to individual files.
slidingTestParams.binSize = 100;
slidingTestParams.binStep = 25;
slidingTestParams.scrambleCount = 100;    % Count of trials to come up with control p values.
slidingTestParams.Omega = 1;       %switch to convert ANOVA curves to Omega curves.
slidingTestParams.target = {'socialInteraction','agents','interaction'};  %Labels which must exist in the stimParamFile associated with the runs. 
slidingTestParams.stimParamFile = stimParamsFilename;
slidingTestParams.outputDir = fullfile(outputDir,'slidingTest');
slidingTestParams.spikeDataFileName = preprocessParams.spikeDataFileName;
slidingTestParams.exportFig = 1; 

noveltyParams.outputDir = fullfile(outputDir,'noveltyAnalysis');


analysisParamFilenameStem = 'batchAnalysisParams.mat';
if ~exist([outputDir '/'])
  mkdir([outputDir '/']);
end
batchAnalysisParamFilename = [outputDir '/' analysisParamFilenameStem];
save(batchAnalysisParamFilename)

end

function swappedString = slashSwap(pathString)
%Swaps direction of slashes to match Unix/Phyzzy, from Windows Path.
  stringParts = split(pathString, '\');
  swappedString = char(join(stringParts, '/'));
end
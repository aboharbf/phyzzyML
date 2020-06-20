function [batchAnalysisParamFilename] = buildBatchAnalysisParamFileSocialVids( varargin )
% Generate a file which specifies the parameters for batch analysis. 

[~, machine] = system('hostname');
machine = machine(~isspace(machine));

switch machine
  case 'Alienware_FA'
    analysisDirectory = slashSwap('D:\DataAnalysis\March2020');
    outputDir = [analysisDirectory '/batchAnalysis'];
    stimParamsFilename = slashSwap('D:\Onedrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\stimParamFileLib\StimParamFileSocialVids_Full.mat');   %#ok
    stimDir = slashSwap('D:\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories');
    eventDataPath = fullfile(stimDir, 'eventData.mat');
    frameMotionDataPath = fullfile(stimDir, 'frameMotion_complete.mat');
    recordingLogxls = 'D:\Onedrive\Lab\ESIN_Ephys_Files\Data\RecordingsMoUpdated.xlsx';
  case 'homeDesktopWork'
    analysisDirectory = 'H:/Analyzed';
    outputDir = [analysisDirectory '/batchAnalysis'];
    stimParamsFilename = slashSwap('C:\Onedrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\stimParamFileLib\StimParamFileSocialVids_Full.mat');   %#ok
    stimDir = slashSwap('C:\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories');
    eventDataPath = fullfile(stimDir, 'eventData.mat');
    frameMotionDataPath = fullfile(stimDir, 'frameMotion_complete.mat');
    recordingLogxls = 'C:\Onedrive\Lab\ESIN_Ephys_Files\Data\RecordingsMoUpdated.xlsx';
    NDTPath = 'C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\NeuralDecodingToolbox';
    NDTAnalysesPath = 'C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses';
end

analysisLabel = 'Basic';
preprocessedDataFilenameStem = 'preprocessedData.mat';
analysisParamFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'

figStruct.saveFig = 1;      % save the figure in its output directory.           
figStruct.closeFig = 0;     % close the figure once it is saved
figStruct.exportFig = 0;    % export figure using export_fig.
figStruct.saveFigData = 0;  % save data with the figure.
figStruct.noOverWrite = 1;      % If a figure is already there, don't make it again.
verbosity = 'INFO';         %other options, 'DEBUG', 'VERBOSE';

%% Switches
calcSwitch.excludeRepeats = 0;
plotSwitch.stimPresCount = 0;         % Figures showing presentation counts across all runs, in development.
plotSwitch.meanPSTH = 0;              % figure showing mean PSTH across all units, MUA, and Unsorted.
plotSwitch.subEventPSTH = 0;          % Analysis of subEvents taking place during stimuli.
plotSwitch.frameFiringRates = 0;      % Figures showing raw, max, mean rates per object depending on viewing during frame.
plotSwitch.novelty = 0;               % 
plotSwitch.slidingWindowANOVA = 0;    % 
plotSwitch.neuralDecodingTB = 1;      % Run the Neural decoding Toolbox

%% Parameters
preprocessParams.spikeDataFileName = 'spikeDataBank'; %File ending in .mat, not included to allow for slicing (e.g. 'spikeDataBank_1.mat'...)
%3 variables below are not used, as defining the variables in a load
%command w/ a cell array (or other structures) doesn't work.
preprocessParams.preprocessedVars = {'spikesByEvent','eventIDs','eventCategories','preAlign','postAlign'}; %Variables extracted from preprocessedData.mat
preprocessParams.analyzedVars = {'analysisParamFilename','dateSubject', 'runNum', 'groupLabelsByImage','psthByImage','attendedObjData'}; %Variables extracted from analyzedData.mat
preprocessParams.analysisParamVars = {'psthParams'}; %Variables extracted from analysisParam.mat

meanPSTHParams.stimInclude = 2;                 % 0 = everything, 1 = Only Animations, 2 = Exclude Animations. 
cellCountParams.excludePhase2 = 0;                                                    % a switch which can be used to remove data from the same neuron collected in subsequent runs. Good for getting accurate counts.
cellCountParams.batchRunxls = fullfile(analysisDirectory,'BatchRunResults.xlsx');     % Batch analysis xlsx produced by processRunBatch.
cellCountParams.recordingLogxls = recordingLogxls;                                    % Used to exclude phase 2 to give accurate unit counts.

meanPSTHParams.outputDir = fullfile(outputDir,'meanPSTH');
meanPSTHParams.stimParamsFilename = stimParamsFilename;
meanPSTHParams.eventData = eventDataPath;
meanPSTHParams.frameMotionDataPath = frameMotionDataPath;
meanPSTHParams.plotHist = 0;
meanPSTHParams.rateThreshold = 0;               % Include only activity with a mean rate of X Hz. 0 for off, over 0 for threshold.

meanPSTHParams.plotTopStim = 0;                 % Only plot stimuli which have been present on at least a certain number of runs.
meanPSTHParams.topStimPresThreshold = 5;        % At least this many stim presentations to be plotted when plotTopStim is on.

meanPSTHParams.broadLabel = 0;                  % Transitions individual stimuli to broad catagory (e.g. chasing).
meanPSTHParams.normalize = 1;                   % Normalizes PSTH values to the recording's fixation period. 1 = Z score.
meanPSTHParams.maxStimOnly = 1;                 % The max value and max index taken from the PSTH is only in the area of the stimulus presentation.
meanPSTHParams.plotLabels = {'chasing','fighting','mounting','grooming','following',...
  'objects','goalDirected','idle','scramble','scene','socialInteraction','animControl','animSocialInteraction','agents','headTurn','allTurn', 'headTurnClassic'}; %If broadLabel is on, all stimuli will have their labels changed to one of the labels in this array.
meanPSTHParams.removeEmpty = 1;         % Used by plotIndex - removes empty columns 
meanPSTHParams.outLogic = 1;            % Used by plotIndex - turns output into logical array of 0s and 1s.
meanPSTHParams.plotLabelSocialInd = [1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0]; %Index for single catagory labels which are social.
meanPSTHParams.socialColor = [240/255 62/255 47/255];
meanPSTHParams.nonSocialColor = [9/255 217/255 107/255];
meanPSTHParams.sortPresCountSort = 1;           % Sorts images based on counts.
meanPSTHParams.fixAlign = 1;                    % For cross catagory comparison lines, shift everything to the mean of the fix period.
meanPSTHParams.topPSTHRunExtract = 3;           % meanPSTH will return a structure of run indices of the top PSTHes by activity (influenced by Z-scoring). This number dictates how many of the top are returned.
meanPSTHParams.type = 'normal';         % options are 'normal', 'baselineSub', 'meanWhite'
meanPSTHParams.latency = 0;
meanPSTHParams.movingWin = [200 5];
meanPSTHParams.smoothingWidth = 25;     % psth smoothing width, in ms
meanPSTHParams.Zscore = 1;              % 0: raw PSTH, 1: pre-trial baseline subtraction Z Scored PSTH
meanPSTHParams.errorType = 2;           % chronux convention: 1 is poisson, 2 is trialwise bootstrap, 3 is across trial std for binned spikes, bootstrap for spike times 
meanPSTHParams.errorRangeZ = 1;         % how many standard errors to show
meanPSTHParams.bootstrapSamples = 100;
meanPSTHParams.sortStim = 1;
meanPSTHParams.sortOrder = {'socialInteraction';'goalDirected';'idle';'objects';'scene';'scramble';'allStim'};
meanPSTHParams.psthColormapFilename = 'cocode2.mat'; % a file with one variable, a colormap called 'map'
load(meanPSTHParams.psthColormapFilename);
meanPSTHParams.colormap = map;
meanPSTHParams.tmpFileName = 'tmpStructPrcSigChange.mat';

meanPSTHParams.removeRewardEpoch = 1;           % Removes the reward period activity when generating plots.
meanPSTHParams.firstXRuns = 0;                % Removes any runs above this number. 0 = Don't remove any.
meanPSTHParams.plotMeanLine = 0;                % For 'All Chasing' plots, include a additional axis as a line plot.
meanPSTHParams.includeMeanTrace = 1;            % For 'All Chasing' plots, include the mean of all traces at the bottom of the PSTH.
meanPSTHParams.traceCountLabel = 1;             % labels on the catagory specific plots include 'n = X' to highlight trace value.

meanPSTHParams.catPSTH = 0;                     % 1.0 - Catagory PSTH Plot - 'All Chasing Stimuli, mean PSTH'
meanPSTHParams.allStimPSTH = 0;                 % 2.0 - All Stimuli means in the same plot.
meanPSTHParams.allRunStimPSTH = 0;              % 3.0 - Stimuli Plot - 'All chasing 1 PSTHs, sorted by...'
meanPSTHParams.lineCatPlot = 0;                 % 4.0 - Line plot with Line per Catagory.
meanPSTHParams.lineBroadCatPlot = 1;            % 5.0 - Means Line plot across broad catagorizations (like Social vs non Social).
meanPSTHParams.splitContrib = 0;                % 5.1 - Mean line plots, split by stimuli.

meanPSTHParams.exportFig = figStruct.exportFig; % Turns on the 'exportFig' feature of saveFigure, which generates .pngs.
meanPSTHParams.plotSizeCatPSTH = [.8 .6];       
meanPSTHParams.plotSizeAllStimPSTH = [.5 1];           
meanPSTHParams.plotSizeAllRunStimPSTH = [1 1];           
meanPSTHParams.plotSizeLineCatPlot = [.5 .6];           
meanPSTHParams.plotSizeLineBroadCatPlot = [.5 .6];           

subEventPSTHParams.outputDir = fullfile(outputDir,'subEventPSTH');
subEventPSTHParams.eventData = eventDataPath;
subEventPSTHParams.stimParamsFilename = stimParamsFilename;
subEventPSTHParams.normalize = 1;                                 % Grab activity from same unit, Z score fixation activity with respect to fixation period activity.
subEventPSTHParams.fixBuffer = 150;                                % normalization acts on the fixation period. Some effects of the fix dot appearance or stimulus onset may be driving neurons away from the true baseline. this number is the millisecond after true fix start, before fix end.
subEventPSTHParams.plotSizeAllRunStimPSTH = [1 1];
subEventPSTHParams.exportFig = figStruct.exportFig;
subEventPSTHParams.saveFig = figStruct.saveFig;
subEventPSTHParams.sparseLabels = 1;                              % In the 'sorted' individual runs, sparse labeling only labels the first entry of that kind in the PSTH.
subEventPSTHParams.frameMotionDataPath = frameMotionDataPath;
subEventPSTHParams.psthPre = 100;                  
subEventPSTHParams.psthImDur = 400;                
subEventPSTHParams.psthPost = 200;             

subEventPSTHParams.allRunStimPSTH = 0;                          % Plot 1 - Individual event PSTHes, stacked
subEventPSTHParams.meanSubEventPSTH = 1;                        % Plot 2 - Mean event PSTHes, line plots
subEventPSTHParams.allRunStimPSTH_extracted = 0;                % Plot 3 - Individual event PSTHes, based on slices extracted from full PSTH data (not collected and avg'd per run).
subEventPSTHParams.meanSubEventPSTH_extracted = 1;              % Plot 4 - Mean event PSTHes, based on slices extracted from full PSTH data (not collected and avg'd per run).

subEventPSTHParams.plotSizeAllRunStimPSTH_extracted = [.5 .6];           
subEventPSTHParams.plotSizeMeanSubEventPSTH_extracted = [.25 .9];               

meanPSTHParams.subEventPSTHParams = subEventPSTHParams; 
meanPSTHParams.subEventPSTHParams.subEventTimes = [200 200];      % psthPre and psthImDur for grabbing events.
meanPSTHParams.subEventPSTHParams.psthPre = 100;                  
meanPSTHParams.subEventPSTHParams.psthImDur = 400;                
meanPSTHParams.subEventPSTHParams.psthPost = 200;                    

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
slidingTestParams.exportFig = figStruct.exportFig;
slidingTestParams.plotSize = [.8 .6];        


% Neural Decoding Toolbox parameters
% Paths
NDTParams.NDTPath = NDTPath;
NDTParams.outputDir = fullfile(outputDir,'NeuralDecodingTB');
NDTParams.rasterDir =  fullfile(NDTParams.outputDir, 'rasterData');
NDTParams.spikeToRasterParams.plotIndParams.stimParamsFilename = stimParamsFilename; % a full path to a phyzzy style stimParamFile.

%  spikeDataBank to rasterData Parameters, should be as inclusive as
%  possible, as changes require deleting previously generated files.
NDTParams.spikeToRasterParams.rasterLabels = {'social', 'agents', 'socialCat'}; % raster labels which are added as fields.
NDTParams.spikeToRasterParams.plotIndParams.plotLabels = {'socialInteraction', 'agents', {'chasing', 'mounting','fighting', 'grooming', 'goalDirected', 'idle', 'objects', 'scene', ...
  'scramble', 'foraging', 'following', 'observing','animControl', 'animSocialInteraction'}}; % a cell array list of labels for inclusion.
NDTParams.spikeToRasterParams.plotIndParams.removeEmpty = 0;
NDTParams.spikeToRasterParams.plotIndParams.outLogic = 0;

% Params for binning of rasters for analyses.
NDTParams.binWidth = 100;
NDTParams.stepSize = 50;

% The following should be specified per analysis
NDTParams.Analyses.SocVNonSoc.label = 'social';                   % should match a field in the raster_labels, and be included in currently generated iteration of the rasterData.
NDTParams.Analyses.SocVNonSoc.num_cv_splits = 100;
NDTParams.Analyses.SocVNonSoc.repeat_each_label_per_split = 1;
NDTParams.Analyses.SocVNonSoc.featPre = [1 0 0];                  % [zscore_nroamlize_FP, select_pvalue_significant_features_FP, select_or_exclude_top_k_features_FP]
NDTParams.Analyses.SocVNonSoc.classifier = [1 0 0];               % [max_correlation_coefficient_CL, poisson_naive_bayes_CL, libsvm_CL], select 1
NDTParams.Analyses.SocVNonSoc.pvalue_threshold = 0.05;            % if select_p FP above, feature needs pVal threshold.
NDTParams.Analyses.SocVNonSoc.k_features = 10;                    % if top_k FP above, feature needs k threshold.
NDTParams.Analyses.SocVNonSoc.label_names_to_use = [];            % the stimuli which all sites must have presented.
NDTParams.Analyses.SocVNonSoc.cross_validator_num_resample = 5;   % Number of times to resample runs for cross validator.
NDTParams.Analyses.SocVNonSoc.plotTitle = 'Social Interaction';   

NDTParams.Analyses.AgentVNonAgent.label = 'agents';
NDTParams.Analyses.AgentVNonAgent.num_cv_splits = 100;
NDTParams.Analyses.AgentVNonAgent.repeat_each_label_per_split = 1;
NDTParams.Analyses.AgentVNonAgent.featPre = [1 0 0];                  % [zscore_nroamlize_FP, select_pvalue_significant_features_FP, select_or_exclude_top_k_features_FP]
NDTParams.Analyses.AgentVNonAgent.classifier = [1 0 0];               % [max_correlation_coefficient_CL, poisson_naive_bayes_CL, libsvm_CL], select 1
NDTParams.Analyses.AgentVNonAgent.pvalue_threshold = 0.05;            % if select_p FP above, feature needs pVal threshold.
NDTParams.Analyses.AgentVNonAgent.k_features = 10;                    % if top_k FP above, feature needs k threshold.
NDTParams.Analyses.AgentVNonAgent.label_names_to_use = [];            % the stimuli which all sites must have presented.
NDTParams.Analyses.AgentVNonAgent.cross_validator_num_resample = 5;   % Number of times to resample runs for cross validator.
NDTParams.Analyses.AgentVNonAgent.plotTitle = 'Agent presence';       

NDTParams.Analyses.socCat.label = 'socialCat';
NDTParams.Analyses.socCat.num_cv_splits = 10;
NDTParams.Analyses.socCat.repeat_each_label_per_split = 1;
NDTParams.Analyses.socCat.featPre = [1 0 0];                      % [zscore_nroamlize_FP, select_pvalue_significant_features_FP, select_or_exclude_top_k_features_FP]
NDTParams.Analyses.socCat.classifier = [1 0 0];                   % [max_correlation_coefficient_CL, poisson_naive_bayes_CL, libsvm_CL], select 1
NDTParams.Analyses.socCat.pvalue_threshold = 0.05;                % if select_p FP above, feature needs pVal threshold.
NDTParams.Analyses.socCat.k_features = 10;                        % if top_k FP above, feature needs k threshold.
NDTParams.Analyses.socCat.label_names_to_use = {'chasing', 'mounting','fighting', 'grooming', 'goalDirected', 'idle'}; % the stimuli which all sites must have presented.
NDTParams.Analyses.socCat.cross_validator_num_resample = 5;       % Number of times to resample runs for cross validator.
NDTParams.Analyses.socCat.plotTitle = 'Social Catagory';

NDTParams.AnalysesDefault.label = 'socialCat';
NDTParams.AnalysesDefault.num_cv_splits = 10;
NDTParams.AnalysesDefault.repeat_each_label_per_split = 1;
NDTParams.AnalysesDefault.featPre = [1 0 0];                      % [zscore_nroamlize_FP, select_pvalue_significant_features_FP, select_or_exclude_top_k_features_FP]
NDTParams.AnalysesDefault.classifier = [1 0 0];                   % [max_correlation_coefficient_CL, poisson_naive_bayes_CL, libsvm_CL], select 1
NDTParams.AnalysesDefault.pvalue_threshold = 0.05;                % if select_p FP above, feature needs pVal threshold.
NDTParams.AnalysesDefault.k_features = 10;                        % if top_k FP above, feature needs k threshold.
NDTParams.AnalysesDefault.label_names_to_use = {'chasing', 'mounting','fighting', 'grooming', 'goalDirected', 'idle'}; % the stimuli which all sites must have presented.
NDTParams.AnalysesDefault.cross_validator_num_resample = 5;       % Number of times to resample runs for cross validator.
NDTParams.AnalysesDefault.plotTitle = '[Analysis File used]'; 

NDTParams.NDTAnalysesPath = NDTAnalysesPath;

% Novelty Analysis
noveltyParams.outputDir = fullfile(outputDir,'noveltyAnalysis');
assert(length(meanPSTHParams.plotLabels) == length(meanPSTHParams.plotLabelSocialInd), 'These two variables must match in length')

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

%#ok<*STRNU> % Gets rid of warning of not use.
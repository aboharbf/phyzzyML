function [batchAnalysisParamFilename] = buildBatchAnalysisLiteParamFileSocialVids( varargin )
% Generate a file which specifies the parameters for batch analysis. 

[~, machine] = system('hostname');
machine = machine(~isspace(machine));

switch machine
  case 'Skytech_FA'
    analysisDirectory = slashSwap('D:\DataAnalysis');
    outputDir = [analysisDirectory '/batchAnalysis'];
    stimParamsFilename = slashSwap('C:\Users\aboha\Onedrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\stimParamFileLib\StimParamFileSocialVids_Full.mat');   %#ok
    stimDir = slashSwap('C:\Users\aboha\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories');
    subEventBatchStructPath = slashSwap(fullfile(analysisDirectory, 'subEventBatchStruct.mat'));
    batchRunxls = fullfile(analysisDirectory,'BatchRunResults.xlsx');
    eventDataPath = fullfile(stimDir, 'eventData.mat');
    frameMotionDataPath = fullfile(stimDir, 'frameMotion_complete.mat');
    recordingLogxls = 'C:\EphysData\Data\RecordingsMoUpdated.xlsx';
    NDTPath = 'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\NeuralDecodingToolbox';
    NDTAnalysesPath = 'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses';
  case 'homeDesktopWork'
    analysisDirectory = 'H:\Analyzed';
    outputDir = [analysisDirectory '/batchAnalysis'];
    stimParamsFilename = slashSwap('C:\Onedrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\stimParamFileLib\StimParamFileSocialVids_Full.mat');   %#ok
    stimDir = slashSwap('C:\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories');
    subEventBatchStructPath = slashSwap(fullfile(analysisDirectory, 'subEventBatchStruct.mat'));
    batchRunxls = fullfile(analysisDirectory,'BatchRunResults.xlsx');
    eventDataPath = fullfile(stimDir, 'eventData.mat');
    frameMotionDataPath = fullfile(stimDir, 'frameMotion_complete.mat');
    recordingLogxls = 'H:\EphysData\Data\RecordingsMoUpdated.xlsx';
    NDTPath = 'D:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\NeuralDecodingToolbox';
    NDTAnalysesPath = 'D:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses';
end

analysisLabel = 'Basic';
preprocessedDataFilenameStem = 'preprocessedData.mat';
analysisParamFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'

figStruct.saveFig = 1;      % save the figure in its output directory.           
figStruct.closeFig = 1;     % close the figure once it is saved
figStruct.exportFig = 0;    % export figure using export_fig.
figStruct.saveFigData = 0;  % save data with the figure.
figStruct.noOverWrite = 1;  % If a figure is already there, don't make it again.
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
spikePathLoadParams.spikePathFileName = 'spikePathBank'; %File ending in .mat, not included to allow for slicing (e.g. 'spikeDataBank_1.mat'...)

preprocessedDataVars = {'spikesByEvent','eventIDs','eventCategories','preAlign','postAlign', 'categoryList'}; %Variables extracted from preprocessedData.mat
analyzedDataVars = {'analysisParamFilename', 'dateSubject', 'runNum', 'groupLabelsByImage','psthByImage','psthErrByImage', ...
                                  'stimStatsTable', 'subEventSigStruct', 'eyeDataStruct','spikesByEventBinned', 'spikesByCategoryBinned', 'psthByCategory', 'psthErrByCategory'}; %Variables extracted from analyzedData.mat
AnalysisParamVars = {'psthParams'}; %Variables extracted from analysisParam.mat


spikePathLoadParams.batchAnalysisOutputName = 'batchAnalyzedData.mat';
spikePathLoadParams.batchAnalysisOutput = fullfile(outputDir, spikePathLoadParams.spikePathFileName);
spikePathLoadParams.files = {'preprocessedData.mat', 'analyzedData.mat', 'AnalysisParams.mat'};
spikePathLoadParams.fileVars = [preprocessedDataVars, analyzedDataVars, AnalysisParamVars];
spikePathLoadParams.fileInds = [ones(size(preprocessedDataVars)), ones(size(analyzedDataVars)) * 2, ones(size(AnalysisParamVars)) * 3];
clear preprocessedDataVars preprocessedDataVars AnalysisParamVars

stimStructParams.spikePathLoadParams = spikePathLoadParams;
stimStructParams.outDir = outputDir;
stimStructParams.figStruct = figStruct;
stimStructParams.eventDataPath = eventDataPath;
stimStructParams.xyStimParams.plotTitle = 'Stimulus Presentation Counts, Sorted by first presentation';
stimStructParams.xyStimParams.XLabel = 'Run Name';
stimStructParams.xyStimparams.YLabel = 'Stimulus Name';
stimStructParams.xyEventparams.plotTitle = 'Events Per Stimulus Count, Sorted by first presentation';
stimStructParams.xyEventparams.XLabel = 'Event Name';
stimStructParams.xyEventparams.YLabel = 'Stimulus Name';

meanPSTHParams.spikePathLoadParams = spikePathLoadParams;
meanPSTHParams.runInclude = 30;                                     % 0 = everything, N = Nth first runs of the stimulus.
meanPSTHParams.stimInclude = 2;                                     % 0 = everything, 1 = Only Animations, 2 = Exclude Animations. 
cellCountParams.excludePhase2 = 0;                                  % a switch which can be used to remove data from the same neuron collected in subsequent runs. Good for getting accurate counts.
cellCountParams.batchRunxls = batchRunxls;                          % Batch analysis xlsx produced by processRunBatch.
cellCountParams.recordingLogxls = recordingLogxls;                  % Used to exclude phase 2 to give accurate unit counts.
cellCountParams.subEventBatchStructPath = subEventBatchStructPath;  % a structure containing info about subEvent selectivity.

meanPSTHParams.outputDir = fullfile(outputDir,'meanPSTH');
meanPSTHParams.stimParamsFilename = stimParamsFilename;
meanPSTHParams.eventData = eventDataPath;
meanPSTHParams.frameMotionDataPath = frameMotionDataPath;
meanPSTHParams.plotHist = 0;
meanPSTHParams.rateThreshold = 0;               % Include only activity with a mean rate of X Hz. 0 for off, over 0 for threshold.
meanPSTHParams.normalize = 0;                   % Normalizes PSTH values to the recording's fixation period. 1 = Z score.

meanPSTHParams.plotTopStim = 0;                 % Only plot stimuli which have been present on at least a certain number of runs.
meanPSTHParams.topStimPresThreshold = 5;        % At least this many stim presentations to be plotted when plotTopStim is on.

meanPSTHParams.broadLabel = 0;                  % Transitions individual stimuli to broad catagory (e.g. chasing).
meanPSTHParams.maxStimOnly = 1;                 % The max value and max index taken from the PSTH is only in the area of the stimulus presentation.
meanPSTHParams.plotLabels = {'chasing','fighting','mounting','grooming','following',...
  'objects','goalDirected','idle','scramble','scene','socialInteraction','animControl','animSocialInteraction','agents','headTurn','allTurn', 'headTurnClassic'}; %If broadLabel is on, all stimuli will have their labels changed to one of the labels in this array.
meanPSTHParams.removeEmpty = 1;                 % Used by plotIndex - removes empty columns 
meanPSTHParams.outLogic = 1;                    % Used by plotIndex - turns output into logical array of 0s and 1s.
meanPSTHParams.plotLabelSocialInd = [1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0]; %Index for single catagory labels which are social.
meanPSTHParams.socialColor = [0/255 162/255 232/255];
meanPSTHParams.nonSocialColor = [244/255 121/255 128/255];
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

% Natural video analysisGroups
meanPSTHParams.analysisGroups.NaturalSocial.socVNonSoc = {'socialInteraction', 'nonSocialInteraction'};

% familiarFace analysisGroups
monkeyIDList = {'Alan', 'Red', 'Otis', 'Pancho', 'Calvin', 'Diego', 'Barney', 'Hobbes'};
monkeyActList = {'IdleC', 'Movement', 'goalDirected', 'IdleS', 'FearGrimace', 'headTurn'};

monkey_act_List = [];
for ii = 1:length(monkeyIDList)
  for jj = 1:length(monkeyActList)
    monkey_act_List = [monkey_act_List; {[monkeyIDList{ii} '_' monkeyActList{jj}]}];
  end
end

meanPSTHParams.analysisGroups.FamilairFace.Identity = monkeyIDList;
meanPSTHParams.analysisGroups.FamilairFace.Action = monkeyActList;
meanPSTHParams.analysisGroups.FamilairFace.stimuli = monkey_act_List;
meanPSTHParams.analysisGroups.FamilairFace.Familair = {'familiar', 'unFamiliar'};

% headTurnCon meanPSTHParams.analysisGroups
meanPSTHParams.analysisGroups.headTurnCon.stimuli = {'chasing1_Turn', 'chasing1_noTurn', 'chasing2_Turn', 'chasing2_noTurn', 'fighting1_noTurn', 'fighting2_Turn', 'fighting2_noTurn', 'grooming1_Turn', 'grooming1_noTurn', 'grooming2_Turn', 'grooming2_noTurn',...
    'idle1_Turn', 'idle1_noTurn', 'idle2_Turn', 'idle2_noTurn', 'mating1_Turn', 'mating1_noTurn', 'mating2_Turn', 'mating2_noTurn', 'objects1_noTurn', 'goalDirected1_Turn', 'goalDirected1_noTurn', 'goalDirected2_Turn'...                      }
    'goalDirected2_noTurn', 'objects2_noTurn', 'scene1_noTurn', 'scene2_noTurn'};
meanPSTHParams.analysisGroups.headTurnCon.socContext = {'socialInteraction', 'nonSocialInteraction'};
meanPSTHParams.analysisGroups.headTurnCon.turnNoTurn = {'noTurn', 'headTurn'};

% headTurnIso analysisGroups
meanPSTHParams.analysisGroups.headTurnIso.turnNoTurn = {'leftFull', 'noTurn', 'rightFull', 'headTurn'};
meanPSTHParams.analysisGroups.headTurnIso.stimuli = {'idle_core_frontCam_Dots_Early', 'idle_core_frontCam_Dots_Late', 'idle_core_frontCam_LowRes_Early', 'idle_core_frontCam_LowRes_Late', 'idle_core_frontCam_Normal_Early', 'idle_core_frontCam_Normal_Late',...          }
        'idle_core_leftCam_Dots_Early', 'idle_core_leftCam_Dots_Late', 'idle_core_leftCam_LowRes_Early', 'idle_core_leftCam_LowRes_Late', 'idle_core_leftCam_Normal_Early', 'idle_core_leftCam_Normal_Late'...
        'idle_core_rightCam_Dots_Early', 'idle_core_rightCam_Dots_Late', 'idle_core_rightCam_LowRes_Early', 'idle_core_rightCam_LowRes_Late', 'idle_core_rightCam_Normal_Early', 'idle_core_rightCam_Normal_Late'...       
        'idle_leftFull_frontCam_Dots_Early', 'idle_leftFull_frontCam_Dots_Late', 'idle_leftFull_frontCam_LowRes_Early', 'idle_leftFull_frontCam_LowRes_Late', 'idle_leftFull_frontCam_Normal_Early'...
        'idle_leftFull_frontCam_Normal_Late', 'idle_leftFull_leftCam_Dots_Early', 'idle_leftFull_leftCam_Dots_Late', 'idle_leftFull_leftCam_LowRes_Early', 'idle_leftFull_leftCam_LowRes_Late'...
        'idle_leftFull_leftCam_Normal_Early', 'idle_leftFull_leftCam_Normal_Late', 'idle_leftFull_rightCam_Dots_Early', 'idle_leftFull_rightCam_Dots_Late', 'idle_leftFull_rightCam_LowRes_Early'...
        'idle_leftFull_rightCam_LowRes_Late', 'idle_leftFull_rightCam_Normal_Early', 'idle_leftFull_rightCam_Normal_Late', 'idle_rightFull_frontCam_Dots_Early', 'idle_rightFull_frontCam_Dots_Late'...
        'idle_rightFull_frontCam_LowRes_Early', 'idle_rightFull_frontCam_LowRes_Late', 'idle_rightFull_frontCam_Normal_Early', 'idle_rightFull_frontCam_Normal_Late', 'idle_rightFull_leftCam_Dots_Early'...
        'idle_rightFull_leftCam_Dots_Late', 'idle_rightFull_leftCam_LowRes_Early', 'idle_rightFull_leftCam_LowRes_Late', 'idle_rightFull_leftCam_Normal_Early', 'idle_rightFull_leftCam_Normal_Late'...
        'idle_rightFull_rightCam_Dots_Early', 'idle_rightFull_rightCam_Dots_Late', 'idle_rightFull_rightCam_LowRes_Early', 'idle_rightFull_rightCam_LowRes_Late', 'idle_rightFull_rightCam_Normal_Early'...
        'idle_rightFull_rightCam_Normal_Late', 'bioMotion_leftFull_frontCam_Normal_Early', 'bioMotion_leftFull_frontCam_Normal_Late', 'bioMotion_leftFull_leftCam_Normal_Early'...
        'bioMotion_leftFull_leftCam_Normal_Late', 'bioMotion_leftFull_rightCam_Normal_Early', 'bioMotion_leftFull_rightCam_Normal_Late'};
meanPSTHParams.analysisGroups.headTurnIso.turnSubj = {'turnToward', 'turnAway'};
meanPSTHParams.analysisGroups.headTurnIso.mesh = {'fullModel_headTurn', 'fullModel_noTurn', 'smoothModel_headTurn', 'smoothModel_noTurn', 'dotModel_headTurn', 'dotModel_noTurn'};

meanPSTHParams.removeRewardEpoch = 1;           % Removes the reward period activity when generating plots.
meanPSTHParams.firstXRuns = 0;                  % Removes any runs above this number. 0 = Don't remove any.
meanPSTHParams.removeFollowing = 1;             % Remove traces related to following due to the ambigious nature of these interactions.
meanPSTHParams.plotMeanLine = 0;                % For 'All Chasing' plots, include a additional axis as a line plot.
meanPSTHParams.includeMeanTrace = 1;            % For 'All Chasing' plots, include the mean of all traces at the bottom of the PSTH.
meanPSTHParams.traceCountLabel = 1;             % labels on the catagory specific plots include 'n = X' to highlight trace value.
meanPSTHParams.addSubEventBars = 0;             % for plot 5.0, add bars underneath for subEvents.

meanPSTHParams.allStimPSTH = 0;                 % 1.0 - All Stimuli means in the same plot.
meanPSTHParams.catPSTH = 0;                     % 2.0 - Catagory PSTH Plot - 'All Chasing Stimuli, mean PSTH'
meanPSTHParams.analysisGroupPSTH = 1;           % 3.0
meanPSTHParams.allRunStimPSTH = 1;              % 3.0 - Stimuli Plot - 'All chasing0001 PSTHs, sorted by...'
meanPSTHParams.lineCatPlot = 1;                 % 4.0 - Line plot with Line per Catagory.
meanPSTHParams.lineBroadCatPlot = 1;            % 5.0 - Means Line plot across broad catagorizations (like Social vs non Social).
meanPSTHParams.splitContrib = 0;                % 5.1 - Mean line plots, split by stimuli.

meanPSTHParams.exportFig = figStruct.exportFig; % Turns on the 'exportFig' feature of saveFigure, which generates .pngs.
meanPSTHParams.plotSizeCatPSTH = [.8 .6];       
meanPSTHParams.plotSizeAllStimPSTH = [.5 1];           
meanPSTHParams.plotSizeAllRunStimPSTH = [1 1];           
meanPSTHParams.plotSizeLineCatPlot = [.6 .7];           
meanPSTHParams.plotSizeLineBroadCatPlot = [.6 .7];           

subEventPSTHParams.spikePathLoadParams = spikePathLoadParams;
subEventPSTHParams.outputDir = fullfile(outputDir,'subEventPSTH');
subEventPSTHParams.eventData = eventDataPath;
subEventPSTHParams.stimParamsFilename = stimParamsFilename;
subEventPSTHParams.normalize = 1;                                   % Grab activity from same unit, Z score fixation activity with respect to fixation period activity.
subEventPSTHParams.fixBuffer = 150;                                 % normalization acts on the fixation period. Some effects of the fix dot appearance or stimulus onset may be driving neurons away from the true baseline. this number is the millisecond after true fix start, before fix end.
subEventPSTHParams.plotSizeAllRunStimPSTH = [1 1];
subEventPSTHParams.exportFig = figStruct.exportFig;
subEventPSTHParams.saveFig = figStruct.saveFig;
subEventPSTHParams.sparseLabels = 1;                              % In the 'sorted' individual runs, sparse labeling only labels the first entry of that kind in the PSTH.
subEventPSTHParams.frameMotionDataPath = frameMotionDataPath;
subEventPSTHParams.psthPre = 100;                  
subEventPSTHParams.psthImDur = 400;                
subEventPSTHParams.psthPost = 200;             

subEventPSTHParams.allRunStimPSTH = 0;                          % Plot 1 - Individual event PSTHes, stacked
subEventPSTHParams.meanSubEventPSTH = 0;                        % Plot 2 - Mean event PSTHes, line plots
subEventPSTHParams.stimPlusEvents_extracted = 0;                % Plot 3 - Show traces of the stimulus on the left (entire trace) + Traces of the subEvent means on the right.  
subEventPSTHParams.eventPsthColorPlots = 0;                     % Plot 3.1 - event PSTH color plots, each individual trace, stacked. Takes length of event into account.
subEventPSTHParams.eventPsthMeanLinePlots = 1;                  % Plot 3.2 - event PSTH mean line plots, shows the full slice of event related activity to the PSTH. ACTIVATES PLOT 3 CODE!
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
slidingTestParams.spikeDataFileName = spikePathLoadParams.spikePathFileName;
slidingTestParams.exportFig = figStruct.exportFig;
slidingTestParams.plotSize = [.8 .6];        

% Neural Decoding Toolbox parameters
% Paths
NDTParams.NDTPath = NDTPath;
NDTParams.outputDir = fullfile(outputDir,'NeuralDecodingTB');
NDTParams.rasterDir =  fullfile(NDTParams.outputDir, 'rasterData');
NDTParams.spikeToRasterParams.plotIndParams.stimParamsFilename = stimParamsFilename; % a full path to a phyzzy style stimParamFile.
NDTParams.spikeToRasterParams.subEventBatchStructPath = subEventBatchStructPath; % a full path to a subEventStruct produced by processRunBatch.
NDTParams.NDTAnalysesPath = NDTAnalysesPath;

%  spikeDataBank to rasterData Parameters, should be as inclusive as
%  possible, as changes require deleting previously generated files.
NDTParams.spikeToRasterParams.rasterLabels = {'social', 'agents', 'socialCat'}; % raster labels which are added as fields.
NDTParams.spikeToRasterParams.plotIndParams.plotLabels = {'socialInteraction', 'agents', {'chasing', 'mounting','fighting', 'grooming', 'goalDirected', 'idle', 'objects', 'scene', ...
  'scramble', 'foraging', 'following', 'observing','animControl', 'animSocialInteraction'}}; % a cell array list of labels for inclusion.
NDTParams.spikeToRasterParams.plotIndParams.removeEmpty = 0;
NDTParams.spikeToRasterParams.plotIndParams.outLogic = 0;

% Params for binning of rasters for analyses.
NDTParams.binWidth = 150;
NDTParams.stepSize = 50;

NDTParams.AnalysesDefault.null_shuffle_count = 20;            % The number of times to randomly shuffle the data to generate a null distribution.
NDTParams.AnalysesDefault.load_data_as_spike_counts = 0;      % loading spike counts is required for poisson_naive_bayes_CL, optional for the rest.
NDTParams.AnalysesDefault.cross_validator_num_resample = 10;  % Number of times to resample runs for cross validator.

NDTParams.figStruct = figStruct;
NDTParams.plotParams = figStruct;
NDTParams.plotParams.points_to_label = [-500, 0, 500, 1000, 1500, 2000, 2500, 3000];
NDTParams.plotParams.points_for_lines = [0, 2800];
NDTParams.plotParams.shift = 800; % prePSTH in the code elsewhere.

NDTParams.p_val_threshold = 0.05;
NDTParams.plot_per_label_acc.p_val_threshold  = NDTParams.p_val_threshold;    % The thrshold to use for determining significant regions
NDTParams.plot_per_label_acc.sig_bar_pos = 'bottom';                          % Determines position of significance bar as top or bottom of plot.
NDTParams.plot_per_label_acc.sig_color = {[232/255 0 5/255]};                 % Determines the color of the bar of the significant regions in 'plot_per_label accuracy'

NDTParams.addTCTSigShading = 1;                                               % Adds shading to the Cross Temporal Decoding matrix
NDTParams.removeSmallPatches = 1;                                             % Removes patches which don't have at least _cutOff elements in them.
NDTParams.removeSmallPatch_cutOff = 6;                                        % The threshold for removing small patches in the Cross temporal decoding matrix.

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
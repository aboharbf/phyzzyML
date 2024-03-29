-function [analysisParamFilename] = buildAnalysisParamFileSocialVids( varargin )    
%buildAnalysisParamFile saves a mat file of parameters, which control the
%behavior of processRun, runAnalysis

% %%%%%%%  USER PARAMETERS, EDIT ROUTINELY %%%%%%%%
% headTurnCon Test - 20201115Mo 001
% naturalSocial Test - 20201117Mo 001

runNum = '001';
dateSubject = '20201115Mo';

assert(~isempty(str2double(runNum)), 'Run number had letters, likely not normal run') %Just for batch runs where unique runs follow unconventional naming scheme.

[~, machine] = system('hostname');
machine = machine(~isspace(machine));

switch machine
  case 'Skytech_FA'
    ephysVolume = slashSwap('D:\EphysData\Data');
    ephysBinVolume = 'D:\EphysDataBin';
    stimulusLogVolume = ephysVolume;
    outputVolume = slashSwap('D:\DataAnalysis');
    phyDir = 'C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\Spike Sorting\npy-matlab';
    spikesDir = 'C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\Spike Sorting\spikes';
    stimParamsFilename = slashSwap('C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\stimParamFileLib\StimParamFileSocialVids_Full.mat');   %#ok
    stimDir = slashSwap('C:\OneDrive\Lab\ESIN_Ephys_Files\Stimuli and Code');
    neuroGLMPath = 'C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\neuroGLM';
  case {'turing.rockefeller.edu','hopper.rockefeller.edu','vera.rockefeller.edu'}
    ephysVolume = '/Freiwald/faboharb/EphysAnalysis/EphysData';
    stimulusLogVolume = ephysVolume;
    stimParamsFilename = 'C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\stimParamFileLib\StimParamFileSocialVids_Full.mat';   %#ok
    outputVolume = '/Freiwald/faboharb/EphysAnalysis/Analyzed';
    stimParamsFilename = '/Freiwald/faboharb/EphysAnalysis/phyzzyML/stimParamFileLib/StimParamFileSocialVids_Full.mat';
    stimDir = '/Freiwald/faboharb/EphysAnalysis/stimDir';
  case 'DataAnalysisPC'
    ephysVolume = slashSwap('\\BlackrockPC\nsp\Data');
    stimulusLogVolume = slashSwap('\\CONTROLLERPC\Monkeylogic Experiments');
    outputVolume = slashSwap('C:\Users\Farid\OneDrive\Lab\ESIN_Ephys_Files\Analysis\Analyzed_Rig');
    stimParamsFilename = slashSwap('C:\Users\Farid\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\stimParamFileLib\StimParamFileSocialVids_Full.mat');                  %#ok
    stimDir = slashSwap('C:\Users\Farid\OneDrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories'); 
end

analysisLabel = 'Basic';
preprocessedDataFilenameStem = 'preprocessedData.mat';
analysisParamFilenameStem = 'AnalysisParams.mat'; %change name should be 'leaf'

figStruct = struct();
figStruct.saveFig = 1;                
figStruct.closeFig = 1;               
figStruct.exportFig = 0;             
figStruct.saveFigData = 0;            
figStruct.figPos = [0 0 .6 0.7];      % Normalized units for figure position

savePreprocessed = 1;             % #ok
verbosity = 'INFO';               % other options, 'DEBUG', 'VERBOSE';

%% Plot switches
plotSwitch.pupilDilation = 0;               % plots image which illustrates continuous values for pupil dilation. 
plotSwitch.eyeStatsAnalysis = 1;            % use ClusterFix to generate a vector characterizing eye movements. used by subEventAnalysis, and subsequently selDir. 
plotSwitch.attendedObject = 1;              % Vectors to distinguish where subject is looking. Required for prefImRasterColorCoded.
plotSwitch.imageEyeMap = 0;                 
plotSwitch.eyeCorrelogram = 0;              % Eye Gram
plotSwitch.orbitPosSel = 0;
plotSwitch.rampingAnalysis = 1;
plotSwitch.latencyAnalysis = 1;

plotSwitch.subEventAnalysis = 1;            % plot traces comparing activity surrounding an event (defined in eventData, generated w/ eventDetectionApp), vs null.
plotSwitch.saccadeRasterPSTH = 1;
plotSwitch.neuroGLM = 0;                    % implements neuroGLM package from jpillow lab/github.

plotSwitch.eyeStimOverlay = 0;              % Visualize eye traces on stimuli. Depends greatly on switches below (may just be used for certain variables).
plotSwitch.spikePupilCorr = 0;              % see the correlation between single trial PSTHes and pupil size.

plotSwitch.clusterOnEyePaths = 0;           % Resort spikes based on distinct eye paths, make "New events".
plotSwitch.stimPSTHoverlay = 0;             % grabs stimuli and plots the PSTH underneath.

plotSwitch.imagePsth = 1;                   % a PSTH for every stimulus in the file.
plotSwitch.categoryPsth = 0;                % a PSTH for every category represented in the file and the categoryList of stimParamFile.

plotSwitch.analysisGroupsPsth = 0;          % a PSTH for every set of analysisGroups defined below.
plotSwitch.prefImRaster = 0;                % Raster, Not color coded.
plotSwitch.prefImRasterColorCoded = 0;      % Raster, uses info from attendedObj switch. 1 is colored spikes, 2 is colored background, 3 is Saccade Image, 4 is pupil img.
plotSwitch.topStimToPlot = 5;               % Determines number of stimuli for which rasters are plotted.
plotSwitch.prefImRasterEvokedOverlay = 0;   % Produces images for MUA and Unsorted even if the same. Relies on sometihng in CatPSTH.
plotSwitch.prefImRasterAverageEvokedOverlay = 0;
plotSwitch.prefImMultiChRasterEvokedOverlay = 0;
plotSwitch.imageTuningSorted = 0;           % Barplot per image, Required for stimPSTHoverlay, sigStruct
plotSwitch.stimPrefBarPlot = 0;             % Per event bar graph.
plotSwitch.stimPrefBarPlotEarly = 0;
plotSwitch.stimPrefBarPlotLate = 0;
plotSwitch.tuningCurves = 0;
plotSwitch.tuningCurvesEarly = 0;
plotSwitch.tuningCurvesLate = 0;
plotSwitch.RF = 0;
plotSwitch.rfEarly = 0;
plotSwitch.rfLate = 0;
plotSwitch.latencyRF = 0;
plotSwitch.evokedPowerRF = 0;
plotSwitch.evokedPsthMuaMultiCh = 0;
plotSwitch.evokedByCategory = 0;
plotSwitch.analogInByItem = 0;
plotSwitch.analogInDerivativesByItem = 0;
plotSwitch.colorPsthEvoked = 0;
plotSwitch.linePsthEvoked = 0;
plotSwitch.runSummary = 0;
plotSwitch.runSummaryImMeanSub = 0;
plotSwitch.runSummaryImMeanSubDiv = 0;
plotSwitch.lfpPowerMuaScatter = 0; 
plotSwitch.lfpPeakToPeakMuaScatter = 0;
plotSwitch.lfpLatencyMuaLatency = 0;
plotSwitch.lfpPowerAcrossChannels = 0;
plotSwitch.lfpPeakToPeakAcrossChannels = 0;
plotSwitch.lfpLatencyShiftAcrossChannels = 0;
plotSwitch.singleTrialLfpByCategory = 0;
plotSwitch.lfpSpectraByCategory = 0; % LFP Comparison 
plotSwitch.lfpAutocorrelTfByItem = 0;
plotSwitch.lfpAutocorrelByItem = 0;
plotSwitch.spikeSpectraByCategory = 0;
plotSwitch.spikeAutocorrelTfByItem = 0;
plotSwitch.spikeAutocorrelByItem = 0;
plotSwitch.spikeSpectraTfByImage = 0;
plotSwitch.lfpSpectraTfByImage = 0;
plotSwitch.couplingPhasesUnwrapped = 0;
plotSwitch.couplingPhasesAsOffsets = 0;
plotSwitch.couplingPhasesPolar = 0;
plotSwitch.tfSpectraByCategory = 1; %Do I want this?
plotSwitch.tfErrs = 0;           

%% Calc switches
calcSwitch.categoryPSTH = 1;
calcSwitch.imagePSTH = 1;
calcSwitch.faceSelectIndex = 0;
calcSwitch.faceSelectIndexEarly = 0;
calcSwitch.faceSelectIndexLate = 0;
calcSwitch.inducedTrialMagnitudeCorrection = 0;
calcSwitch.evokedSpectra = 0;
calcSwitch.inducedSpectra = 0;
calcSwitch.evokedImageTF = 0;
calcSwitch.inducedImageTF = 0;
calcSwitch.evokedCatSpikeTF = 0; %Required for one of the above plot switches to actually produce the figure, but crashes @ "spikesByItemBinned = spikesByCategoryBinned;" in the 2k lines.
calcSwitch.inducedCatSpikeTF = 0;
calcSwitch.evokedCatLFPTF = 0; %Required for one of the above plot switches to actually produce the figure, but crashes @ "spikesByItemBinned = spikesByCategoryBinned;" in the 2k lines.
calcSwitch.inducedCatLFPTF = 0;
calcSwitch.evokedCoupling = 0;
calcSwitch.inducedCoupling = 0;
calcSwitch.meanEvokedTF = 0;
calcSwitch.trialMeanSpectra = 0;
calcSwitch.coherenceByCategory = 0;
calcSwitch.spikeTimes = 0;
calcSwitch.useJacknife = 0;      

%% Parameters
channels2Read = 1:32;

% parameters preprocessSpikes and preprocessLFP, see functions for details
ephysParams.needLFP = 1;
ephysParams.needSpikes = 1;
autoChannelDetect = 1;
ephysParams.spikeChannels = channels2Read; %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
ephysParams.lfpChannels = channels2Read;
ephysParams.channelNames = arrayfun(@(x) {sprintf('Ch%d', x)}, channels2Read);
ephysParams.lfpChannelScaleBy = repmat(8191/32764, [length(channels2Read), 1]); %converts raw values to microvolts
ephysParams.offlineSorted = 0;        % Checks for a '*.xls' Structure in the folder, with resorted spikes.
ephysParams.spikeSort = 1;            % Using other spike sorting packages. 1 = waveClus, 2 = previously sorted on Phy.
ephysParams.paramHandle = @set_parameters_batchFeb2020; %Function which produces param struct for wave_clus. in wave_clus folder.
ephysParams.waveClusReclass = 0;      % Reclassify clusters (as defined by mean waveform proximity to threshold) to MUA.
ephysParams.waveClusMUAThreshold = 1.25; %scaling for reclassification of clusters as MUA. 1 = 0 scaling = no reclassification of clusters.
ephysParams.waveClusProjPlot = 1;     % Plots all the clusters in higher dimensional space (defined in type and number in wave_clus parameters).
ephysParams.waveClusClear = 1;        % 1 deletes all files associated with analysis (leaves processed NSX files), 2 deletes the entire associated '_parsed' folder.
ephysParams.commonRef = [0];          % not yet implemented; will allow software re-refrence across headstages
ephysParams.stimulationChannels = []; % not yet implemented; will read stimulation currents recorded at headstage
ephysParams.cPtCal = 1/30;            % conversion from spike sample indices to timestep of decimated LFP
ephysParams.decimateFactorPass1 = 6;  % note: product of the two decimate factors should be 30, if 1 khz samples desired
ephysParams.decimateFactorPass2 = 5;
ephysParams.samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to LFP (should be raw rate/productOfDecimateFactors)
%note: use Blackrock indexing for unitsToUnsort and unitsToDiscard, so unsorted is 0, first defined unit is 1, etc.
ephysParams.unitsToUnsort = cell(length(channels2Read),1); %these units will be re-grouped with u0
ephysParams.unitsToDiscard = cell(length(channels2Read),1); %these units will be considered noise and discarded
ephysParams.spikeWaveformPca = 0;
ephysParams.plotSpikeWaveforms = 0; %0, 1 to build then close, 2 to build and leave open
ephysParams.spikeWaveformsColors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483]];
ephysParams.shiftSpikeWaveforms = 0;
% see http://www.mathworks.com/help/signal/examples/filter-design-gallery.html
hp1Hz = designfilt('highpassiir', 'FilterOrder',8,'PassbandFrequency',1, ...
  'StopbandAttenuation',100,'PassbandRipple',0.5,'SampleRate',1000);     %#ok
% note: with these specifications, returns a 48th order butterworth filter
butter1Hz200Hz_v1 = designfilt('bandpassiir','DesignMethod','butter','PassbandFrequency1',1,'PassbandFrequency2',200,...
  'SampleRate',1000,'MatchExactly','passband','StopbandFrequency1',0.67,'StopbandFrequency2',250);
[tmp1,tmp2] = butter(4,[1/500,200/500]);
butter1Hz200Hz_v2 = [tmp1,tmp2];        %#ok
clear('tmp1','tmp2')
ephysParams.filter = butter1Hz200Hz_v1; % if filtering desired, ephysFilter is a digitalFilter
ephysParams.plotFilterResult = 0; 
ephysParams.outDir = sprintf('%s/%s/%s/%s/',outputVolume, dateSubject,analysisLabel,runNum);
ephysParams.saveFig = figStruct.saveFig;
% ephysParams.phyParams.phyPath = phyDir;
% ephysParams.phyParams.spikesDir = spikesDir;
% ephysParams.phyParams.ephysBinVolume = ephysBinVolume;

% parameters preprocessAnalogIn, see function for details
analogInParams.needAnalogIn = 1;
analogInParams.analogInChannels = [129,130,137]; 
analogInParams.channelNames = {'eyeX','eyeY','pupil'};
analogInParams.channelUnits = {'dva','dva','au'};
analogInParams.analogInChannelScaleBy = [5/32764 5/32764 5/32764]; %converts raw values to volts
analogInParams.decimateFactorPass1 = 1; 
analogInParams.decimateFactorPass2 = 1;
analogInParams.filterPad = 0;
analogInParams.samPerMS = 1; %THIS IS AFTER DECIMATION, and applies to analogIn (should be raw rate/productOfDecimateFactors)
% see http://www.mathworks.com/help/signal/examples/filter-design-gallery.html
butter200Hz_v1 = designfilt('lowpassiir', 'PassbandFrequency', 120, 'StopbandFrequency', 480, 'PassbandRipple', 1,...
  'StopbandAttenuation', 60, 'SampleRate', 1000, 'MatchExactly', 'passband');  %this returns a 3rd order iir filter
nyqfrq = 1000 ./ 2;
%lowPass30Hz = fir2(60,[0,30./nyqfrq,30./nyqfrq,1],[1,1,0,0]); %60th order, 30 Hz low pass filter
lowPass30Hz = designfilt('lowpassfir', 'FilterOrder', 3, 'PassbandFrequency', 30, 'StopbandFrequency', 31, 'SampleRate', 1000);
analogInParams.filters = {0, 0, 0}; %{lowPass30Hz, lowPass30Hz, 0};%{butter200Hz_v1;butter200Hz_v1;butter200Hz_v1}; %filter channel i if filters{i} is digital filter or 1x2 numeric array
analogInParams.plotFilterResult = 1;

photodiodeParams.needStrobe = 1;
photodiodeParams.inputDataType = 'blackrockFilename';
photodiodeParams.dataChannel = 131;
photodiodeParams.dataLoader = [];         % Incase you're using something besides a raw array or a blackrock file.
photodiodeParams.peakTimeOffset = 0;      % this is the offset, in ms, of the peak from the event it's coupled to (note: will be subtracted, so should be > 0 if peak follows event, type: numeric)
photodiodeParams.strobeTroughs = 1;       % Strobe causes troughs.
photodiodeParams.peakFreq = 85;           % approximate number of peaks per second
photodiodeParams.levelCalibType = 'auto'; % 'hardcode', 'hardcodeAndPlot', 'hardcodeAndCheck', 'auto', 'autoAndPlot','autoAndCheck', 'manual'
photodiodeParams.numLevels = 1;
photodiodeParams.saveCalibFile = 0;
photodiodeParams.centerCornerOffset = 5.3;
photodiodeParams.stimulusTriggerChannel = [];
photodiodeParams.peaksToPlot = 50;        % number of peaks to show in calibration plots
photodiodeParams.cleanPeaks = 1;
photodiodeParams.useRisingEdge = 0;       % 0 = peaks, 1 = rising edge.
photodiodeParams.numLevels = 3;           % Number of levels the strobe is acting at.
photodiodeParams.levelHigh = 5350;
photodiodeParams.levelMid = 4000;
photodiodeParams.levelLow = 1250;
photodiodeParams.checkHighLowAlternation = 0;
photodiodeParams.minPeakNumInLevel = 30;
%photodiodeParams.samplingFreq = 30000; %Defined in Blackrock
photodiodeParams.hardcodeFromFile = 0; %(only needed for hardcode calib types)
photodiodeParams.saveCalibFile = 0;
photodiodeParams.saveFigures = 0;
photodiodeParams.displayStats = 1;
lineNoiseTriggerParams.needStrobe = 0;

%     - cleanPeaks: if true, keep peaks only if the signal dropped below the low-peak threshold between them;
%                   for a series of peaks without a sub-threshold trough, keep only the highest (type: logical)
%     - useRisingEdge: define threshold using peaks, then use rising edge
%       threshold crossing for trigger times; optional, default 0 (type: logical)
%     for trigger times
%     - checkHighLowAlternation: (optional, default zero, only used if numLevels == 2)
%     - numLevels: (type: int)
%     - minPeakNumInLevel: min number of peaks required to define a level (type: int)
%     - samplingFreq: not needed if inputDataType is blackrockFilename, or
%                     a custom handle that sets sampling freq (numeric)
%     - peaksToPlot: number of peaks to show in calibration plots; suggested value > 10 (int)
%     - peakFreq: approximate number of peaks per second
%     - hardcodeFromFile: (type: logical) (only needed for hardcode calib types)
%     - saveCalibFile: (type: logical)
%     - inputCalibrationFile: filename ending in .mat, contains field lowThreshold (numeric)
%                             if numLevels >= 2, also contains field highThreshold (numeric)
%                             if numLevels == 3, also contains field midThreshold (numeric)
%                             (required only if hardcodeFromFile): (type: string)    
%     - outputCalibrationFile: filename ending in .mat (required only if saveCalibFile): (type:string)
%     - saveFigures: (type: logical)
%     - calibFigFname: figure filename, without .fig stem; required only if saveFigures (type:string)
%     - triggersFigFname: figure filename, without .fig stem; required only if saveFigures (type:string)
%     - dateSubject: MMDDYYNAME or similar; required only if saveFigures (type: string)
%     - runNum: eg 002; required only if saveFigures (type: string)
%     - outDir: directory for figure and calib. output; include the final '/' in the path (type: string) 
%               (required only if saveFigures or saveCalibFile is 1)
%

% parameters preprocessLogFile, see function for details
stimSyncParams.logProcessor = @preprocessLogFileMonkeyLogic;
stimSyncParams.tryPreparsedLogFile = 1;
stimSyncParams.keepTriggersAndSubTriggers = 0;
stimSyncParams.subTriggerArrayFilenames = {'socialSceneConcatSubTriggers.mat'};
stimSyncParams.usePhotodiode = 0;
stimSyncParams.plotFitShifts = 0;   % Plots a histogram of differences b/t trial start times on Blackrock vs on MonkeyLogic, following a shift to the Blackrock clock.
stimSyncParams.stimDir = stimDir;
stimSyncParams.outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum); 
%
eyeCalParams.needEyeCal = 1;
eyeCalParams.method = 'monkeyLogic'; %'zeroEachFixation', 'monkeyLogic'
eyeCalParams.monkeyLogicShift = 1;
eyeCalParams.makePlots = 0;
eyeCalParams.eyeXChannelInd = 1;
eyeCalParams.eyeYChannelInd = 2;
eyeCalParams.eyeDChannelInd = 3;
eyeCalParams.samplingRate = 120; %Sampling rate in Hz of equipment
eyeCalParams.gainX = 0;
eyeCalParams.gainY = 0;
eyeCalParams.flipX = 0;
eyeCalParams.flipY = 0; 
eyeCalParams.offsetX = 0;
eyeCalParams.offsetY = 0; 
eyeCalParams.calFile = ''; %note: needed only when method = fromFile
eyeCalParams.fixOutLag = 10; 
eyeCalParams.minFixZeroTime = 1000; 
eyeCalParams.analogInParams = analogInParams; %#ok<*STRNU>

accelParams.needAccelCal = 0;
accelParams.accelChannels = {[4;5;6]};
accelParams.channelGains = {[1/.666 1/.666 1/.666]};
accelParams.calMethods = {'hardcode'}; %other option is 'calFile'; calibration method
accelParams.calFiles = {''}; %if method is 'calFile', an ns2 filename

% parameters for excludeStimuli, see function for details
%excludeStimParams.excludeTrials = @excludeTrialsMonkeyLogic;
excludeStimParams.excludeProcessor = @excludeTrialsMonkeyLogic;
excludeStimParams.needExcludeTrials = 1;
excludeStimParams.excludeFailed = 1;
excludeStimParams.excludeAfterFailed = 0;
excludeStimParams.frameDropThreshold = 5;
excludeStimParams.fixPre = 100; %ms
excludeStimParams.fixPost = 100; %ms
excludeStimParams.flashPre = 0;  %ms
excludeStimParams.flashPost = 0; %ms
excludeStimParams.juicePre = 0; % optional, ms
excludeStimParams.juicePost = 0; % optional, ms
excludeStimParams.DEBUG = 0; % makes exclusion criterion plots if true
% additional optional excludeStimParams: accel1, accel2, minStimDur (ms)

% TW=3 with T=.2, then W = 15 Hz (5 tapers)
% TW=1.5 with T=.1, then W = 15 Hz (2 tapers)
% TW = 1.5 with T=.2, then W = 7.5 Hz (2 tapers)
chronuxParams.tapers = [3 5]; %[3 5] is chronux default; 
chronuxParams.pad = 1;
chronuxParams.fs = 1;
chronuxParams.trialave = 1;
chronuxParams.err = [1 .05];  %note: first entry will be automatically switched to 2 if calcSwitch.useJacknife == 1
chronuxParams.fpass = [0 .1]; 
tfParams.movingWin = [300 5]; 
tfParams.specgramRowAve = 0;

psthParams.type = 'normal'; %options are 'normal', 'baselineSub', 'meanWhite'
psthParams.psthPre = 1200; % if e.g. +200, then start psth 200ms before trial onset; 
psthParams.psthImDur = 2800;  % only need to set this for variable length stim runs; else, comes from log file
psthParams.psthPost = 800;
psthParams.latency = 0;
psthParams.movingWin = tfParams.movingWin;
psthParams.smoothingWidth = 40;  %psth smoothing width, in ms
psthParams.Zscore = 0;  % 0: raw PSTH, 1: pre-trial baseline subtraction Z Scored PSTH, 2: whole trial baseline subtracted Z Scored PSTH
psthParams.errorType = 3; % chronux convention: 1 is poisfStimson, 2 is trialwise bootstrap, 3 is across trial std for binned spikes, bootstrap for spike times 
psthParams.errorRangeZ = 1; % how many standard errors to show
psthParams.bootstrapSamples = 10;
psthParams.sortStim = 1;

% sortOrder is referenced via .(taskData.paradigm).sortOrder. 
psthParams.headTurnIso.sortOrder = {'leftFull', 'noTurn', 'rightFull', 'headTurn'};
psthParams.familiarFace.sortOrder = {'Alan', 'Red', 'Otis', 'Pancho', 'Calvin', 'Diego', 'Barney', 'Hobbes', 'Empty'};
psthParams.headTurnCon.sortOrder = {'socialInteraction';'goalDirected';'idle';'objects';'scene'};
psthParams.naturalSocial.sortOrder = {'socialInteraction';'goalDirected';'idle';'objects';'scene';'scramble';'headTurn';'subEvents'};
psthParams.psthColormapFilename = 'cocode2.mat'; % a file with one variable, a colormap called 'map'
load(psthParams.psthColormapFilename);
psthParams.colormap = map;

eyeStatsParams.psthPre = psthParams.psthPre;
eyeStatsParams.psthImDur = psthParams.psthImDur;
eyeStatsParams.stimDir = stimDir;
eyeStatsParams.lfpPaddedBy = tfParams.movingWin(1)/2;
eyeStatsParams.outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
eyeStatsParams.clusterFixLPFilterIn = 25;                 % Filter used for saccade detection internally in ClusterFix.

eyeStimOverlayParams.stimDir = stimDir;
eyeStimOverlayParams.spikeOverlay = 1;         % Adding 3D Gaussian blurred spiking activity. This function is UNSTABLE due to MASSIVE MEMORY USAGE.
eyeStimOverlayParams.videoOutput = false;      % Switch for producing video output. If off, videos are not produced, but down sampled eye signal and other structures are still produced.
eyeStimOverlayParams.shapeOverlay = 1;         % Switch for shape overlay, visualizing frame motion data if available.
eyeStimOverlayParams.trialNumberOverlay = 1;   % Switch for trial number overlay on eye signal.
eyeStimOverlayParams.colorCodedTrace = 1;      % Trace changes colors corresponding to what is being looked at. 1 = per target of gaze, 2 = Saccade vs fixation.
eyeStimOverlayParams.computeEyeSpikeMat = 0;   % Compute the eyeSpikeMat. It is massive, and will likely exhaust system memory.
eyeStimOverlayParams.spikeOverlay = 0;         % Overlay spikes onto the video
eyeStimOverlayParams.frameCountOverlay = 1;    % Places the frame number in the lower left.
eyeStimOverlayParams.eyeSig.shape = 'Circle';  % How to visualize the eye signal
eyeStimOverlayParams.outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
eyeStimOverlayParams.eyeSig.color = {{[1. 0. 0.];[0 .4 1.];[.1 .8 .1];[.1 .8 .1];[0 0 0];[1 .4 0]; ...
  [.7 0 0];[0 0 .7];[0 .5 0];[0 .5 0];[0 0 0];[1 .6 0]; [1 1 1]}; {[0 .653 .91]; [1 0 0]; [1 1 0]}}; %Faces, Bodies, Hands, Obj, bkg, then Fix v Saccade v Blink
  
subEventAnalysisParams.outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
subEventAnalysisParams.preAlign = 300;
subEventAnalysisParams.postAlign = 400;
subEventAnalysisParams.nullAllStim = 1;
subEventAnalysisParams.RewardEvent = 1;
subEventAnalysisParams.rewardAntTime = 200;       % A time to look prior to the mean reward time. This window is compared to the period after reward delivery. PreReward > Post Reward = Reward Anticipation Neuron.
subEventAnalysisParams.FixEvent = 1;
subEventAnalysisParams.StimOnOffEvent = 1;
subEventAnalysisParams.nullSampleMult = 10;       % For every n stimuli with the event, sample all other stimuli this number of times for the null distribution.
subEventAnalysisParams.psthParams = psthParams;
subEventAnalysisParams.psthParams.psthPre = 200;
subEventAnalysisParams.psthParams.psthImDur = 200;
subEventAnalysisParams.psthParams.psthPost = 0;
subEventAnalysisParams.psthParams.smoothingWidth = 10;
subEventAnalysisParams.psthParams.movingWin = psthParams.movingWin;
subEventAnalysisParams.stimPlotParams.psthPre = psthParams.psthPre;
subEventAnalysisParams.stimPlotParams.psthImDur = psthParams.psthImDur;
subEventAnalysisParams.stimPlotParams.psthPost = psthParams.psthPost;
subEventAnalysisParams.stimDir = stimDir;
subEventAnalysisParams.spikeTimes = calcSwitch.spikeTimes;
subEventAnalysisParams.genPlots = 1;                                    % Asks if you want to generate plots.
subEventAnalysisParams.specSubEvent = 0;                                % Analyze individual instances of subEvents in eventData.
subEventAnalysisParams.preSaccOffset = 200;                             % Offset to use when determining pre-saccade period.
subEventAnalysisParams.possibleEvents = {'headTurn_right', 'headTurn_left', 'bodyTurn', 'eyeContact', 'turnToward', 'turnAway',...
                                         'pre-saccades', 'saccades', 'pre-saccadesNonStim', 'saccadesNonStim', 'blinks', 'fixation',...
                                         'reward', 'reward_soc', 'rewardAbsent', 'rewardAnt', 'stimOnset', 'stimOffset'};
subEventAnalysisParams.possibleEventsPlotNames = {'Head turn, right', 'Head turn, left', 'Body turn', 'Eye contact', 'turnToward', 'turnAway',...
                                         'Saccade', 'Saccade', 'Saccade', 'Saccade', 'Blink', 'Fixation',...
                                         'Reward', 'Reward', 'Reward', 'Reward', 'Stimulus Onset', 'Stimulus Offset'};
subEventAnalysisParams.testPeriodPerEvent = [[0 200]; [0 200]; [0 200]; [0 200]; [0 200]; [0 200];...
                                             [-200 0]; [0 100]; [-200 0]; [0 100]; [-50 150]; [0 200];...
                                             [0 200]; [0 200]; [0 200]; [-subEventAnalysisParams.rewardAntTime 0]; [0 200]; [0 200]];
subEventAnalysisParams.nonParametric = 0;                                 % Use non parametric test.

saccadeRasterParams.psthParams = psthParams;
saccadeRasterParams.preAlign = 400;
saccadeRasterParams.postAlign = 400;
saccadeRasterParams.psthParams.psthPre = 200;
saccadeRasterParams.psthParams.psthImDur = 100;
saccadeRasterParams.psthParams.psthPost = 0;  
saccadeRasterParams.spikeTimes = 0;
saccadeRasterParams.psthParams.movingWin(1) = 240;
saccadeRasterParams.smoothingWidth = 20;  %psth smoothing width, in ms

% Variables for creating saccade rasters and PSTHes.
saccadeStackParams.preEventTime = 300;
saccadeStackParams.postEventTime = 300;
saccadeStackParams.preStimTime = saccadeStackParams.preEventTime;
saccadeStackParams.postStimTime = 100;
saccadeStackParams.smoothingWidth = 20;

correlParams.maxShift = []; % a number, or empty
correlParams.matchTimeRanges = 1;
correlParams.timeDifferenceBound = [0,200];
correlParams.normalize = 1;
correlParams.useJacknife = 0;
correlParams.jacknifeDraws = 100;
switch machine
  case 'laptop'
    correlParams.jacknifeParallelWorkers = 0;    
  case 'hopper'
    correlParams.jacknifeParallelWorkers = 0;   
  case 'turing'
    correlParams.jacknifeParallelWorkers = 20;    
end
spikeCorrelSmoothingWidth = 5; %ms
filterPoints = -20*spikeCorrelSmoothingWidth:20*spikeCorrelSmoothingWidth;
smoothingFilter = exp(-1*filterPoints.^2/(2*spikeCorrelSmoothingWidth^2));
correlParams.smoothingFilter = smoothingFilter/sum(smoothingFilter); %#ok
%
lfpAlignParams.samPerMS = 1; % because this is after decimation
lfpAlignParams.msPreAlign = psthParams.psthPre + tfParams.movingWin(1)/2; 
lfpAlignParams.msPostAlign = psthParams.psthImDur + psthParams.psthPost + tfParams.movingWin(1)/2;
%
spikeAlignParams.preAlign = psthParams.psthPre + 3*psthParams.smoothingWidth;
spikeAlignParams.postAlign = psthParams.psthImDur + psthParams.psthPost+3*psthParams.smoothingWidth;   %#ok
% for lfps, constrain first and (optional) last [n m] samples to 0 mean
useDCSUB = 0;
if useDCSUB
  %lfpAlignParams.DCSUB_SAM = [lfpAlignParams.msPreAlign, lfpAlignParams.msPreAlign+10; 0, 0 ]; % 0-th order 
  lfpAlignParams.DCSUB_SAM = [lfpAlignParams.msPreAlign, lfpAlignParams.msPreAlign+10;lfpAlignParams.msPreAlign, lfpAlignParams.msPreAlign+10 ]; % 1st order 
else
  lfpAlignParams.DCSUB_SAM = 0;
end
% firing rate calculation epochs. Can provide either time (ms from stim onset),
% or function handle, which will receive the minimum stimulus duration in
% the run as an input.
% frEpochsCell = {{60, @(stimDur) stimDur+60};...
%                 {-800, 60}; ...
%                 {100, 1100}}; %#ok
%               
frEpochsCell = {{60, @(stimDur) stimDur+60}};...
%                 {-800, 60}; ...
%                 {@(stimDur) stimDur+60, @(stimDur) stimDur+460}}; %#ok
epochLabels = {'Presentation'};%,'Fixation','Reward'};

% epochStats Params, more are below analysisGroups
preFix = [-500 0];      % Calculated on a fixation dot aligned activity trace.
Fix = [-800 0];
stimEarly = [0 500];
stimLate = [500 psthParams.psthImDur];
reward = [psthParams.psthImDur psthParams.psthImDur + 350];

if psthParams.psthPre == 500
  error('The Fixation period above should be defined for this paradigm, and it isnt');
end

epochTargParams.stimParamsFilename = stimParamsFilename;
epochTargParams.nonParametric = 0;                      % switch to run non Parametric tests.
 
epochTargParams.times = [preFix; Fix; stimEarly; stimLate; reward];
epochTargParams.labels = {'preFix', 'Fix', 'stimEarly', 'stimLate', 'reward'};    

epochTargParams.naturalSocial.targNames = {'socVNonSoc', 'categories'};
epochTargParams.naturalSocial.targ = {{'agents', 'socialInteraction'}, {'chasing', 'fighting', 'grooming', 'mounting', 'idle', 'goalDirected', 'objects', 'scene'}};
epochTargParams.naturalSocial.targetEpochs = [0 1 1 1 1; 0 1 1 1 1];           % Which of the labeled time bins to do the comparison for, per group, defined in analysisGroups.stimulusLabelGroups.groups, where first element is target.
epochTargParams.naturalSocial.oneVsAll = [1 0]; % Switch to do a target v baseline (targ = whichever group gets id'd as 1) or an anova across all groups (should use for categories (0).

epochTargParams.headTurnCon.targNames = {'socVNonSoc', 'categories'};
epochTargParams.headTurnCon.targ = {{'agents', 'socialInteraction'}, {'chasing', 'fighting', 'grooming', 'mounting', 'idle', 'goalDirected', 'objects', 'scene'}};
epochTargParams.headTurnCon.targetEpochs = [0 1 1 1 1; 0 1 1 1 1];           % Which of the labeled time bins to do the comparison for, per group, defined in analysisGroups.stimulusLabelGroups.groups, where first element is target.
epochTargParams.headTurnCon.oneVsAll = [1 0];

epochTargParams.headTurnIso.targNames = {'model', 'headTurn', 'headvArms', 'turnToward',};
epochTargParams.headTurnIso.targ = {{'fullModel', 'smoothModel', 'dotModel'}, {'headTurn', 'noTurn'}, {'headIso', 'bioMotion'}, {'turnAway', 'turnToward'}};
epochTargParams.headTurnIso.targetEpochs = [[0 1 1 1 1]; [0 1 1 1 1]; [0 1 1 1 1]; [0 1 1 1 1]];           % Which of the labeled time bins to do the comparison for, per group, defined in analysisGroups.stimulusLabelGroups.groups, where first element is target.
epochTargParams.headTurnIso.oneVsAll = [1 1 1 1];

% model using ANOVAs to detect encoding of variable features.
epochSWparams.stimParamsFilename = stimParamsFilename;
epochSWparams.binSize = 150;
epochSWparams.binStep = 25;
% epochSWparams.startTime = -psthParams.psthPre;
epochSWparams.startTime = -psthParams.psthPre;

% epochSWparams.naturalSocial.testLabel = {'broadCat', 'categories', 'broadCatEyes', 'categoriesEyes'};
% epochSWparams.naturalSocial.testCategoryLabels = {{'socialInteraction', 'idle', 'goalDirected', 'objects', 'scene'}, {'chasing', 'fighting', 'grooming', 'mounting', 'idle', 'goalDirected', 'objects', 'scene'}, ...
%                                                   {'socialInteraction', 'idle', 'goalDirected', 'objects', 'scene'}, {'chasing', 'fighting', 'grooming', 'mounting', 'idle', 'goalDirected', 'objects', 'scene'}};
% epochSWparams.naturalSocial.includeEyes = [false false true true];

epochSWparams.naturalSocial.testLabel = {'categories', 'categoriesEyes'};
epochSWparams.naturalSocial.testCategoryLabels = {{'chasing', 'fighting', 'grooming', 'mounting', 'idle', 'goalDirected', 'objects', 'scene'}, ...
                                                  {'chasing', 'fighting', 'grooming', 'mounting', 'idle', 'goalDirected', 'objects', 'scene'}};
epochSWparams.naturalSocial.includeEyes = [false true];

epochSWparams.headTurnCon = epochSWparams.naturalSocial;

epochSWparams.headTurnIso.testLabel = {'Model', 'headTurn', 'headvArms'};
epochSWparams.headTurnIso.testCategoryLabels = {{'fullModel', 'smoothModel', 'dotModel'}, {'headTurn', 'noTurn'}, {'headIso', 'bioMotion'}};
epochSWparams.headTurnIso.includeEyes = [false false false false];

epochSWparams.targNames = {'socInt', 'headTurn', 'fullModel', 'turnToward'};       % The names which end up in the table row names.
epochSWparams.targ = {'socialInteraction', 'headTurn', 'fullModel', 'turnToward'}; % The group members to be targeted for comparison against the rest.

% neuroGLMParams.neuroGLMPath = neuroGLMPath;
% neuroGLMParams.psthPre = psthParams.psthPre; % if e.g. +200, then start psth 200ms before trial onset; 
% neuroGLMParams.psthImDur = psthParams.psthImDur;  % only need to set this for variable length stim runs; else, comes from log file
% neuroGLMParams.psthPost = psthParams.psthPost;
% neuroGLMParams.visualizeDesignMat = 1;        % Visualizes a random set of 5 trials from the design matrix.

%% Plotting Params
assert(length(frEpochsCell) == length(epochLabels), 'Epoch time bins and epochLabel lengths must match')
%%%% note: all analysisGroups cell arrays are nx1, NOT 1xn

% familiarFace analysisGroups
monkeyIDList = {'Alan', 'Red', 'Otis', 'Pancho', 'Calvin', 'Diego', 'Barney', 'Hobbes'};
monkeyActList = {'IdleC', 'Movement', 'goalDirected', 'IdleS', 'FearGrimace', 'headTurn'};

monkey_act_List = [];
for ii = 1:length(monkeyIDList)
  for jj = 1:length(monkeyActList)
    monkey_act_List = [monkey_act_List; {[monkeyIDList{ii} '_' monkeyActList{jj}]}];
  end
end

analysisGroups.analysisGroupPSTH.familiarFace.Identity = monkeyIDList;
analysisGroups.analysisGroupPSTH.familiarFace.Action = monkeyActList;
analysisGroups.analysisGroupPSTH.familiarFace.stimuli = monkey_act_List;
analysisGroups.analysisGroupPSTH.familiarFace.Familair = {'familiar', 'unFamiliar'};

% headTurnCon analysisGroups
analysisGroups.analysisGroupPSTH.headTurnCon.stimuli = {'chasing1_Turn', 'chasing1_noTurn', 'chasing2_Turn', 'chasing2_noTurn', 'fighting1_noTurn', 'fighting2_Turn', 'fighting2_noTurn', 'grooming1_Turn', 'grooming1_noTurn', 'grooming2_Turn', 'grooming2_noTurn',...
    'idle1_Turn', 'idle1_noTurn', 'idle2_Turn', 'idle2_noTurn', 'mating1_Turn', 'mating1_noTurn', 'mating2_Turn', 'mating2_noTurn', 'objects1_noTurn', 'goalDirected1_Turn', 'goalDirected1_noTurn', 'goalDirected2_Turn'...                      }
    'goalDirected2_noTurn', 'objects2_noTurn', 'scene1_noTurn', 'scene2_noTurn'};
analysisGroups.analysisGroupPSTH.headTurnCon.socContext = {'socialInteraction', 'nonSocialInteraction'};
analysisGroups.analysisGroupPSTH.headTurnCon.turnNoTurn = {'noTurn', 'headTurn'};

% headTurnIso analysisGroups
analysisGroups.analysisGroupPSTH.headTurnIso.turnNoTurn = {'leftFull', 'noTurn', 'rightFull', 'headTurn'};
analysisGroups.analysisGroupPSTH.headTurnIso.stimuli = {'idle_core_frontCam_Dots_Early', 'idle_core_frontCam_Dots_Late', 'idle_core_frontCam_LowRes_Early', 'idle_core_frontCam_LowRes_Late', 'idle_core_frontCam_Normal_Early', 'idle_core_frontCam_Normal_Late',...          }
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
analysisGroups.analysisGroupPSTH.headTurnIso.turnSubj = {'turnToward', 'turnAway'};
analysisGroups.analysisGroupPSTH.headTurnIso.mesh = {'fullModel_headTurn', 'fullModel_noTurn', 'smoothModel_headTurn', 'smoothModel_noTurn', 'dotModel_headTurn', 'dotModel_noTurn'};
  
%Defined for Groups of 2, A-B/A+B type index.
analysisGroups.selectivityIndex.groups = {{'socialInteraction';'nonInteraction'},{'socialInteraction';'agents'}};

%Barplots showing average activity across all members of a catagory
analysisGroups.stimPrefBarPlot.groups = {{{'noTurn', 'headTurn'}, {'chasing', 'grooming', 'mounting', 'fighting', 'idle', 'goalDirected', 'objects', 'scene'}}};
analysisGroups.stimPrefBarPlot.colors  = {{{[0.55 0.13 0.16];[0.93 .2 0.15];[.98 0.65 0.13];[0 0.55 0.25];[0.15, 0.20, 0.5];[0.15, 0.20, 0.5]; [.5 0 .5]; [0 0.5 0.5]}, {[0.55 0.13 0.16];[0.93 .2 0.15];[.98 0.65 0.13];[0 0.55 0.25];[0.15, 0.20, 0.5];[0.15, 0.20, 0.5]; [.5 0 .5]; [0 0.5 0.5]}}};
analysisGroups.stimPrefBarPlot.names = {'Barplots per label'};
analysisGroups.stimPrefBarPlot.groupDepth = 2; %2 subplots, each figure is defined by a cell array in the first item (groups).

%
analysisGroups.stimulusLabelGroups.groups = {{'socialInteraction'; 'goalDirected'; 'idle'; 'objects'; 'scene'; 'scramble'},{'headTurn', 'noTurn'}, {'fullModel', 'dotModel', 'smoothModel'}, {'turnToward', 'turnAway'}};
analysisGroups.stimulusLabelGroups.names = {'Soc V Non Soc', 'Head Turn', 'Model', 'Turn Dir'};
analysisGroups.stimulusLabelGroups.colors = {{[0.55 0.13 0.16];[0.93 .2 0.15];[.98 0.65 0.13];[0 0.55 0.25];[0.15, 0.20, 0.5];[0.15, 0.20, 0.5]; [.5 0 .5]}, {[0.55 0.13 0.16];[0.93 .2 0.15]}, {[0.55 0.13 0.16];[0.93 .2 0.15];[.98 0.65 0.13]}, {[0.55 0.13 0.16];[0.93 .2 0.15]}};

%Essentially LFP selectivity/strength/quality
analysisGroups.evokedPotentials.groups = {{'socialInteraction';'nonInteraction';'objects'};{'socialInteraction';'nonInteraction'}};
analysisGroups.evokedPotentials.names = {'socVobj';'Social v Non-Soc'};
analysisGroups.evokedPotentials.colors = {{[0.1 0.1 1];[0.1 .7 0.1];[1 0.1 0.1]};{[0.1 0.1 1];[0.1 .7 0.1];[1 0.1 0.1]}};

%Looks like evokedpotentials, but pulls from a different analog channels
%(like pupul or eye).
analysisGroups.analogInPotentials.groups = {};
analysisGroups.analogInPotentials.channels = {[1; 2]};
analysisGroups.analogInPotentials.names = {'eyePositions,fobPlus'};
analysisGroups.analogInPotentials.units = {'degrees visual angle'};
analysisGroups.analogInPotentials.colors = {{[0.1 0.1 1];[1 0.1 0.1];'k';[0.1 .7 0.1];'m';[1 0.1 0.1];'k'}};

%Same thing, but derivatives.
analysisGroups.analogInDerivatives.groups = {};
analysisGroups.analogInDerivatives.channels = {[1; 2]};
analysisGroups.analogInDerivatives.names = {'eyeVelocity,fobPlus'};
analysisGroups.analogInDerivatives.units = {'degrees visual angle/sec'};
analysisGroups.analogInDerivatives.colors = {{[0.1 0.1 1];[0.1 .7 0.1];[1 0.1 0.1];[0.1 .7 0.1];'m';[1 0.1 0.1];'k'}};

%Makes subplots w/ PSTH on top and evoked potential on the bottom
analysisGroups.colorPsthEvoked.groups = {{'socialInteraction';'nonInteraction';'objects';'landscapes'}};
analysisGroups.colorPsthEvoked.names = {'socVobj'};
analysisGroups.colorPsthEvoked.colors = {{[0.1 0.1 1];[0 .6 0];[1 0.1 0.1]}};

%same as above, but shows error bars, harder to see catagory selectivity
%though
analysisGroups.linePsthEvoked.groups = {{'socialInteraction';'nonInteraction';'objects'}};
analysisGroups.linePsthEvoked.names = {'socVobj'};
analysisGroups.linePsthEvoked.colors = {{[0.1 0.1 1];[0 .6 0];[1 0.1 0.1]}};

%Same as above, but everything ontop of eachother in 1 panel.
analysisGroups.evokedPsthOnePane.groups = {};
analysisGroups.evokedPsthOnePane.names = {'faceVnon'};

%Creates tuning curves for units, you should have some meaningful numeric
%descriptor.
analysisGroups.tuningCurves.groups = {}; %can be images or categories
analysisGroups.tuningCurves.paramValues = {[-90 -45 0 45 90], [-90 -45 0 45 90]};
analysisGroups.tuningCurves.paramLabels = {'viewing angle (degrees)','viewing angle (degrees)'};
analysisGroups.tuningCurves.names = {'Human face view','Monkey face view'};

%Power spectra of SPIKE TIMES (i.e. 1's at spike times, flat elsewhere).
%Used for both Spikes and LFPs.
analysisGroups.spectraByCategory.groups = {{'interaction';'scrambles';}};  %todo: add spectra diff?
analysisGroups.spectraByCategory.names = {'Interaction V Scrambles'};
analysisGroups.spectraByCategory.colors = {{[1 0.1 0.1];[0.1 0.1 1]}};

%Calculates the same spectra, but w/ sliding windows. sometimes the Power
%spectrum changes overtime.
analysisGroups.tfSpectraByCategory.groups = {{'socialInteraction';'nonInteraction';}};%{'object'};{'body'}      %todo: add tf spectra diff?
analysisGroups.tfSpectraByCategory.names = {{'socialInt';'Int'}};%'nonface';'object';'body'

%Evoked potential plot, shows individual traces for a bunch of trials.
analysisGroups.lfpSingleTrialsByCategory.groups = {{'socialInteraction';'nonInteraction';}};
analysisGroups.lfpSingleTrialsByCategory.names = {'SocialVNonSocial'};

%Coherence between LFP time series and spike time series w/i single
%channel. Does this on a trial by trial basis, and then averages across all
%members of each group. 
analysisGroups.coherenceByCategory.groups = {{'interaction';'scrambles';}};
analysisGroups.coherenceByCategory.colors = {{[1 0.1 0.1];[0.1 0.1 1]}}; 
analysisGroups.coherenceByCategory.names = {'Interaction V Scrambles'}; %'fob';'slimCats'

%Calculates the same as above but in sliding windows.
analysisGroups.tfCouplingByCategory.groups = {{'socialInteraction'};{'nonInteraction';};{'objects'};{'goalDirected'}}; %#ok
%%%%%

% Do a quick check on analysisGroups just defined to conform to proper
analysisGroupFields = fields(analysisGroups);
for analysisGroup_i = 1:length(analysisGroupFields)
  field = analysisGroupFields{analysisGroup_i};
  
  % If the field doesn't conform to this checker's syntax, or is more
  % complicated, then do not process it in this way.
  hasGroups = isfield(analysisGroups.(field), 'groups');
  if ~hasGroups || (hasGroups && isempty(analysisGroups.(field).groups))
    continue
  end
  if isfield(analysisGroups.(field),'groupDepth') && analysisGroups.(field).groupDepth > 1
    continue
  end
  
  % Check that all the elements of the analysisGroup have the same length
  assert(length(unique(structfun(@(x) length(x), analysisGroups.(field)))) == 1, "analysisGroups.%s has uneven group members - all subfields should have the same length", field)
  
end

if calcSwitch.useJacknife
  chronuxParams.err(1) = 2; %#ok
end

%% set paths and directories, EDIT RARELY %%
analogInFilename = sprintf('%s/%s/%s%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);
lfpFilename = sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum);        
spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
taskFilename = sprintf('%s/%s/%s%s.mat',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance 
[photodiodeFilename, lineNoiseTriggerFilename] = deal(lfpFilename);                %#ok
outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
analysisParamFilename = strcat(outDir,analysisParamFilenameStem);
preprocessedDataFilename = strcat(outDir,preprocessedDataFilenameStem);                     %#ok

% In case difference logfile is being used.
if ~logical(exist(taskFilename,'file'))
  [A, B, C] = fileparts(taskFilename);
  switch C
    case '.mat'
      taskFilename = [A '/' B '.bhv2'];
    case '.bhv2'
      taskFilename = [A '/' B '.mat'];
  end
end

% Check if the key file types exist.
assert(logical(exist(analogInFilename,'file')),'The analog input file you requested does not exist.');
assert(logical(exist(taskFilename,'file')),'The log file you requested does not exist.');

% Auto channel detect - if set, find the parsed file (or parse it), and
% reassign spike channels on this.
if autoChannelDetect
  %If the file is parsed, retrieve channels present
  parsedFolderName = sprintf('%s/%s/%s%s_parsed',ephysVolume,dateSubject,dateSubject,runNum);
  if exist(parsedFolderName,'dir')
    [ephysParams.spikeChannels, ephysParams.lfpChannels, ephysParams.channelNames] = autoDetectChannels(parsedFolderName);
  else
    error("folder doesn't exist, change to non-autoChannelDetect");
  end
end

if ~exist(outDir,'dir')
  mkdir(outDir);
end
if isempty(varargin) 
  save(analysisParamFilename);
elseif strcmp(varargin,'saveNoPreprocParams')
  save(analysisParamFilename,'calcSwitch','analysisGroups','plotSwitch','-append');
end
end

function swappedString = slashSwap(pathString)
%Swaps direction of slashes to match Unix/Phyzzy, from Windows Path.
  stringParts = split(pathString, '\');
  swappedString = char(join(stringParts, '/'));
end
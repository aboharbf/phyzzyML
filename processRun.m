function [runAnalysisInputs, analysisOutFilename] = processRun( varargin )
%processRun loads, preprocesses, and aligns task events, lfps, muas, units, and analog signals (eg eye, accelerometer,photodiodes)
%   - handles single-channel and multi-channel sessions
%   - relies only on raw task log and physiology fiels: currently visiko (.log) and blackrock (.ns5,.ns2)
%   - aligns task events to ephys system clock (preprocessLogFile.m)
%   - excludes trials with broken fixation, fix spot flash, and (optionally)
%     juice delivery and high head acceleration(via exludeStimuli)
%   - decimates and (optionally) filters raw lfp data (default, to 1 kHz)
%   - calibrates eye signals (via preprocessEyeSignals)
%   - calibrates accelerometer signals (via processAccelSignals)
%   - calculates spike isolation measures (via preprocessSpikes)
%   - shows task event summary (via excludeTrials)
%   - optionally calls an analysis routine; default is runAnalyses.m
%   Inputs:
%   - varargin consists of name-value argument pairs:
%       - 'paramBuilder', paramBuilderFunctionName
%       - 'paramFile', paramFilename (must be .mat)
%       - 'analyzer', analyzerFunctionName
%       - 'preprocessed', preProcessedDataFilename (must be .mat), or '-d'
%             for default, in which case preprocessedDataFilename will be
%             extracted from the dateSubj+runNum in buildAnalysisParamFile,
%             paramBuilder, or paramFile. Unless paramFile is supplied,
%             will overwrite the calcSwitch, plotSwitch, and analysisGroups
%             fields of the saved analysisParamFile, based on their current
%             values in the param file builder. 
%       - 'preprocessedCC', preProcessedDataFilename (must be .mat), or '-d'
%             identical behavior to preprocessed, except that no variables
%             in analysisParamFile will be updated or overwritten. This is
%             useful to reproduce an analysis, or to run identically specified analyses
%             with a different or modified analyzer.
%       - 'keepItemsNotPresented', logical, default 0. If 0, events that
%           appear in the stimulus param files, but that have no valid trials,
%           do not appear in output data structures. 
%       Note: functionNames are strings, and do not include the trailing '.m'; see 'feval' docs
%   Notes:
%   Depends:
%   - contents of 'dependencies' folder.
%   - R2016a (or later) if joint psth-evoked potential plots desired
%   - Signal Processing Toolbox (for dpss taper calculation, LFP filters)

addpath(genpath('dependencies'));
rmpath(genpath('dependencies/mvgc_v1.0')); %note: Not sure if this is appropriate replacement for genpath_exclude. previous line caused issues in parallel runs.
addpath('buildAnalysisParamFileLib');

%% load analysis parameters

%Default values
usePreprocessed = 0;
preprocessedCC = 0;
defaultPreprocessed = 0;
keepItemsNotPresented = 0;

%Parse inputs
assert(mod(nargin,2) == 0, 'processRun takes arguments as name-value pairs; odd number of arguments provided');
for argPair_i=1:nargin/2
  argName = varargin{1+2*(argPair_i-1)};
  argVal = varargin{2+2*(argPair_i-1)};
  switch argName
    case 'paramBuilder'
      paramBuilder = argVal;
    case 'analyzer'
      analyzer = argVal;
    case 'paramFile'
      analysisParamFilename = argVal;
    case 'keepItemsNotPresented'
      keepItemsNotPresented = argVal;
    case 'preprocessed'
      if strcmp(argVal,'-d')
        defaultPreprocessed = 1;
      else
        preprocessedDataFilename = argVal;
      end
      usePreprocessed = 1;
    case 'preprocessedCC'
      if strcmp(argVal,'-d')
        defaultPreprocessed = 1;
      else
        preprocessedDataFilename = argVal;
      end
      usePreprocessed = 1;
      preprocessedCC = 1;
    otherwise
      error('Invalid field %s provided to processRun',argName);
  end
end

%% ParamBuilder and preprocessed file handling

if usePreprocessed
  if defaultPreprocessed
    if preprocessedCC 
      if ~exist('analysisParamFilename','var')
        if exist('paramBuilder','var')
          analysisParamFilename = feval(paramBuilder,'noSave');
        else
          analysisParamFilename = buildAnalysisParamFile('noSave');
        end
      end
    else
      if ~exist('analysisParamFilename','var')
        if exist('paramBuilder','var')
          analysisParamFilename = feval(paramBuilder,'saveNoPreprocParams');
        else
          analysisParamFilename = buildAnalysisParamFile('saveNoPreprocParams');
        end
      end
    end
    load(analysisParamFilename,'preprocessedDataFilename');
  end
  load(preprocessedDataFilename);
else
  if ~exist('analysisParamFilename','var')
    if exist('paramBuilder','var')
      analysisParamFilename = feval(paramBuilder);
    else
      analysisParamFilename = buildAnalysisParamFile();
    end
  end
end

if ~usePreprocessed
  checkAnalysisParamFile(analysisParamFilename);
  load(analysisParamFilename);
  
  % set verbose level
  if strcmp(verbosity,'VERBOSE')
    Output.level(Output.DISP_VERBOSE);
  end
  if strcmp(verbosity,'DEBUG')
    Output.level(Output.DISP_DEBUG);
  end

  %% Preprocess inputs
  
  analogInData =                                          preprocessAnalogIn(analogInFilename, analogInParams); 
  [spikesByChannel, taskTriggers, channelUnitNames] =     preprocessSpikes(spikeFilename, ephysParams);
  lineNoiseTriggers =                                     preprocessStrobe(lineNoiseTriggerFilename, lineNoiseTriggerParams);
  lfpData =                                               preprocessLFP(lfpFilename, ephysParams, lineNoiseTriggers);
  diodeTriggers =                                         preprocessStrobe(photodiodeFilename, photodiodeParams);
  [taskData, stimTiming] =                                preprocessLogFile(taskFilename, taskTriggers, diodeTriggers, stimSyncParams); % load visual stimulus data and transform its timestamps to ephys clock reference
% [pre, post] =                                           spikeBackground(spikesByChannel, taskData, ephysParams.spikeChannels, params ) 
  [analogInData, taskData] =                              preprocessEyeSignals(analogInData,taskData,eyeCalParams);
  analogInData =                                          preprocessAccelSignals(analogInData, accelParams); 
  
  % determine the duration of ephys data collection, in ms
  if ~isempty(lfpData)
    ephysDuration = size(lfpData,2); 
  elseif ~isempty(analogInData)
    ephysDuration = size(analogInData,2); 
  elseif ~isempty(spikesByChannel)
    ephysDuration = 0;
    for channel_i = 1:length(spikesByChannel)
      ephysDuration = max(ephysDuration, max(spikesByChannel{channel_i}.times)); 
    end
    lastTaskTriggerTime = 1000*max(taskTriggers.TimeStampSec); 
    if lastTaskTriggerTime > ephysDuration
      ephysDuration = lastTaskTriggerTime;
      warning('Using final digital trigger time to define end of ephys data collection, because it followed the final spike, and because no LFP or analog in data available. May bias results');
    else
      warning('Using final spike time to define end of ephys data collection, because it followed the final digital trigger, and because no LFP or analog in data available. May bias results');
    end
  end
  excludeStimParams.ephysDuration = ephysDuration;
  
  taskDataAll = taskData;
  taskData = excludeTrials(taskData, excludeStimParams); %exclude trials for lost fixation etc. 
  
  %% find Stimulus length
  
  if stimTiming.shortest == stimTiming.longest
    
    spikeAlignParams.postAlign = psthParams.psthImDur + psthParams.psthPost + 3*psthParams.smoothingWidth;
    lfpAlignParams.msPostAlign = psthParams.psthImDur + psthParams.psthPost + tfParams.movingWin(1)/2;
    if length(lfpAlignParams.DCSUB_SAM) > 1 && any(lfpAlignParams.DCSUB_SAM(2,:))  % todo: handle case where lfp analysis window doesn't cover full ISI
      lfpAlignParams.DCSUB_SAM(2,:) = (psthParams.psthImDur+stimTiming.ISI)+lfpAlignParams.DCSUB_SAM(2,:);
    end
    
  else
    
    assert(psthParams.psthImDur < stimTiming.longest, 'psthImDur is longer than longest trial; nothing to analyze');
    assert(psthParams.psthImDur > 0, 'psthImDur = 0; nothing to analyze. Likely cause: unexpected variable stimulus length');
    fprintf('Variable stimulus length run. Excluding trials shorter than %d\n', psthParams.psthImDur);
    
  end

%% Sort trials by image and image category, align data appropriately.
  %loads variables paramArray, categoryLabels,eventLabels, and extract eventIDs
  tmp = load(stimParamsFilename); 
  eventCategories = tmp.paramArray;
  categoryList = tmp.categoryLabels;
  refStruct = tmp.refStruct;
  try
    eventLabels = tmp.eventLabels;
  catch
    %implements back-compatibility with previous variable names
    eventLabels = tmp.pictureLabels; 
  end
  %per the stimParamFile spec, the first entry of each array is the event ID
  eventIDs = cellfun(@(x) x(1), eventCategories);   

  % Identify offsets and onsets for each eventID in the stimParam file
  % during this run.
  [onsetsByEvent, offsetsByEvent, trialIDsByEvent, fixOnsetsByEvent] = deal(cell(size(eventIDs)));
  for event_i = 1:length(eventIDs)
    wholeRunEventInd = strcmp(taskData.taskEventIDs, eventIDs{event_i});
    onsetsByEvent{event_i} = taskData.taskEventStartTimes(wholeRunEventInd);
    offsetsByEvent{event_i} = taskData.taskEventEndTimes(wholeRunEventInd);
    trialIDsByEvent{event_i} = find(wholeRunEventInd);
    fixOnsetsByEvent{event_i} = taskData.taskEventFixDur(wholeRunEventInd);
  end
  eventsObserved = ~cellfun('isempty', onsetsByEvent);
  assert(any(eventsObserved), 'No Events observed. Confirm correct stimulus Param file is in use.');
  
  % Check that the events observed match the reported events from the log
  % preprocessor. Throw errors for unrepresented events.
  uniqueEventsInRun = unique(taskData.taskEventIDs);
  eventsByOnset = eventIDs(eventsObserved);
  mismatchEvents = setdiff(uniqueEventsInRun, eventsByOnset);
  if ~isempty(mismatchEvents)
    disp(mismatchEvents)
    error('Events missing from stimParamFile');
  end

  % todo: need defense against image with onset but no offset? 
  % todo: add similar defense for rf map locations here?
  if ~taskData.RFmap
    %disp('No presentations of the following images survived exclusion:');
    %disp(eventIDs(eventsObserved));
  end
  
  % Once eventIDs is finished, rearrange things in taskData which are
  % dependent on this order, just to be sure
  
%   assert(all(eventIDResort' == 1:length(eventIDResort)), 'taskEventIDs from taskData not in the same order as eventIDs - consider investigating further')
  
  eventsObserved = find(eventsObserved);
  eventIDsTmp = eventIDs(eventsObserved);
  [~, eventIDResort] = ismember(taskData.taskEventList, eventIDsTmp);
  eventsObserved = eventsObserved(eventIDResort);
  
  if ~keepItemsNotPresented
    refStruct = tmp.refStruct(eventsObserved);
    onsetsByEvent = onsetsByEvent(eventsObserved);
    offsetsByEvent = offsetsByEvent(eventsObserved);
    trialIDsByEvent = trialIDsByEvent(eventsObserved);
    fixOnsetsByEvent = fixOnsetsByEvent(eventsObserved);
    eventIDs = eventIDs(eventsObserved);
    eventLabels = eventLabels(eventsObserved);
    eventCategories = eventCategories(eventsObserved);
  end
  
  % Identify onsets and offsets by Catagory
  [onsetsByCategory, offsetsByCategory, trialIDsByCategory] = deal(cell(length(categoryList),1));
  for cat_i = 1:length(categoryList)
    [catOnsets, catOffsets, catTrialIDs] = deal([]);
    eventsWithCat = find(cellfun(@(x) any(strcmp(x, categoryList{cat_i})), eventCategories)); % Find events w/ current catagory label.
    
    for event_i = eventsWithCat
      catOnsets = vertcat(catOnsets, onsetsByEvent{event_i});
      catOffsets = vertcat(catOffsets, offsetsByEvent{event_i});
      catTrialIDs = vertcat(catTrialIDs, trialIDsByEvent{event_i});
    end
    
    onsetsByCategory{cat_i} = catOnsets;
    offsetsByCategory{cat_i} = catOffsets;
    trialIDsByCategory{cat_i} = catTrialIDs;
    
    Output.DEBUG('numel cat onsets');
    Output.DEBUG(numel(catOnsets));
  end
  catsObserved = ~cellfun('isempty', onsetsByCategory);
  
  if ~keepItemsNotPresented
    onsetsByCategory = onsetsByCategory(catsObserved);
    offsetsByCategory = offsetsByCategory(catsObserved);
    categoryList = categoryList(catsObserved);
    trialIDsByCategory = trialIDsByCategory(catsObserved);
  end
    
  if taskData.RFmap
    jumpsByImage = cell(size(eventIDs));
    for i = 1:length(eventIDs)
      jumpsByImage{i} = taskData.stimJumps(strcmp(taskData.taskEventIDs,eventIDs{i}),:);
    end
  else
    jumpsByImage = [];
  end
  
  % align spikes by trial, and sort by image and category
  spikeAlignParamsToCoverMovingWin.preAlign = lfpAlignParams.msPreAlign; % need to include spikes throughout the window chronux will use to calculate spectra
  spikeAlignParamsToCoverMovingWin.postAlign = lfpAlignParams.msPostAlign;
  spikeAlignParamsToCoverMovingWin.refOffset = 0;
  [spikesByEvent, psthEmptyByEvent] = alignSpikes(spikesByChannel, onsetsByEvent, ephysParams.spikeChannels, spikeAlignParamsToCoverMovingWin);
  [spikesByCategory, psthEmptyByCategory] = alignSpikes(spikesByChannel, onsetsByCategory, ephysParams.spikeChannels, spikeAlignParamsToCoverMovingWin);
  
  % Generate a second array of spikesAligned to the fixation cue rather than the stimulus - mostly for fixation analysis and baseline.
  spikeAlignParamsToCoverMovingWin.preAlign = lfpAlignParams.msPreAlign; % need to include spikes throughout the window chronux will use to calculate spectra
  spikeAlignParamsToCoverMovingWin.postAlign = lfpAlignParams.msPostAlign;
  spikeAlignParamsToCoverMovingWin.refOffset = 0;
  [spikesByEventFixAlign, ~] = alignSpikes(spikesByChannel, fixOnsetsByEvent, ephysParams.spikeChannels, spikeAlignParamsToCoverMovingWin);
  [spikesByCategoryFixAlign, ~] = alignSpikes(spikesByChannel, onsetsByCategory, ephysParams.spikeChannels, spikeAlignParamsToCoverMovingWin);

  % align spikes again, but this time reference time to lfp sample number (required for chronux TF, even spike-spike)
  spikeAlignParamsTF.preAlign = lfpAlignParams.msPreAlign;
  spikeAlignParamsTF.postAlign = lfpAlignParams.msPostAlign;
  spikeAlignParamsTF.refOffset = -lfpAlignParams.msPreAlign;
  spikesByEventForTF = alignSpikes( spikesByChannel, onsetsByEvent, ephysParams.spikeChannels, spikeAlignParamsTF );
  spikesByCategoryForTF = alignSpikes( spikesByChannel, onsetsByCategory, ephysParams.spikeChannels, spikeAlignParamsTF );
  
  % align LFP data  by trial, sort by image and category, and possibly remove DC and linear components
  lfpByEvent = alignLFP(lfpData, onsetsByEvent, ephysParams.lfpChannels, lfpAlignParams);
  lfpByCategory = alignLFP(lfpData, onsetsByCategory, ephysParams.lfpChannels, lfpAlignParams);
  analogInByEvent = alignAnalogIn(analogInData, onsetsByEvent, analogInParams.analogInChannels, lfpAlignParams);
  analogInByCategory = alignAnalogIn(analogInData, onsetsByCategory, analogInParams.analogInChannels, lfpAlignParams);
  
  %   for cat_i = 1:length(categoryList)
  %     Output.VERBOSE(categoryList{cat_i});
  %     Output.VERBOSE(size(lfpByCategory{cat_i}));
  %   end
  
  if savePreprocessed
    save(preprocessedDataFilename,'analysisParamFilename', 'spikesByChannel', 'lfpData', 'analogInData', 'taskData', 'taskDataAll', 'refStruct', 'categoryList', 'eventLabels', 'eventIDs', ...
      'jumpsByImage', 'spikesByEvent', 'psthEmptyByEvent', 'spikesByEventFixAlign', 'spikesByCategory', 'psthEmptyByCategory', 'spikesByCategoryFixAlign', 'spikesByEventForTF', 'spikesByCategoryForTF', ...
      'lfpByEvent', 'lfpByCategory', 'analogInByEvent', 'analogInByCategory', 'channelUnitNames', 'stimTiming', 'eventCategories', 'onsetsByEvent', 'offsetsByEvent', 'onsetsByCategory', 'offsetsByCategory', ...
      'trialIDsByEvent', 'trialIDsByCategory', 'excludeStimParams');
  end
end

%% Package Analysis inputs, and run Analysis
runAnalysisInputs.analysisParamFilename = analysisParamFilename;
runAnalysisInputs.spikesByChannel = spikesByChannel;
runAnalysisInputs.lfpData = lfpData;
runAnalysisInputs.analogInData = analogInData;
runAnalysisInputs.taskData = taskData;
runAnalysisInputs.taskDataAll = taskDataAll;
runAnalysisInputs.refStruct = refStruct;
runAnalysisInputs.categoryList = categoryList;
runAnalysisInputs.eventLabels = eventLabels;
runAnalysisInputs.eventIDs = eventIDs;
runAnalysisInputs.jumpsByImage = jumpsByImage;
runAnalysisInputs.spikesByEvent = spikesByEvent;
runAnalysisInputs.spikesByEventFixAlign = spikesByEventFixAlign;
runAnalysisInputs.psthEmptyByEvent = psthEmptyByEvent;
runAnalysisInputs.spikesByCategory = spikesByCategory;
runAnalysisInputs.spikesByCategoryFixAlign = spikesByCategoryFixAlign;
runAnalysisInputs.psthEmptyByCategory = psthEmptyByCategory;
runAnalysisInputs.spikesByEventForTF = spikesByEventForTF;
runAnalysisInputs.spikesByCategoryForTF = spikesByCategoryForTF;
runAnalysisInputs.lfpByEvent = lfpByEvent;
% runAnalysisInputs.lfpByCategory = lfpByCategory;
runAnalysisInputs.analogInByEvent = analogInByEvent;
runAnalysisInputs.analogInByCategory = analogInByCategory;
runAnalysisInputs.channelUnitNames = channelUnitNames;
runAnalysisInputs.stimTiming = stimTiming;
runAnalysisInputs.eventCategories = eventCategories;
runAnalysisInputs.onsetsByEvent = onsetsByEvent;
runAnalysisInputs.offsetsByEvent = offsetsByEvent;
runAnalysisInputs.trialIDsByEvent = trialIDsByEvent;
runAnalysisInputs.onsetsByCategory = onsetsByCategory;
runAnalysisInputs.offsetsByCategory = offsetsByCategory;
runAnalysisInputs.trialIDsByCategory = trialIDsByCategory;
runAnalysisInputs.excludeStimParams = excludeStimParams;

if ~exist('analyzer','var')
  [analysisOutFilename] = runAnalyses(runAnalysisInputs);
else
  feval(analyzer,runAnalysisInputs);
end
end


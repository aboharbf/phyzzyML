function [ taskData, stimTiming ] = preprocessLogFileMonkeyLogic(logfile, taskTriggers, diodeTriggers, params)
%% Reads in Black
%  - in presentation order, stimulus filenames, start/end times and jump RF mapping positions
%  - task event times: fixation in/out, fixspot flash start/end, juice delivery start/end 
%  - currently, uses digital IO packets for synchronization. In future
%    implementation, will optionally refine this sync with photodiode traces
%   Inputs:
%   - logFilename: an xml file generated by Visiko (or, in other
%     implementations, some other stimulation software)
%   - taskTriggers: digital IO data from blackrock NEV file
%   - params: currently must be struct of form params.usePhotodiode = 0
%   Outputs:
%   - taskData: struct with fields:
%       - taskEventIDs: nTrials x 1 cell array of alignment point identifiers, e.g. stimulus filenames. 
%         Used to retrieve information from stimulusParamFile.
%       - taskEventStartTimes: nTrials x 1 array of task event (e.g. stimulus) start times in ms
%       - taskEventEndTimes: nTrials x 1 array of task event (e.g. stimulus) end times in ms
%       _ stimParams: struct containing information about stimulus size etc. Not currently used elsewhere, so may be left empty
%       - RFmap: 1 for runs where stimulus position varies, else 0
%       - stimJumps: if RFmap, nTrialsx2 array or stim x and y positions, in degrees of visual angle, else may be empty
%           - note: position relative to center of jump grid, not absolute position
%       - gridPointsDegX: if RFmap, 1xN array of unique stimulus jump X locations, else may be empty 
%       - gridPointsDegY: if RFmap, 1xN array of unique stimulus jump Y locations, else may be empty
%       - fields likely required for excludeStimuli and runSummary plots (may be left empty if not needed):
%         - stimFramesLost: nTrials x 1, number of frames lost during trial
%         - fixationInTimes: nTrials x 1, times in ms when fixation epochs begin
%         - fixationOutTimes: nTrials x 1, times in ms when fixation epochs end
%         - juiceOnTimes: nTrials x 1, times in ms when juice delivery begins
%         - juiceOffTimes: nTrials x 1,times in ms when juice delivery ends
%         - fixSpotFlashStartTimes: nTrials x 1, times in ms when flash begins
%         - fixSpotFlashEndTimes: nTrials x 1, times in ms when flash ends
%   Dependencies:
%   - Statistics and Machine Learning Toolbox (for synchronizaton)
%   - xml2struct (from Matlab fileExchange)

%% Process MonkeyLogic File
%Parse the Log file
disp('Loading MonkeyLogic log file');
assert(logical(exist(logfile,'file')),'The logfile you requested does not exist.');
[data,MLConfig,TrialRecord] = mlread(logfile);
assert(length(unique(TrialRecord.ConditionsPlayed)) < ceil((length(TrialRecord.ConditionsPlayed))/2) , 'MonkeyLogic file reports each condition was not repeated at least 2 times')

%Find stimulus timing range - in ML, we have no range within valid trials.
stimTiming.shortest = 2800; %Making this 0 for now
stimTiming.longest = 2800; %setting this to 2800 for now, but maybe save as an editable variable.
stimTiming.ISI = MLConfig.InterTrialInterval;

%Construct the translation table from the logfile
%Find the Conditions that need to be connected to codes, and loop through
%the trials until each of those numbers has a stim attached.
allConditions = unique(TrialRecord.ConditionsPlayed); %Pull unique members of this list
translationTable = cell(length(allConditions),1); %Initialize translation table
trial_ind = 1;

while sum(find(cellfun('isempty', translationTable))) ~= 0
  trialHolder = data(trial_ind);
  stimName = trialHolder.TaskObject.Attribute(2).Name; %Pull the string containing the stimulus name.
  translationTable{trialHolder.Condition} = stimName(~isspace(stimName));
  trial_ind = trial_ind + 1;
end

logFileStimTable = strtrim((TrialRecord.TaskInfo.Stimuli)); %Clean up the stim table in the logfile.

%Check to see if what you collected by cycling through the stim matches
%this table.
assert(sum(strcmp(logFileStimTable,translationTable)) == length(translationTable), 'Stim table from MonkeyLogic doesnt match collected table');

%Pull Behavioral codes for other events
behavioralCodes = TrialRecord.TaskInfo.BehavioralCodes;

fixCueMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Fix Cue'));
trialStartMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Start trial'));
stimStartMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Stimuli On'));
stimEndMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Stimuli Off')); %Also means fixation off.
rewardMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Reward'));
frameSkipMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Frame skipped'));
trialEndMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'End trial'));
manualRewardMarker = behavioralCodes.CodeNumbers(strcmp(behavioralCodes.CodeNames,'Manual reward'));
fixFailMarker = 4;
stimFailMarker = 3;

%MonkeyLogic's output has an absolute trial start (time within entire run)
%and additional time stamps within trials, relative to that start. Collect
%all the start times, then create subsequent ones through the addition of
%the right amount. Step below removes failed trials.
errorArray = TrialRecord.TrialErrors';
correctTrialArray = (errorArray == 0);
mklTrialStarts = [data.AbsoluteTrialStartTime];
tmpStruct = [data.BehavioralCodes];
rwdTimes = [data.RewardRecord];

%initialize everything
[taskEventStartTimesLog, taskEventEndTimesLog, juiceOnTimesLog, juiceOffTimesLog, stimFramesLost] = deal(zeros(length(correctTrialArray), 1));

%Collect trialEventIDs.
taskEventIDsLog = TrialRecord.ConditionsPlayed';
  
%Note - the below setup pulls trials which fail during the stimuli, but not
%during the fixation period. for those trials, both start and end are NaN,
%for stim failures, the numbers are meaningful.
for ii = 1:length(mklTrialStarts)
  if ~isempty(tmpStruct(ii).CodeTimes(tmpStruct(ii).CodeNumbers == stimStartMarker))
    taskEventStartTimesLog(ii) = mklTrialStarts(ii) + tmpStruct(ii).CodeTimes(tmpStruct(ii).CodeNumbers == stimStartMarker);
    taskEventEndTimesLog(ii) = mklTrialStarts(ii) + tmpStruct(ii).CodeTimes(tmpStruct(ii).CodeNumbers == stimEndMarker);
  else
    taskEventStartTimesLog(ii) = nan;
    taskEventEndTimesLog(ii) = nan;
  end
  if ~isempty(rwdTimes(ii).StartTimes)
    juiceOnTimesLog(ii) = mklTrialStarts(ii) + rwdTimes(ii).StartTimes(end);
    juiceOffTimesLog(ii) = mklTrialStarts(ii) + rwdTimes(ii).EndTimes(end);
  else
    juiceOnTimesLog(ii) = nan;
    juiceOffTimesLog(ii) = nan;
  end
end

% Process the Eye data into a structure which can be used later to line up
% with Blackrock signals.
tmpEye = [data.AnalogData];
tmpEye = rmfield(tmpEye, {'SampleInterval', 'EyeExtra','Joystick','Mouse','PhotoDiode','General','Button'});

%Package the information necessary for proper eye/stimuli relationship.
tmp = split(MLConfig.Resolution,' ');
screenStats.screenSize = {str2double(tmp{1}), str2double(tmp{3})};
screenStats.PixelsPerDegree = abs(MLConfig.PixelsPerDegree);
screenStats.DiagonalSize = MLConfig.DiagonalSize;

%Translations - Function changes the name of stimuli from an Old convention
translationTable = updateTranslationTable(translationTable);

%Create Frame lost array
for ii = 1:length(data)
    stimFramesLost(ii) = sum(data(ii).BehavioralCodes.CodeNumbers == frameSkipMarker);
end

% Behavioral summary of performance during recording
taskData = struct();
taskData = behaviorsummaryPhyzzy(logfile, taskData);
[~, filename, ~] = fileparts(logfile);
savefig(sprintf('%sBehavSum_%s',params.outDir,filename));

%% Process the NEV derived data (Blackrock's Eventmarkers)
disp('parsing serisal IO packets');
packetTimes = double(taskTriggers.TimeStampSec)*1000; %convert seconds to milliseconds
packetData = double(taskTriggers.UnparsedData);

if isempty(packetData) % Means Blackrock/MKL Communication was not correctly connected.
  error('The Blackrock Digital inputs are empty. Digital inputs may have not been plugged in.');
else % Means Blackrock/MKL Communication was intact, correct
  %Code to get rid of any markers prior to the first trial beginning, or after the final trial end.
  %a defense against double 9's
  trueStart = find(packetData == trialStartMarker);
  if length(trueStart) > 1
    trueStartTrials = (diff(trueStart) > 1);
    trueStart = trueStart(trueStartTrials);
  end
  trueStart = trueStart(1);
  trueEnd = find(packetData == trialEndMarker,1,'last');
  
  packetData = packetData(trueStart:trueEnd);
  packetTimes = packetTimes(trueStart:trueEnd);
  
  %to later figure out how to shape data, we need to know what we saw and
  %got rid of.
  
  %Below are things which should happen every trial  
  trialStartInds = find(packetData == trialStartMarker); 
  trialEndInds = find(packetData == trialEndMarker);
  stimStartInds = packetData == stimStartMarker;
  condNumbers = packetData(packetData > 100)-100;
  fixStartInds = find(packetData == fixCueMarker);
  stimFixEndInds = find(packetData == stimEndMarker);
  
  assert(length(trialStartInds) == length(condNumbers), 'Every trial doesnt have a condition number - this may be due to changes in when the marker is sent.')
  
  %Now, use the correct trial array from MKL to pick out only the correct
  %trials from the whole array. Assign them them to the correct vectors. 
  [taskEventStartTimesBlk, taskEventEndTimesBlk, juiceOnTimesBlk, juiceOffTimesBlk, trialStartTimesBlk] = deal(nan(size(trialStartInds)));
  
  %Construct Error array from Blackrock files only, useful in the case that
  %incomplete trials lead to unpaired start and stop triggers.
  errorArray = zeros(sum(packetData == trialStartMarker), 1);
  for ii = 1:length(trialStartInds) %for every trial
    %Go through packet data, starting at that stim, and see if you hit a
    %40, 3, or 4 first.
    found = 0;
    stepsAhead = 1;
    while ~found
      switch packetData(trialStartInds(ii)+stepsAhead)
        case fixFailMarker
          errorArray(ii) = fixFailMarker;
          found = 1;
        case stimFailMarker
          errorArray(ii) = stimFailMarker;
          found = 1;
        case rewardMarker
          errorArray(ii) = 0;
          found = 1;
        otherwise
          stepsAhead = stepsAhead + 1;
      end
    end
  end
  
  %packetData is directly referenced for events which only happen
  %sometimes.
  taskEventIDsBlk = condNumbers;
  trialStartTimesAll = packetTimes(packetData == trialStartMarker);
  trialStartTimesBlk(errorArray ~= 4) = trialStartTimesAll(errorArray ~= 4);
  taskEventStartTimesBlk(errorArray ~= 4) = packetTimes(packetData == stimStartMarker); %Stores every trial where stim started (no fail during fix);
  taskEventEndTimesBlk(errorArray ~= 4) = packetTimes(stimFixEndInds(errorArray ~= 4));
  juiceOnTimesBlk(errorArray == 0) = packetTimes(packetData == rewardMarker);
  juiceOffTimesBlk(errorArray == 0) = packetTimes(trialEndInds(errorArray == 0));
  
  %use the condition number to fill out the taskEventsIDs using the
  %translation table.
  taskEventIDs = translationTable(condNumbers);
end

%% Run some checks comparing data collected from Blackrock, MonkeyLogic.
if (length(taskEventIDsLog) ~= length(taskEventIDsBlk))
  warning('Blackrock and MonkeyLogic report different numbers of correct trials. Attempting to find best alignment')
  %Odd situation, likely due to hitting runAnalyses prior to the MKL being
  %stopped. Seems like MKL doesn't have the complete record, Blackrock
  %does. In these cases, use the shorter of the two.
  lag = finddelay(taskEventIDsBlk, taskEventIDsLog);
  %Switch below written for a single case in July, unlikely to happen
  %again, needs more extensive testing.
  %If find delay has an answer - this might invalidate the length checking method, but it seems to not always work.
  %This means shifting things forward improves the pairwise match - means
  %something extra at the start.
  if (lag > 0)
    taskEventIDsLog = taskEventIDsLog((lag+1):end);
    taskEventStartTimesLog = taskEventStartTimesLog((lag+1):end);
    taskEventEndTimesLog = taskEventEndTimesLog((lag+1):end);
    juiceOnTimesLog = juiceOnTimesLog((lag+1):end);
    juiceOffTimesLog = juiceOffTimesLog((lag+1):end);
    tmpEye = tmpEye((lag+1):end);
  elseif (lag < 0)
    %Means shifting back improves match,
    lag = abs(lag);
    taskEventIDsBlk = taskEventIDsBlk((lag+1):end);
    taskEventStartTimesBlk = taskEventStartTimesBlk((lag+1):end);
    taskEventEndTimesBlk = taskEventEndTimesBlk((lag+1):end);
    juiceOnTimesBlk = juiceOnTimesBlk((lag+1):end);
    juiceOffTimesBlk = juiceOffTimesBlk((lag+1):end);
    taskEventIDs = taskEventIDs((lag+1):end);
  end
  %Recalculate lengths
  LogLen = length(taskEventIDsLog);
  BlkLen = length(taskEventIDsBlk);
  if LogLen < BlkLen %monkeyLogic stopped storing things early.
    %Chop the Blackrock events down to the size of what Mkl has stored
    taskEventIDsBlk = taskEventIDsBlk(1:LogLen);
    if sum(taskEventIDsBlk == taskEventIDsLog) == length(taskEventIDsLog) %If They are now a match
      taskEventStartTimesBlk = taskEventStartTimesBlk(1:LogLen);
      taskEventEndTimesBlk = taskEventEndTimesBlk(1:LogLen);
      stimFramesLost = stimFramesLost(1:LogLen);
      errorArray = errorArray(1:LogLen); %Unique here because it is a stored variable from MKL.
      juiceOnTimesBlk = juiceOnTimesBlk(1:LogLen);
      juiceOffTimesBlk = juiceOffTimesBlk(1:LogLen);
      taskEventIDs = taskEventIDs(1:LogLen);
    else
      error("Simple alignment failed - Blackrock and MonkeyLogic logs still don't match")
    end
  else %Blackrock stopped early.
    %Blackrock shut off before monkeyLogic was done.
    taskEventIDsLog = taskEventIDsLog(1:BlkLen);
    if sum(taskEventIDsBlk == taskEventIDsLog) == length(taskEventIDsLog) %If They are now a match
      %do the same chopping to monkeyLogic related structures
      %Collect trialEventIDs.
      taskEventStartTimesLog = taskEventStartTimesLog(1:BlkLen);
      taskEventEndTimesLog = taskEventEndTimesLog(1:BlkLen);
      stimFramesLost = stimFramesLost(1:BlkLen);
      juiceOnTimesLog = juiceOnTimesLog(1:BlkLen);
      juiceOffTimesLog = juiceOffTimesLog(1:BlkLen);
      tmpEye = tmpEye(1:BlkLen);
    else
      error("Simple alignment failed - Blackrock and MonkeyLogic logs still don't match")
    end
  end
end
assert(sum(taskEventIDsLog == taskEventIDsBlk) == length(taskEventIDsLog), 'Blackrock and MonkeyLogic do not have matching Event numbers.')

%% The strobe is the most reliable marker that the stimuli has begun - use this for event start times.
%We now have the start of trials, Blackrock.Eventmarkers. We want to move that to Blackrock.Strobe, which we 
%believe is more accurate. This will be done by finding the transitions for
%the strobe from High (white) to low (black) closest to the time stamp.

%What is the discrepency?
offsets = zeros(sum(packetData == rewardMarker),1);
taskEventStartTimesBlkPreStrobe = taskEventStartTimesBlk; %Saving these for Eye signal related processing.

%Using Photodiode to redefine start times of trials and stimulus.
if params.usePhotodiode
  % Identify which transitions represent the trial starts.
  if length(diodeTriggers.low) < length(diodeTriggers.mid)
    startTimeTrans = diodeTriggers.highToMid;
    endTimeTrans = diodeTriggers.midToHigh;
  else
    startTimeTrans = diodeTriggers.highToLow;
    endTimeTrans = diodeTriggers.lowToHigh;
  end
  
  % for every time, grab the closest and replace the recorded time w/ the
  % closest transition of the correct kind.
  for ii = 1:length(taskEventStartTimesBlk)
    
    % Trial Start times
%     if ~isnan(trialEventStartTimesBlk(ii))
%       [offsets(ii), ind_trueStart(ii)] = min(abs(endTimeTrans-trialEventStartTimesBlk(ii)));
%       taskEventStartTimesBlk(ii) = endTimeTrans(ind_trueStart(ii));
%     end
    % Seems Trial start times don't line up w/ transitions as well, mean 20 ms discrepencies.
  
    % Stim Start times
    if ~isnan(taskEventStartTimesBlk(ii))
      [offsets(ii), ind_trueStart] = min(abs(startTimeTrans-taskEventStartTimesBlk(ii)));
      taskEventStartTimesBlk(ii) = startTimeTrans(ind_trueStart);
    end
    
    % Stim End Times
    if ~isnan(taskEventStartTimesBlk(ii))
      [offsets(ii), ind_trueEnd] = min(abs(endTimeTrans-taskEventEndTimesBlk(ii)));
      taskEventEndTimesBlk(ii) = endTimeTrans(ind_trueEnd);
    end
    
  end
end

%% In cases where the stimuli are .avi's, we should check for the appropriate frameMotion data (representing)
if any(strfind(translationTable{1},'.avi'))
  %assumes filename below sitting in directory with stimuli
  frameMotionFile = dir([params.stimDir '/**/frameMotion_complete.mat']);
  load([frameMotionFile(1).folder filesep frameMotionFile(1).name],'frameMotionData');
  %go through the translation table, comparing to frameMotionData.stimVid
  frameMotionNames = [{frameMotionData(:).stimVid}'];
  tmpFrameMotionData = struct('stimVid',[],'objNames',[],'objShapes',[],'objRadii',[],'vidB_XShift',[],'objLoc',[],'frameCount',[],'timePerFrame',[],'fps',[],'width',[],'height',[]);
  for table_i = 1:length(translationTable)
    dataInd = find(strcmp(frameMotionNames, translationTable{table_i}));
    if ~isempty(dataInd)
      tmpFrameMotionData(table_i) = frameMotionData(dataInd);
    else
      tmpFrameMotionData(table_i) = frameMotionData(1);
      tmpFrameMotionData(table_i).objNames = [];
      tmpFrameMotionData(table_i).objShapes = [];
      tmpFrameMotionData(table_i).objRadii = [];
      tmpFrameMotionData(table_i).objLoc = [];
      tmpFrameMotionData(table_i).vidB_XShift = [];
      tmpFrameMotionData(table_i).stimVid = translationTable{table_i};
    end
  end
end

%% Now, calculate stimulation log to ephys clock conversion
%Make the Model by comparing Blackrock event start times to monkeylogic.
logVsBlkModel = fitlm(taskEventStartTimesLog, taskEventStartTimesBlkPreStrobe);
taskEventStartTimesFit = predict(logVsBlkModel, taskEventStartTimesLog); %Shift over all the events, using the calculated line.

% Gauge the quality of the fit, check if its too much adjustment.
eventTimeAdjustments = taskEventStartTimesFit-taskEventStartTimesBlkPreStrobe;
assert(max(eventTimeAdjustments) < 10, 'over 10 ms adjustment for strobe, check')
disp(strcat('Max magnitude fit residual, msec: ',num2str(max(abs(eventTimeAdjustments)))));

if params.plotFitShifts
  figure('Name','Sync adjustments from model fitting','NumberTitle','off')
  histogram(eventTimeAdjustments,40);
  title('Sync adjustments from model fitting');
  xlabel('offset (adjusted - original) (ms)');
  ylabel('count');
  disp(median(eventTimeAdjustments));
  [~,adjSortInds] = sort(abs(eventTimeAdjustments-median(eventTimeAdjustments)), 'descend');
  disp('worst alignments, log file times');
  disp(taskEventStartTimesLog(adjSortInds(1:min(5,length(adjSortInds)))));
  disp('worst alignments, adjusted times');
  disp(taskEventStartTimesBlk(adjSortInds(1:min(5,length(adjSortInds)))));
  disp('worst alignments, adjustment values (ms)');
  disp(eventTimeAdjustments(adjSortInds(1:min(5,length(adjSortInds)))));
end

%Assuming you're happy with this model, shift values only in the log to
%Blackrock time.
juiceOnTimesBlk = predict(logVsBlkModel, juiceOnTimesLog);
juiceOffTimesBlk = predict(logVsBlkModel, juiceOffTimesLog);
%Not doing this for event start times, since I trust the photodiode more
%than the eventmarker as a metric of true stim start time.

fprintf('average offset %s ms\n', num2str(mean(offsets)));
fprintf('range of offset %d ms - %d ms \n', [min(offsets), max(offsets)])
  
%% If en eventData file is present, search it, and generate composite subEvents.
% Code below was part of attempt to encorperate eventData into processRun
% code. Decision was made to move it all to runAnalyses section.
% eventDataFile = dir([params.stimDir '/**/eventData.mat']);
if 0%~isempty(eventDataFile)
  load(fullfile(eventDataFile(1).folder, eventDataFile(1).name), 'eventData');
  [eventDataRun, taskData.eventDataRun] = deal(eventData(translationTable, :));
  eventsInEventData = eventDataRun.Properties.VariableNames;
  % Cycle through events, generating startTimes.
  for event_i = 1:length(eventsInEventData)
    subEventRunData = eventDataRun(:, eventsInEventData{event_i});
    emptyIndTmp = rowfun(@(x) isempty(x{1}), subEventRunData);
    emptyInd = ~emptyIndTmp.Var1;
    subEventRunData = subEventRunData(emptyInd,:);
    stimSourceArray = subEventRunData.Properties.RowNames;
    for subEvent_i = 1:length(stimSourceArray)
      eventTable = subEventRunData{stimSourceArray(subEvent_i),:}{1};
      entryCount = size(eventTable,1);
      if subEvent_i == 1 && event_i == 1
        eventStimTable = [repmat(eventsInEventData(event_i), [entryCount, 1]), repmat(stimSourceArray(subEvent_i), [entryCount, 1]), eventTable(:,'startFrame'), eventTable(:,'endFrame')];
      else
        eventStimTable = [eventStimTable; [repmat(eventsInEventData(event_i), [entryCount, 1]), repmat(stimSourceArray(subEvent_i), [entryCount, 1]), eventTable(:,'startFrame'), eventTable(:,'endFrame')]];
      end
    end
  end
  
  % After Event compilation, begin generating the necessary output structures
  for event_i = 1:length(eventsInEventData)
    eventData = eventStimTable(strcmp(eventStimTable.Var1, eventsInEventData{event_i}),2:end);
    stimWithEvent = unique(eventData.Var2);
    for stim_i = 1:length(stimWithEvent)
      % Get the appropriate start and end frames, convert to 
      eventsInStim = table2array(eventData(strcmp(eventData.Var2, stimWithEvent(stim_i)), 2:3));
      frameMotionDataInd = strcmp({frameMotionData.stimVid}, stimWithEvent(stim_i));
      if ~any(frameMotionDataInd)
        stimVid = dir([params.stimDir '/**/' stimWithEvent{stim_i}]);
        vidObj = VideoReader(fullfile(stimVid(1).folder, stimVid(1).name));
        msPerFrame = (vidObj.Duration/vidObj.NumFrames) * 1000;
        clear vidObj
      else
        msPerFrame = frameMotionData(frameMotionDataInd).timePerFrame;
      end
      eventsInStimTimes = eventsInStim .* msPerFrame;
      
      % Find out when the stimulus happens
      stimInd = strcmp(taskEventIDs, stimWithEvent(stim_i));
      eventErrorArray = errorArray(stimInd);
      eventStartTimesTmp = taskEventStartTimesBlk(stimInd) + eventsInStimTimes(1);
      eventStartTimesTmp = eventStartTimesTmp(eventErrorArray == 0);
      eventEndTimesTmp = taskEventEndTimesBlk(stimInd) + eventsInStimTimes(2);
      eventEndTimesTmp = eventEndTimesTmp(eventErrorArray == 0);
      eventArrayTmp = repmat(eventsInEventData(event_i),[length(eventStartTimesTmp),1]);
      
      if stim_i == 1 && event_i == 1
        eventStartTimes = eventStartTimesTmp;
        eventEndTimes = eventEndTimesTmp;
        eventNames = eventArrayTmp;
      else
        eventStartTimes = [eventStartTimes; eventStartTimesTmp];
        eventEndTimes = [eventEndTimes; eventEndTimesTmp];
        eventNames = [eventNames; eventArrayTmp];
      end
    end
  end
  
  % Combine and sort output structures
  taskEventStartTimesBlk = [taskEventStartTimesBlk; eventStartTimes];
  taskEventEndTimesBlk = [taskEventEndTimesBlk; eventEndTimes];
  taskEventIDs = [taskEventIDs; eventNames];
  errorArray = [errorArray; zeros(length(eventNames),1)];
  stimFramesLost = [stimFramesLost; zeros(length(eventNames),1)];
  juiceOnTimesBlk = [juiceOnTimesBlk; nan(length(eventNames),1)];
  juiceOffTimesBlk = [juiceOffTimesBlk; nan(length(eventNames),1)];
end

%% Addtion in Aug 2020 - Realization of diversity of fixation times due to 'punishment' dynamic (failed trial = +50 msec on subsuquent fix time to initiate next trial)
  % Cycle through trials and collect pre-stimulus fixation duration
  tmp = [data.ObjectStatusRecord];
  taskEventFixDur = zeros(length(tmp),1);
  for ii = 1:length(tmp)
    taskEventFixDur(ii) = tmp(ii).SceneParam(1).AdapterArgs{3}{2,2};
  end
  
  % Establish a metric for the duration of the trial pre stimulus. Will be
  % useful for shifting things to a trial aligned, instead of stimulus
  % aligned, analysis.
  taskEventFixDur = taskEventStartTimesBlk - trialStartTimesBlk;
  assert(length(taskEventFixDur) == length(taskEventIDs), 'Problem w/ fixation times')
  
  % Generate a vector of reward times relative to stimulus onset.
  rewardTimePerTrial = juiceOnTimesBlk - taskEventStartTimesBlk;

%% Output
%Adding random numbers to these - they aren't relevant for my current task,
%nor are they directly recorded by MKL.
firstNonNaN = find(~isnan(taskEventStartTimesBlk),1);
fixSpotFlashStartTimesBlk = taskEventStartTimesBlk(firstNonNaN);
fixSpotFlashEndTimesBlk = taskEventEndTimesBlk(firstNonNaN);
fixationInTimesBlk = taskEventStartTimesBlk(firstNonNaN);
fixationOutTimesBlk = taskEventEndTimesBlk(firstNonNaN);

% finally, build the output structure:
% Note: **If something is missing in RunAnalyses** - likely due to taskData
% being overwritten following excludeTrial function.

taskData.taskDataSummary.TrialRecord = TrialRecord;
taskData.errorArray = errorArray;
taskData.taskEventIDs = taskEventIDs;
taskData.taskEventList = translationTable;
taskData.frameMotionData = tmpFrameMotionData';
taskData.stimFramesLost = stimFramesLost;

% Timing variables
taskData.taskEventStartTimes = taskEventStartTimesBlk;
taskData.taskEventEndTimes = taskEventEndTimesBlk;
taskData.taskEventFixDur = taskEventFixDur;
taskData.rewardTimePerTrial = rewardTimePerTrial;
taskData.trialStartTimesMkl = mklTrialStarts';
taskData.logVsBlkModel = logVsBlkModel;
taskData.fixationInTimes = fixationInTimesBlk;
taskData.fixationOutTimes = fixationOutTimesBlk;
taskData.juiceOnTimes = juiceOnTimesBlk;
taskData.juiceOffTimes = juiceOffTimesBlk;
taskData.fixSpotFlashStartTimes = fixSpotFlashStartTimesBlk;
taskData.fixSpotFlashEndTimes = fixSpotFlashEndTimesBlk;
taskData.stimParams = 0;
taskData.RFmap = 0;
taskData.eyeData = tmpEye;
taskData.screenStats = screenStats;
taskData.eyeCal.PixelsPerDegree = MLConfig.PixelsPerDegree;
taskData.eyeCal.origin = MLConfig.EyeTransform{1,2}.origin;
taskData.eyeCal.gain = MLConfig.EyeTransform{1,2}.gain;

end
%

function translationTable = updateTranslationTable(translationTable)
%Function takes the translation table drawn from the MonkeyLogic file, and
%replaces entries with Old names with their new names.

RosettaStone = {...
{'monkeyObserving_1191.avi', 'monkeyObserving_1171.avi'}
{'monkeyObserving_1192.avi', 'monkeyObserving_1172.avi'}
{'monkeyObserving_1193.avi', 'monkeyObserving_1173.avi'}
{'monkeyObserving_1194.avi', 'monkeyObserving_1174.avi'}
{'monkeyObserving_1195.avi', 'monkeyObserving_1175.avi'}

{'monkeyGoalDir_1201.avi', 'monkeyGoalDir_1101.avi'}
{'monkeyGoalDir_1202.avi', 'monkeyGoalDir_1102.avi'}
{'monkeyGoalDir_1203.avi', 'monkeyGoalDir_1103.avi'}
{'monkeyGoalDir_1204.avi', 'monkeyGoalDir_1104.avi'}
{'monkeyGoalDir_1205.avi', 'monkeyGoalDir_1105.avi'}

{'humanGoalDir_1201.avi', 'humanGoalDir_1101.avi'}
{'humanGoalDir_1202.avi', 'humanGoalDir_1102.avi'}
{'humanGoalDir_1203.avi', 'humanGoalDir_1103.avi'}
{'humanGoalDir_1204.avi', 'humanGoalDir_1104.avi'}
{'humanGoalDir_1205.avi', 'humanGoalDir_1105.avi'}

{'Cut_obj_interact1_1.avi', 'objects_2101.avi'}; ...
{'Cut_obj_interact1_2.avi', 'objects_2102.avi'}; ...
{'Cut_Order1Interaction_1.avi', 'monkeyChasing_1111.avi'}; ...
{'Cut_Order1Interaction_2.avi', 'monkeyFighting_1121.avi'}; ...
{'Cut_Order1Interaction_3.avi', 'monkeyGrooming_1141.avi'}; ...
{'Cut_Order1Interaction_4.avi', 'monkeyGrooming_1142.avi'}; ...
{'Cut_Order1Interaction_5.avi', 'monkeyMounting_1131.avi'}; ...
{'Cut_Order2Interaction_1.avi', 'monkeyFighting_1124.avi'}; ...
{'Cut_Order2Interaction_2.avi', 'monkeyGrooming_1146.avi'}; ...
{'Cut_Order2Interaction_3.avi', 'monkeyMounting_1133.avi'}; ...
{'Cut_Order2Interaction_4.avi', 'monkeyGrooming_1147.avi'}; ...
{'Cut_Order2Interaction_5.avi', 'monkeyMounting_1136.avi'}; ...
{'Cut_Order2LandscapesFull_1.avi', 'landscape_4001.avi'}; ...
{'Cut_Order2LandscapesFull_2.avi', 'landscape_4002.avi'}; ...

{'Dephased_Order1Interactionchasing.avi', 'scramble_3001.avi'}; ...
{'Dephased_Order1Interactionfighting.avi', 'scramble_3002.avi'}; ...

{'Dephased_Order1Interaction_1.avi', 'scramble_3001.avi'}; ...
{'Dephased_Order1Interaction_2.avi', 'scramble_3002.avi'}; ...
{'Dephased_Order1Interaction_3.avi', 'scramble_3003.avi'}; ...
{'Dephased_Order1Interaction_4.avi', 'scramble_3004.avi'}; ...
{'Dephased_Order1Interaction_5.avi', 'scramble_3005.avi'}; ...

{'Dephased_Order2Interaction_1.avi', 'scramble_3006.avi'}; ...
{'Dephased_Order2Interaction_2.avi', 'scramble_3007.avi'}; ...
{'Dephased_Order2Interaction_3.avi', 'scramble_3008.avi'}; ...
{'Dephased_Order2Interaction_4.avi', 'scramble_3009.avi'}; ...
{'Dephased_Order2Interaction_5.avi', 'scramble_30010.avi'}; ...

{'Cut_human_interaction_1.avi', 'humanChasing_1111.avi'}; ...
{'Cut_human_interaction_2.avi', 'humanFollowing_1161.avi'}; ...
{'Cut_human_interaction_3.avi', 'humanGrooming_1141.avi'}; ...
{'Cut_human_interaction_4.avi', 'humanFighting_1121.avi'}; ...
{'Cut_human_interaction_5.avi', 'humanMounting_1131.avi'};...

{'Order1GoalLeft_1_Order1GoalRight_1.avi', 'monkeyGoalDir_1101.avi'};...
{'Order1GoalLeft_2_Order1GoalRight_2.avi', 'monkeyGoalDir_1102.avi'};...
{'Order1GoalLeft_3_Order1GoalRight_3.avi', 'monkeyGoalDir_1103.avi'};...
{'Order1GoalLeft_4_Order1GoalRight_4.avi', 'monkeyGoalDir_1104.avi'};...
{'Order1GoalLeft_5_Order1GoalRight_5.avi', 'monkeyGoalDir_1105.avi'};...

{'Order1MonkLeft_1_Order1MonkRight_1.avi', 'monkeyIdle_1301.avi'};...
{'Order1MonkLeft_2_Order1MonkRight_2.avi', 'monkeyIdle_1302.avi'};...
{'Order1MonkLeft_3_Order1MonkRight_3.avi', 'monkeyIdle_1303.avi'};...
{'Order1MonkLeft_4_Order1MonkRight_4.avi', 'monkeyIdle_1304.avi'};...
{'Order1MonkLeft_5_Order1MonkRight_5.avi', 'monkeyIdle_1305.avi'};...

{'human_alone1_1_1_human_alone2_1_1.avi', 'humanIdle_1301.avi'};...
{'human_alone1_1_2_human_alone2_1_2.avi', 'humanIdle_1302.avi'};...
{'human_alone1_1_3_human_alone2_1_3.avi', 'humanIdle_1303.avi'};...
{'human_alone1_1_4_human_alone2_1_4.avi', 'humanIdle_1304.avi'};...
{'human_alone1_1_5_human_alone2_1_5.avi', 'humanIdle_1305.avi'};...

{'human_goal1_1_4_human_goal2_1_4.avi', 'humanGoalDir_1101.avi'};
{'human_goal1_1_5_human_goal2_1_5.avi', 'humanGoalDir_1102.avi'};
{'human_goal1_1_1_human_goal2_1_1.avi', 'humanGoalDir_1103.avi'};
{'human_goal1_1_2_human_goal2_1_2.avi', 'humanGoalDir_1104.avi'};
{'human_goal1_1_3_human_goal2_1_3.avi', 'humanGoalDir_1105.avi'};

{'Cut_Order1Chasing1.avi', 'monkeyChasing_1113.avi'}; ...
{'Cut_Order1Chasing2.avi', 'monkeyChasing_1112.avi'}; ...
{'Cut_Order1Chasing3.avi', 'monkeyChasing_1111.avi'}; ...
{'Cut_Order1Chasing4.avi', 'monkeyChasing_1114.avi'}; ...
{'Cut_Order1Chasing5.avi', 'monkeyChasing_1115.avi'}; ...

{'Cut_Order1Fighting1.avi', 'monkeyFighting_1123.avi'}; ...
{'Cut_Order1Fighting2.avi', 'monkeyFighting_1122.avi'}; ...
{'Cut_Order1Fighting3.avi', 'monkeyFighting_1121.avi'}; ...
{'Cut_Order1Fighting4.avi', 'monkeyFighting_1124.avi'}; ...
{'Cut_Order1Fighting5.avi', 'monkeyFighting_1125.avi'}; ...

{'Cut_Order1Mounting1.avi', 'monkeyMounting_1131.avi'}; ...
{'Cut_Order1Mounting2.avi', 'monkeyMounting_1132.avi'}; ...
{'Cut_Order1Mounting3.avi', 'monkeyMounting_1133.avi'}; ...
{'Cut_Order1Mounting4.avi', 'monkeyMounting_1134.avi'}; ...
{'Cut_Order1Mounting5.avi', 'monkeyMounting_1135.avi'}; ...

{'monkeyGrooming_11413.avi', 'monkeyGrooming_11413A.avi'}; ...
{'monkeyGrooming_11413_C1.avi', 'monkeyGrooming_11413A_C1.avi'}; ...

};

oldHalf = cell(size(RosettaStone));
newHalf = cell(size(RosettaStone));

for event_i = 1:length(RosettaStone)
  oldHalf{event_i} = RosettaStone{event_i}{1};
  newHalf{event_i} = RosettaStone{event_i}{2};
end
 
% {'Cut_Order1Fighting1.avi', 'socialInteraction', 'fighting','highMotion'}; ...
% {'Cut_Order1Fighting2.avi', 'socialInteraction', 'fighting','highMotion'}; ...
% {'Cut_Order1Fighting3.avi', 'socialInteraction', 'fighting','highMotion'}; ...
% {'Cut_Order1Fighting4.avi', 'socialInteraction', 'fighting','highMotion'}; ...
% {'Cut_Order1Fighting5.avi', 'socialInteraction', 'fighting','highMotion'}; ...
% {'Cut_Order1Grooming5.avi', 'socialInteraction', 'grooming'}
% 
% {'Dephased_Order1Interactionchasing.avi', 'scramble', 'highMotion'}; ...
% {'Dephased_Order1Interactionfighting.avi', 'scramble', 'highMotion'}; ...
% };

for ii = 1:length(translationTable)
  if any(strcmp(translationTable{ii}, oldHalf))
    translationTable{ii} = newHalf{strcmp(translationTable{ii}, oldHalf)};
  end
end

end


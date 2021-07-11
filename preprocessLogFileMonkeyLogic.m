function [ taskData, stimTiming ] = preprocessLogFileMonkeyLogic(logfile, taskTriggers, diodeTriggers, params)
%% Reads in Blackrock
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
% FOR TESTING:
% load('monkeyLogicPreProcessVars');

%Parse the Log file
disp('Loading MonkeyLogic log file');
assert(logical(exist(logfile,'file')),'The logfile you requested does not exist.');
[data, MLConfig, TrialRecord] = mlread(logfile);

% Check if there are duplicates (in this case, just check for more than 1
% 10 in the code.
if sum(data(1).BehavioralCodes.CodeNumbers == 10) > 1
  data = removeDups(data);
end

% assert(length(unique(TrialRecord.ConditionsPlayed)) < ceil((length(TrialRecord.ConditionsPlayed))/2) , 'MonkeyLogic file reports each condition was not repeated at least 2 times')

% Find out experiment length
tmp = datevec(TrialRecord.TaskInfo.EndTime - TrialRecord.TaskInfo.StartTime);
runTimeMins = tmp(end-1) + (tmp(end)/60);

%Construct the translation table from the logfile
%Find the Conditions that need to be connected to codes, and loop through
%the trials until each of those numbers has a stim attached.
allConditions = unique(TrialRecord.ConditionsPlayed); %Pull unique members of this list
translationTable = repmat({''}, [max(allConditions),1]); %Initialize translation table

% Find out which TaskObject is the Fixation point (this changes)
trialHolder = data(1);
stimName = {trialHolder.TaskObject.Attribute.Name}; %Pull the string containing the stimulus name.
stimInd = find(~strcmp(stimName, 'Fixation Point'));

% Extract the scale factor
try
  scaleFactor = [data(1).ObjectStatusRecord.SceneParam.Scale];
  scaleFactor = scaleFactor(stimInd);
catch
  error('SCALE FACTOR NOT FOUND!') % Temporary, just for the first run.
  scaleFactor = 1;
end

% Use that to generate a table of stimuli
for trial_ind = 1:length(data)
  trialHolder = data(trial_ind);
  stimName = trialHolder.TaskObject.Attribute(stimInd).Name; %Pull the string containing the stimulus name.
  translationTable{trialHolder.Condition} = stimName(~isspace(stimName));
end

if strcmp(TrialRecord.TaskInfo.Stimuli{1}(2), ':') % new version of ML might show full path?
  logFileStimTable = cell(length(TrialRecord.TaskInfo.Stimuli),1);
  for ii = 1:length(logFileStimTable)
    [~, B, C] = fileparts(TrialRecord.TaskInfo.Stimuli{ii});
    logFileStimTable{ii} = [B C];
  end
else
  logFileStimTable = strtrim((TrialRecord.TaskInfo.Stimuli)); %Clean up the stim table in the logfile.
end
%Check to see if what you collected by cycling through the stim matches
%this table.
% assert(sum(strcmp(logFileStimTable,translationTable)) == length(translationTable), 'Stim table from MonkeyLogic doesnt match collected table');

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
codesAndTimes = [data.BehavioralCodes];
rwdTimes = [data.RewardRecord];

%initialize everything
[fixOnset, taskEventStartTimesLog, taskEventEndTimesLog, juiceOnTimesLog, juiceOffTimesLog, stimFramesLost] = deal(zeros(length(correctTrialArray), 1));

%Collect trialEventIDs.
taskEventIDsLog = TrialRecord.ConditionsPlayed';
  
for ii = 1:length(mklTrialStarts)
  trialData = codesAndTimes(ii);
  
  % Fixation cue onset.
  fixOnset(ii) = trialData.CodeTimes(trialData.CodeNumbers == fixCueMarker);
  
  % Identify stimulus start and stop times (if they take place).
  if ~isempty(trialData.CodeTimes(trialData.CodeNumbers == stimStartMarker))
    taskEventStartTimesLog(ii) = mklTrialStarts(ii) + trialData.CodeTimes(trialData.CodeNumbers == stimStartMarker);
    taskEventEndTimesLog(ii) = mklTrialStarts(ii) + trialData.CodeTimes(trialData.CodeNumbers == stimEndMarker);
  else
    taskEventStartTimesLog(ii) = nan;
    taskEventEndTimesLog(ii) = nan;
  end
  
  % Reward delivery times (if it happens).  
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

assert(~isempty(packetData), 'The Blackrock Digital inputs are empty. Digital inputs may have not been plugged in.')

% Check for 8s (in some runs, the 1 bit seems to be malfunctioning, other
% times something more wide ranging). Other numbers are signs of other
% malfunctions.
if any(packetData == 8) || any(packetData == 22) || any(packetData == 39)
  warning('Every trial doesnt have a condition number - this may be due to changes in when the marker is sent. Seeing if we can use MKL data to recreate')
  % This code will now take the monkeyLogic signals as ground truth
  [packetData, packetTimes] = replacePacketData(data, TrialRecord, packetData, packetTimes);
end

%Code to get rid of any markers prior to the first trial beginning, or after the final trial end.
%a defense against double 9's
trueStart = find(packetData == trialStartMarker);
if length(trueStart) > 1
  trueStartTrials = (diff(trueStart) > 1);
  trueStart = trueStart(trueStartTrials);
end
trueStart = trueStart(1);
trueEnd = find(packetData == trialEndMarker, 1, 'last');

packetData = packetData(trueStart:trueEnd);
packetTimes = packetTimes(trueStart:trueEnd);

%to later figure out how to shape data, we need to know what we saw and
%got rid of.

%Below are things which should happen every trial
trialStartInds = find(packetData == trialStartMarker);
trialEndInds = find(packetData == trialEndMarker);
condNumbers = packetData(packetData > 100)-100;
stimFixEndInds = find(packetData == stimEndMarker);

%Now, use the correct trial array from MKL to pick out only the correct
%trials from the whole array. Assign them them to the correct vectors.
[trialStartTimesBlk, taskEventStartTimesBlk, taskEventEndTimesBlk, juiceOnTimesBlk, juiceOffTimesBlk, rewardTrialInds] = deal(nan(size(trialStartInds)));

% Identify which trials were rewarded by looking for the correct series
% of eventMarkers.
trialIntervals = [find(packetData == 9), find(packetData == 18)]';
for trial_i = 1:length(juiceOnTimesBlk)
  trial_int = trialIntervals(:, trial_i);
  rewardTrialInds(trial_i) = any(packetData(trial_int(1):trial_int(2)) == rewardMarker);
end

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
      case {rewardMarker, trialEndMarker}
        errorArray(ii) = 0;
        found = 1;
      otherwise
        stepsAhead = stepsAhead + 1;
    end
  end
end

%packetData is directly referenced Below
stimOnTrials = errorArray ~= 4;

taskEventIDsBlk = condNumbers;
trialStartTimesAll = packetTimes(packetData == trialStartMarker);
trialStartTimesBlk(stimOnTrials) = trialStartTimesAll(stimOnTrials);

% Save all the trials fixation dot (happens on 100% of trials).
fixOnsetTimesBlk = packetTimes(packetData == fixCueMarker);

taskEventStartTimesBlk(stimOnTrials) = packetTimes(packetData == stimStartMarker); %Stores every trial where stim started (no fail during fix);
taskEventEndTimesBlk(stimOnTrials) = packetTimes(stimFixEndInds(stimOnTrials));

% Juice processing needs to account for paradigms where reward doesn't
% come w/ every trial.
juiceOnTimesBlk(rewardTrialInds == 1) = packetTimes(packetData == rewardMarker);
juiceOffTimesBlk(rewardTrialInds == 1) = packetTimes(trialEndInds(rewardTrialInds == 1));

%use the condition number to fill out the taskEventsIDs using the
%translation table.
taskEventIDs = translationTable(condNumbers);

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
  % Recalculate lengths
  LogLen = length(taskEventIDsLog);
  BlkLen = length(taskEventIDsBlk);
  
  % If The lengths are still different, one of them stopped early.
  
  if LogLen < BlkLen %monkeyLogic stopped storing things early.
    %Chop the Blackrock events down to the size of what Mkl has stored
    taskEventIDsBlk = taskEventIDsBlk(1:LogLen);
    if sum(taskEventIDsBlk == taskEventIDsLog) == length(taskEventIDsLog) % If They are now a match
      trialStartTimesBlk = trialStartTimesBlk(1:LogLen);
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
  load(fullfile(frameMotionFile(1).folder, frameMotionFile(1).name), 'frameMotionData');
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
  
%% If eventData file is present, search it, and generate composite subEvents.
% Code below generates a table used by neuroGLM function.
eventDataFile = dir([params.stimDir '/**/eventData.mat']);
if ~isempty(eventDataFile)
  load(fullfile(eventDataFile(1).folder, eventDataFile(1).name), 'eventData');
  if ~isempty(intersect(eventData.Properties.RowNames, translationTable))
    eventDataRun = eventData(translationTable, :);
    eventsInEventData = eventDataRun.Properties.VariableNames;
    
    % Generate an output table
    eventDataStarts = cell2table(cell(size(eventDataRun)));
    eventDataStarts.Properties.RowNames = translationTable;
    eventDataStarts.Properties.VariableNames = eventsInEventData;
    
    for stim_i = 1:length(translationTable)
      % Converting frames to times.
      timePerFrame = tmpFrameMotionData(stim_i).timePerFrame;
      
      for event_i = 1:length(eventsInEventData)
        events2pull = eventDataRun{translationTable{stim_i}, eventsInEventData{event_i}}{1};
        if ~isempty(events2pull)
          eventDataStarts{translationTable{stim_i}, eventsInEventData{event_i}} = {events2pull.startFrame*timePerFrame};
        end
      end
      
    end
    
    % Initialize output struct, remove empty columns, load variables
    taskData.eventData = struct();
    column2keepInd = logical(sum(~cellfun('isempty', table2cell(eventDataStarts))));
    taskData.eventData.events = eventsInEventData(column2keepInd);
    taskData.eventData.eventDataStarts = eventDataStarts(:,column2keepInd);
  end
end

%% Reward processing
% Generate a vector of reward times relative to stimulus onset.
rewardTimePerTrial = juiceOnTimesBlk - taskEventStartTimesBlk;

% Identify if rewards are missing on correct trials (rewardParadigm).
% Missing reward on completed trials
missingReward = rewardTimePerTrial(TrialRecord.TrialErrors == 0);

if any(isnan(missingReward))
  rewardParadigm = true;
else
  rewardParadigm = false;
end

%% Jan 2021 - Add paradigm tag, based on stimuli present

if any(contains(translationTable, 'headTurnCon'))
  taskData.paradigm = 'headTurnCon';
elseif any(contains(translationTable, 'headTurnIso'))
  taskData.paradigm = 'headTurnIso';
elseif any(contains(translationTable, 'monkey'))
  taskData.paradigm = 'naturalSocial';
elseif any(contains(translationTable, 'Barney'))
  taskData.paradigm = 'familiarFace';
end

%Find stimulus timing range - in ML, we have no range within valid trials.
switch taskData.paradigm
  case {'headTurnIso', 'familiarFace'}
    stimTiming.shortest = 1000; %Making this 0 for now
    stimTiming.longest = 1000; %setting this to 2800 for now, but maybe save as an editable variable.
  case {'headTurnCon', 'naturalSocial'}
    stimTiming.shortest = 2800; %Making this 0 for now
    stimTiming.longest = 2800; %setting this to 2800 for now, but maybe save as an editable variable.
  otherwise
    error('stimTiming Missing');
end
stimTiming.ISI = MLConfig.InterTrialInterval;

%% March 2021
% Add in fixation time per trial - this varies due to punishment dynamic.
% It is 95%+ the typical value, but other times it is more.

fixTime = nan(length(data),1);
for ii = 1:length(data)
  fixTime(ii) = data(ii).ObjectStatusRecord.SceneParam(1).AdapterArgs{3}{2,2};
end

% April 2021
[~, fName, ~] = fileparts(TrialRecord.DataFile);
dateSubj = extractBefore(fName, '00');

if strcmp(dateSubj, '20201115Mo')
  gridHole = 'C9';
  recDepth = '0';
else
  gridHole = data(1).VariableChanges.gridHole;
  recDepth = data(1).VariableChanges.recordingDepth;
end

% Merge IDs for the animated paradigms
taskEventIDsMerged = eventIDMerge(taskEventIDs);

%% Create a fixation time per stimulus
taskEventFixDurBlk = fixOnsetTimesBlk;
assert(length(taskEventFixDurBlk) == length(taskEventIDs), 'Problem w/ fixation times')

%% Output
%Adding random numbers to these - they aren't relevant for my current task,
%nor are they directly recorded by MKL.
firstNonNaN = find(~isnan(taskEventStartTimesBlk),1);
fixSpotFlashStartTimesBlk = taskEventStartTimesBlk(firstNonNaN);
fixSpotFlashEndTimesBlk = taskEventEndTimesBlk(firstNonNaN);
fixationInTimesBlk = taskEventStartTimesBlk(firstNonNaN);
fixationOutTimesBlk = taskEventEndTimesBlk(firstNonNaN);

% finally, build the output structure:
% Note: **If something is missing in runAnalyses** - likely due to taskData
% being overwritten following excludeTrial function.

taskData.taskDataSummary.TrialRecord = TrialRecord;
taskData.errorArray = errorArray;
taskData.taskEventIDsMerged = taskEventIDsMerged;
taskData.taskEventIDs = taskEventIDs;
taskData.taskEventList = translationTable;
taskData.frameMotionData = tmpFrameMotionData';
taskData.stimFramesLost = stimFramesLost;
taskData.scaleFactor = scaleFactor;
taskData.gridHole = gridHole;
taskData.recDepth = recDepth;
taskData.rewardParadigm = rewardParadigm;

% Timing variables
taskData.runTime = runTimeMins;
taskData.taskEventStartTimes = taskEventStartTimesBlk;
taskData.taskEventEndTimes = taskEventEndTimesBlk;
taskData.taskEventFixDur = fixOnsetTimesBlk;
taskData.rewardTimePerTrial = rewardTimePerTrial;
taskData.trialStartTimesMkl = mklTrialStarts';
taskData.logVsBlkModel = logVsBlkModel;
taskData.fixTime = fixTime;
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
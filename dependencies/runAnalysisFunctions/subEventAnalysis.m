function [subEventSigStruct, specSubEventStruct, selTable] = subEventAnalysis(eyeBehStatsByStim, spikesByChannel, taskData, eventIDs, onsetsByEventStim, subEventParams, selTable, psthParams, figStruct)
% subEventAnalysis
% Description - looks through spikesByEvent, calculates PSTH for activity
% aligned to a specific event as well as a null distribution from the
% remaining stimuli. Produces statistical tests for the traces.
% Parameters
% eyeBehStatsbyStim - a structure generated by earlier eye analysis
% functions. {stim}{trial}
% spikesByChannel - a relatively raw stream of spikes, produced in
% processRun. time stamps are used to pull spikes from it.
% subEventParams - has path to eventData, stimDir, and variables used to
% plot PSTHes.

saccadeSplit = true;           % Split subEvent in the stimulus depending on whether they were followed by a saccade.
keepNonSplit = true;           % if splitting subevents by saccade or non-saccade, decide if you keep the original. 

% Unpack Variables, some preprocessing. Generate an index for individual
% trials that take place during the run, matching them to the stimuli as
% present in the taskData.eventIDs field.
frameMotionData = taskData.frameMotionData;
taskEventIDs = taskData.taskEventIDs;
% specSubEvent = subEventParams.specSubEvent;
specSubEvent = 0;
saccadeWindow = 200;    % The amount of time after a subEvent during which a saccade must occur to classify the trial as 'saccade'.

taskEventIDInd = zeros(length(taskEventIDs),1);
for ev_i = 1:length(eventIDs)
  taskEventIDInd = taskEventIDInd + (strcmp(taskEventIDs, eventIDs{ev_i}) * ev_i);
end
assert(max(taskEventIDInd) == length(eventIDs))

% Event Type 1 - Find events which take place in the stimuli presented,
% defined in the eventData file.
eventDataFile = dir([subEventParams.stimDir '/**/eventData.mat']);
% If we find a file, load it.
load(fullfile(eventDataFile(1).folder, eventDataFile(1).name), 'eventData');
eventsInData = intersect(eventIDs, eventData.Row);

% If all the events in the run are in the eventData
eventDataRun = eventData(eventIDs, :);
emptyEntryInd = cellfun(@(x) size(x,1) ~= 0, table2cell(eventDataRun));
presentEventInd = any(emptyEntryInd);
emptyEntryInd = emptyEntryInd(:,presentEventInd);
eventsInEventData = eventDataRun.Properties.VariableNames(presentEventInd);
eventStimTable = cell2table(cell(0,4), 'VariableNames', {'eventName', 'stimName', 'startFrame', 'endFrame'});

% If our run has eventData specific events
subEventNames = {};
subEventNames = [subEventNames; eventsInEventData'];

% Cycle through events, generating table to be used as reference for
% collecting event specific times (has startTime, endTime, stim name).
for event_i = 1:length(eventsInEventData)
  
  subEventRunData = eventDataRun(:, eventsInEventData{event_i});
  subEventRunData = subEventRunData(emptyEntryInd(:, event_i),:);
  stimSourceArray = subEventRunData.Properties.RowNames;
  
  for stim_i = 1:length(stimSourceArray)
    
    %Identify time per frame, shift event frames to event times.
    frameMotionDataInd = strcmp({frameMotionData.stimVid}, stimSourceArray(stim_i));
    
    if ~any(frameMotionDataInd)
      stimVid = dir([subEventParams.stimDir '/**/' stimSourceArray{stim_i}]);
      vidObj = VideoReader(fullfile(stimVid(1).folder, stimVid(1).name));
      msPerFrame = (vidObj.Duration/vidObj.NumFrames) * 1000;
      clear vidObj
    else
      msPerFrame = frameMotionData(frameMotionDataInd).timePerFrame;
    end
    
    % Add events to larger table
    eventTable = subEventRunData{stimSourceArray(stim_i),:}{1};
    eventTable.startFrame = eventTable.startFrame * msPerFrame;
    eventTable.endFrame = eventTable.endFrame * msPerFrame;
    entryCount = size(eventTable,1);
    entryStarts = table2cell(eventTable(:,'startFrame'));
    entryEnds =  table2cell(eventTable(:,'endFrame'));
    
    eventStimTable = [eventStimTable; [repmat(eventsInEventData(event_i), [entryCount, 1]), repmat(stimSourceArray(stim_i), [entryCount, 1]), entryStarts, entryEnds]];
    
  end
end

% Generate a vector of spikeTimes in the conventional structure.
[onsetsByEvent, onsetsByEventNull, stimSourceByEvent] = deal(cell(size(eventsInEventData')));
stimEventMat = zeros(length(eventIDs), subEventParams.stimPlotParams.psthImDur, length(eventsInEventData));

% Event Type 1 - subEvents hand labeled, which occur at fixed times in the
% stimulus.
for event_i = 1:length(eventsInEventData)
  % Identify spaces with the events in the presented stimuli.
  eventData = eventStimTable(strcmp(eventStimTable.eventName, eventsInEventData{event_i}),2:end);
  stimWithEvent = unique(eventData.stimName);
  
  % Identify trials of stimuli without this event, for null sampling
  stimWEvent = ismember(eventIDs, stimWithEvent);
  stimWOEventInd = find(~stimWEvent);
  stimWEventInd = find(stimWEvent)';
  taskEventNullSampling = ismember(taskEventIDInd, stimWOEventInd);
  
  [eventStartTimes, eventStartTimesNull, eventStimSource] = deal([]);
  for stim_i = stimWEventInd
    % Get the appropriate start and end frames, convert to
    stimTrialInd = (taskEventIDInd == stim_i);
    stimEventData = eventData(strcmp(eventData.stimName, eventIDs{stim_i}), :);
    
    % Cycle through and find event occurance times
    for tab_i = 1:size(stimEventData,1)
      startTime = round(stimEventData.startFrame(tab_i));
      eventTitle = {sprintf('%s | %s - %d', eventsInEventData{event_i}, extractBefore(eventIDs{stim_i}, '.'), startTime)};
      
      % Label the stimEventMat to highlight sampled region.
      startEventInd = round(stimEventData.startFrame(tab_i));
      endEventInd = round(stimEventData.endFrame(tab_i));
      stimEventMat(stim_i, startEventInd:endEventInd, event_i) = deal(1);
      
      % Find out when the stimulus happens, and add these structures to
      % the larger ones
      eventStartTimes = [eventStartTimes; taskData.taskEventStartTimes(stimTrialInd) + startTime];
      eventStimSource = [eventStimSource; repmat(eventTitle, [length(taskData.taskEventStartTimes(stimTrialInd)), 1])];
      
      % Generate null times by using other trials where the event doesn't
      % take place.
      nullTimesStart = taskData.taskEventStartTimes(taskEventNullSampling) + startTime;
      eventStartTimesNull = [eventStartTimesNull; nullTimesStart];
      
    end
    %eventStimSource = [eventStimSource; repmat(stimWithEvent(stim_i), [length(eventStartTimes), 1])];
  end
  onsetsByEvent{event_i} = eventStartTimes;
  onsetsByEventNull{event_i} = eventStartTimesNull;
  stimSourceByEvent{event_i} = eventStimSource;
end

% Event Type 1.5 - If we want, we can parse the events into individual instances, incase
% there is some effect taking place for 'head turns during idle'
if specSubEvent
  
  % Create larger vectors
  onsetsByEventAll = vertcat(onsetsByEvent{:});
  stimSourceByEventAll = vertcat(stimSourceByEvent{:});
  uniqueSubEvents = unique(stimSourceByEventAll);

  % Grab times of the individual events.
  [onsetsByEventSpec, onsetsByEventSpecNull] = deal(cell(size(uniqueSubEvents)));
  onsetsByEventSpecCounts = zeros(size(uniqueSubEvents));
  for ii = 1:length(uniqueSubEvents)
    onsetsByEventSpec{ii} = onsetsByEventAll(strcmp(stimSourceByEventAll, uniqueSubEvents{ii}));
    onsetsByEventSpecNull(ii) = onsetsByEventNull(strcmp(subEventNames, extractBefore(uniqueSubEvents{ii}, ' |')));
    onsetsByEventSpecCounts(ii) = length(onsetsByEventSpec{ii});
    
  end
  
  onsetsByEvent = [onsetsByEvent; onsetsByEventSpec];
  onsetsByEventNull = [onsetsByEventNull; onsetsByEventSpecNull];
  subEventNames = [subEventNames; uniqueSubEvents];
end

% Event Type 2 - Generate/Organize Event times for eye movements here, concatonate onto
% the 'onsetsByEvent' cue.
if ~isempty(eyeBehStatsByStim)
  trialsPerStim = cellfun('length', eyeBehStatsByStim);
  stimuliStartTimes = onsetsByEventStim;     % Find all the start times for the stimulus
  
  [stimSaccadeArrays, stimSaccadeNullArrays] = generateAdvSaccadeNullTimes(eyeBehStatsByStim, psthParams);
  
  [blinkTimes, saccadeTimes, saccadeNonStimTimes, saccadeNullTimes, saccadeNonStimNullTimes] = deal([]);                          % Initialize vectors
  for stim_i = 1:length(eyeBehStatsByStim)
    for trial_i = 1:trialsPerStim(stim_i)
      % Find Absolute eye event start times.
      eyeTrial = eyeBehStatsByStim{stim_i}{trial_i};
      
      if ~isempty(eyeTrial.blinktimes)
        blinkTimes = [blinkTimes; (eyeTrial.blinktimes(1,:)' + stimuliStartTimes{stim_i}(trial_i))];
      end
      
      if ~isempty(eyeTrial.saccadetimes)
        saccadeBins = find(stimSaccadeArrays{stim_i}(trial_i,:)) - psthParams.psthPre;
        saccadeNullBins = find(stimSaccadeNullArrays{stim_i}(trial_i,:)) - psthParams.psthPre;
        nonStimSaccInd = saccadeBins < 0 | saccadeBins > psthParams.psthImDur;
        nonStimNullSaccInd = saccadeNullBins < 0 | saccadeNullBins > psthParams.psthImDur;
        
        saccadeTimes = [saccadeTimes; (saccadeBins' + stimuliStartTimes{stim_i}(trial_i))];
        saccadeNonStimTimes = [saccadeNonStimTimes; (saccadeBins(nonStimSaccInd)' + stimuliStartTimes{stim_i}(trial_i))];
        saccadeNullTimes = [saccadeNullTimes; (saccadeNullBins' + stimuliStartTimes{stim_i}(trial_i))];
        saccadeNonStimNullTimes = [saccadeNonStimNullTimes; (saccadeNullBins(nonStimNullSaccInd)' + stimuliStartTimes{stim_i}(trial_i))];
      end
      
    end
  end
  
  % Generate pre-saccade times - *Pre-sacc test period is -200 to 0, so the
  % event times are the same as for saccades.
  preSaccadeTimes = saccadeTimes;% - subEventParams.preSaccOffset;
  preSaccadeNonStimTimes = saccadeNonStimTimes;% - subEventParams.preSaccOffset;
  preSaccadeNullTimes = saccadeNullTimes - subEventParams.preSaccOffset;
  preSaccadeNonStimNullTimes = saccadeNonStimNullTimes - subEventParams.preSaccOffset;
  
  % Sort these lists for the next processing steps
%   eyeEventTimes = {blinkTimes; saccadeTimes; preSaccadeTimes}; % Putting saccadeTimes in twice for pre-saccade times.
  eyeEventTimes = {blinkTimes; saccadeTimes; preSaccadeTimes; saccadeNonStimTimes; preSaccadeNonStimTimes}; % Putting saccadeTimes in twice for pre-saccade times.
  
  % Generate Null times, add a number (100 - 1000), shuffle them and check
  % if its close to a number on the list. If not, its a null value.
  
  blinkTimesNull = createNullScrambleTimes({blinkTimes}, 500, 300);
  eyeEventNullTimes = [blinkTimesNull; saccadeNullTimes; preSaccadeNullTimes; saccadeNonStimNullTimes; preSaccadeNonStimNullTimes];
  
  if saccadeSplit
    % Event Type 2.5 - subEvents occuring in the stimuli, split by whether a
    % subsequent saccade took place.
    saccadeTimes = sort(saccadeTimes);
    subEventsNew = cell(length(subEventNames), 2);
    subEventNull = [onsetsByEventNull, onsetsByEventNull];
    subEventNamesNew = [strcat(subEventNames, '_sacc'), strcat(subEventNames, '_saccNone')];
    for event_i = 1:length(subEventNames)
      
      % Pull event times
      eventOnsets = onsetsByEvent{event_i};
      
      % See how they compare to previously collect saccade times
      eventOnsets2Check = repmat(saccadeTimes, [1, length(eventOnsets)]);
      eventOnsetsOffsetBySacc = eventOnsets2Check - eventOnsets';
      
      % If any number is between 0 and the previously defined window, it means
      % a saccade occured close enough after.
      eventsWSaccades = any(eventOnsetsOffsetBySacc > 0 & eventOnsetsOffsetBySacc < saccadeWindow, 1);
      
      % Store them as distinct events
      subEventsNew{event_i, 1} = eventOnsets(eventsWSaccades);
      subEventsNew{event_i, 2} = eventOnsets(~eventsWSaccades);
      
    end
    
    % Rearrange and overwrite originals
    if keepNonSplit
      onsetsByEvent = [onsetsByEvent; reshape(subEventsNew', [], 1)];
      onsetsByEventNull = [onsetsByEventNull; reshape(subEventNull', [], 1)];
      subEventNames = [subEventNames; reshape(subEventNamesNew', [], 1)];
    else
      onsetsByEvent = reshape(subEventsNew', [], 1);
      onsetsByEventNull = reshape(subEventNull', [], 1);
      subEventNames = reshape(subEventNamesNew', [], 1);
    end
  end
  
    % Remove empty ones - Dont so that subsequent table can stack.
%   keepInd = ~cellfun('isempty', onsetsByEvent);
%   onsetsByEvent = onsetsByEvent(keepInd);
%   subEventNames = subEventNames(keepInd);
%   onsetsByEventNull = onsetsByEventNull(keepInd);
  
  % Add Previously generate eye focused events to array
  subEventNames = [subEventNames; ["blinks"; "saccades"; "pre-saccades"; "saccadesNonStim"; "pre-saccadesNonStim"]]; % Presaccade will have the same times as saccades, but the comparison window for the t test will stretch back.
  onsetsByEvent = [onsetsByEvent; eyeEventTimes];
  onsetsByEventNull = [onsetsByEventNull; eyeEventNullTimes];
  
end

% Quick cleaning step - for saccade split events, you don't want 2 data
% points compared to 1000. for each event, make sure the number of Null
% events is proportional.
rng(1)
eventsCount = cellfun('length', onsetsByEvent) * 10;
eventsNullCount = cellfun('length', onsetsByEventNull);
nullDownSampInd = find(eventsCount < eventsNullCount)';
for ii = nullDownSampInd
  onsetsByEventNull{ii} = onsetsByEventNull{ii}(randperm(eventsNullCount(ii), eventsCount(ii)));
end

% Event Type 3 - add the Reward as an event
if subEventParams.RewardEvent
  
  % Establish times of reward delivery
  rewardTimes = taskData.juiceOnTimes(~isnan(taskData.juiceOnTimes));
  rewardTimesNull = createNullScrambleTimes({rewardTimes}, 300, 200);  % Generate Null times for the rewardTimes event.
    
  % Establish reward anticipation vector based on 'rewardAntTime' param.
  % Reward anticipation testing window is 0 - 200 ms, so times are the same
  % as reward.
  % Because the test window is 0 to rewardAntTime, adding it in just makes it normal reward times vs the precending x ms.
  
  rewardAntTimes = rewardTimes; 
  rewardAntNull = rewardTimes + subEventParams.rewardAntTime;  
  
  % Stick onto larger structures - rewardAnt is t tested for the period 
  subEventNames = [subEventNames; "reward"; "rewardAnt"];
  onsetsByEvent = [onsetsByEvent; rewardTimes; rewardAntTimes];
  onsetsByEventNull = [onsetsByEventNull; rewardTimesNull; rewardAntNull];
  
  % Generate a 2nd reward comparison, contraining reward to to unrewarded trials.
  if taskData.rewardParadigm
    unRewardedTrials = isnan(taskData.juiceOnTimes);
    unRewardedStimStart = taskData.taskEventStartTimes(unRewardedTrials);
    
    % When rewards would have taken place, normally.
    unRewardedStimRewardTimes = unRewardedStimStart + round(nanmean(taskData.rewardTimePerTrial));
    rewardSampling = randi(length(unRewardedStimStart), [length(unRewardedStimStart)*10, 1]);
    unRewardedStimRewardTimes = unRewardedStimRewardTimes(rewardSampling);
    
    % Stick onto larger structures
    subEventNames = [subEventNames; 'rewardAbsent'];
    onsetsByEvent = [onsetsByEvent; rewardTimes];
    onsetsByEventNull = [onsetsByEventNull; unRewardedStimRewardTimes];
    
  end
  
end

% Event Type 4 - add Fixation as an event
if subEventParams.FixEvent
  fixOnTimes = taskData.taskEventFixDur;
  fixOnTimesNull = createNullScrambleTimes({fixOnTimes}, 400, 200);  % Generate Null times for the rewardTimes event.

  subEventNames = [subEventNames; 'fixation'];
  onsetsByEvent = [onsetsByEvent; fixOnTimes];
  onsetsByEventNull = [onsetsByEventNull; fixOnTimesNull];
end

% Event Type 5 - stimulus onset and offset!
subEventParams.StimOnOffEvent = 1;
if subEventParams.StimOnOffEvent
  % Compare stimulus onset to the 200 ms preceding
  stimOnTimes = taskData.taskEventStartTimes;
  stimOnNullTimes = stimOnTimes - 200;
  
  % Compare stimulus offset to the 200 ms preceding it. 
  % for -cohensD, you'd call it a stimulus offset inhibition. 
  stimOffTimes = taskData.taskEventEndTimes;
  stimOffNulllTimes = stimOffTimes - 200;
  
  subEventNames = [subEventNames; 'stimOnset'; 'stimOffset'];
  onsetsByEvent = [onsetsByEvent; stimOnTimes; stimOffTimes];
  onsetsByEventNull = [onsetsByEventNull; stimOnNullTimes; stimOffNulllTimes];
end

% Sample the null times, expand each onsertsByEventNull cell by the
% parameter below (subEventParams.nullSampleMult).
% onsetsByEventNull = cellfun(@(x) repmat(x, [subEventParams.nullSampleMult,1]), onsetsByEventNull, 'UniformOutput', false);
if ~isempty(onsetsByEvent)
  
  subEventParams.refOffset = 0;
  [spikesBySubEvent, spikesEmptyBySubEvent] = alignSpikes(spikesByChannel, onsetsByEvent, ones(length(spikesByChannel),1), subEventParams);
  [spikesBySubEventNull, spikesEmptyBySubEventNull] = alignSpikes(spikesByChannel, onsetsByEventNull, ones(length(spikesByChannel),1), subEventParams);
  
  % follow with putting spikesBySubEvent into calcPSTH.
  if ~subEventParams.spikeTimes
    spikesBySubEventBinned = calcSpikeTimes(spikesBySubEvent, subEventParams.psthParams);
    spikesBySubEventNullBinned = calcSpikeTimes(spikesBySubEventNull, subEventParams.psthParams);
    
    [psthBySubEvent, psthErrBySubEvent] = calcStimPSTH(spikesBySubEventBinned, spikesEmptyBySubEvent, subEventParams.spikeTimes, subEventParams.psthParams, subEventParams);
    [psthBySubEventNull, psthErrBySubEventNull] = calcStimPSTH(spikesBySubEventNullBinned, spikesEmptyBySubEventNull, subEventParams.spikeTimes, subEventParams.psthParams, subEventParams);
  else
    [psthBySubEvent, psthErrBySubEvent] = calcStimPSTH(spikesBySubEvent, spikesEmptyBySubEvent, subEventParams.spikeTimes, subEventParams.psthParams, subEventParams);
    [psthBySubEventNull, psthErrBySubEventNull] = calcStimPSTH(spikesBySubEventNull, spikesEmptyBySubEventNull, subEventParams.spikeTimes, subEventParams.psthParams, subEventParams);
  end
  
  if specSubEvent
    % Remove the subEvents which are specific instances, they need to be processed differently.
    subEventInd = contains(subEventNames, '|');
    
    % Copy the structure generated
    specSubEventPSTH = psthBySubEvent;
    for chan_i = 1:length(specSubEventPSTH)
      for unit_i = 1:length(specSubEventPSTH{chan_i})
        % Go to the base layer, segregate info for specific subEvents and general ones.
        specSubEventPSTH{chan_i}{unit_i} = specSubEventPSTH{chan_i}{unit_i}(subEventInd,:);
        psthBySubEvent{chan_i}{unit_i} = psthBySubEvent{chan_i}{unit_i}(~subEventInd,:);
        psthErrBySubEvent{chan_i}{unit_i} = psthErrBySubEvent{chan_i}{unit_i}(~subEventInd,:);
        psthBySubEventNull{chan_i}{unit_i} = psthBySubEventNull{chan_i}{unit_i}(~subEventInd,:);
        psthErrBySubEventNull{chan_i}{unit_i} = psthErrBySubEventNull{chan_i}{unit_i}(~subEventInd,:);
      end
    end
    
    % Modify other structures where event is at the top layer
    spikesBySubEvent = spikesBySubEvent(~subEventInd);
    subEventNames = subEventNames(~subEventInd);

  end
  
  % Statistics - for every event
  % Method 1 - perform T test on spike rates for the period assigned.
  selArray = ones(size(selTable,1), length(subEventNames));
  diffArray = zeros(size(selTable,1), length(subEventNames));
  
  [spikeCounts, spikeCountsNull] = deal(cell(length(subEventNames), 1));
  chanCount =  length(spikesBySubEvent{1});
  unitCount = cellfun('length', spikesBySubEvent{1});
    
  for event_i = 1:length(subEventNames)
    
    unitInd = 1;
    
    % Find the test window for this event
    eventName = subEventNames{event_i};
    testPeriodInd = strcmp(eventName, subEventParams.possibleEvents);
    if ~any(testPeriodInd)
      testPeriodInd = find(cellfun(@(x) contains(eventName, x), subEventParams.possibleEvents), 1);
    end
    testPeriod = subEventParams.testPeriodPerEvent(testPeriodInd, :);
    assert(~isempty(testPeriod), 'Event %s is missing from array', subEventNames{event_i})
    
    % Count the spikes
    [spikeCounts{event_i}, ~, ~] = spikeCounter(spikesBySubEvent(event_i), testPeriod(1), testPeriod(2));
    [spikeCountsNull{event_i}, ~, ~] = spikeCounter(spikesBySubEventNull(event_i), testPeriod(1), testPeriod(2));
    
    % Initialize test arrays
    if event_i == 1
      [testResults, cohensD] = deal(initNestedCellArray(spikeCounts{1}, 'ones', [1 1], 3));
    end
        
    for chan_i = 1:chanCount
      for unit_i = 1:unitCount(chan_i)
        
        % Collect relevant spikes
        eventSpikes = [spikeCounts{event_i}{chan_i}{unit_i}{1}.rates];
        nullSpikes = [spikeCountsNull{event_i}{chan_i}{unit_i}{1}.rates];
        
        % Perform the test, and store results.
        if isempty(eventSpikes) || isempty(nullSpikes)
          pVal = 1;
        else
          if subEventParams.nonParametric
            [pVal, ~, ~] = ranksum(eventSpikes,  nullSpikes);
          else
            [~, pVal, ~] = ttest2(eventSpikes,  nullSpikes);
          end
        end
        
        % Save the pValue and Cohen's D.
        [selArray(unitInd, event_i), testResults{chan_i}{unit_i}{event_i}] = deal(pVal);
        [diffArray(unitInd, event_i), cohensD{chan_i}{unit_i}{event_i}] = deal(mean(eventSpikes) - mean(nullSpikes));
        
        % Increment units
        unitInd = unitInd + 1;
        
      end
    end
  end
  
  % Plotting
  if subEventParams.genPlots    
    figPerEvent = 1;    % Save a figure for each event. Only save for units since that is what we're using in the paper.
        
    % If we're plotting, we need to remove all the individual subEvent
    % instances, the plot becomes unreadable with them, and we are
    % predominantly interested in them statistically, not to visualize
    dataInd = find(~contains(subEventNames, '|'))';
    
    % Create a plotting array which can be referenced
    [~, plotNameInd] = ismember(subEventNames(dataInd), subEventParams.possibleEvents);
    if ~all(plotNameInd)
      missingLabels = find(plotNameInd == 0)';
      for label_i = missingLabels
        plotNameInd(label_i) = find(cellfun(@(x) contains(subEventNames(label_i), x), subEventParams.possibleEvents), 1);
      end
    end
    plotLabels = subEventParams.possibleEventsPlotNames(plotNameInd);
    
    % Generate an array for referencing subEvents correctly.
    eventsInEventDataPlot = strrep(subEventNames, '_',' ');
    if exist('uniqueSubEvents', 'var') && specSubEvent
      tmp = cellfun(@(x) strsplit(x, '|'), uniqueSubEvents, 'UniformOutput', 0);
      specSubEventNamesSorted = cellfun(@(x) x{1}, tmp, 'UniformOutput', 0);
      %   eventsInEventDataPlot = strrep(subEventNames, '_',' ');
    end
    
    for chan_i = 1:length(psthBySubEvent)
      unitCount = length(psthBySubEvent{chan_i});
      for unit_i = 2:unitCount-1 % skip unsorted and MUA.
        
        % Create per unit figure
        if ~isempty(figStruct.channelUnitNames{chan_i})
          
          chUnitStr = sprintf('%s%s', figStruct.channelNames{chan_i}, figStruct.channelUnitNames{chan_i}{unit_i});
          
          % If one plot for all events, make it here
          if ~figPerEvent
            figTitleSave = sprintf('SubEvent comparison - %s', chUnitStr);
            figH = figure('Name',figTitleSave,'NumberTitle','off','units','normalized', 'position', figStruct.figPos);
            sgtitle(figTitle)
          end
          
          % Cycle through events
          for event_i = 1:length(dataInd)
            
            % If not significant, don't plot
            pValNum = testResults{chan_i}{unit_i}{dataInd(event_i)};
            if pValNum > 0.05
              continue
            end
            
            % Generate the title
            eventName = subEventNames{dataInd(event_i)};
            testPeriodInd = strcmp(eventName, subEventParams.possibleEvents);
            if ~any(testPeriodInd)
              testPeriodInd = find(cellfun(@(x) contains(eventName, x), subEventParams.possibleEvents), 1);
            end
            testPeriod = subEventParams.testPeriodPerEvent(testPeriodInd, :);
            
            % Generate strings for plotting
            pValStr = num2str(pValNum, 2);
            eventCount = length(spikeCounts{dataInd(event_i)}{chan_i}{unit_i}{1}.rates);
            eventNamePlot = plotLabels{event_i};
            eventName = eventsInEventDataPlot{dataInd(event_i)};

            % Open up the figure
            if figPerEvent
              % plotTitle = sprintf('%s, %s (p = %s, %d - %d ms, cD = %s, N = %d)', chUnitStr, eventName, pValStr, testPeriod(1), testPeriod(2), num2str(cohensD{chan_i}{unit_i}{dataInd(event_i)}, 3), eventCount);
              figTitle = sprintf('%s, %s (p = %s, %d - %d ms, cD = %s, N = %d)', chUnitStr, eventName, pValStr, testPeriod(1), testPeriod(2), num2str(cohensD{chan_i}{unit_i}{dataInd(event_i)}, 3), eventCount);
              figH = figure('Name', figTitle, 'NumberTitle','off','units','normalized', 'position', [0.3792    0.1528    0.3771    0.3213]);
              figTitleSave = sprintf('subEventPSTH_%s_%s', chUnitStr, eventName);
              plotTitle = sprintf('%s Aligned Activity', eventNamePlot);
            else
              plotTitle = sprintf('%s (p = %s, %d - %d ms, cD = %s, N = %d)', eventName, pValStr, testPeriod(1), testPeriod(2), num2str(cohensD{chan_i}{unit_i}{dataInd(event_i)}, 3), eventCount);
              axesH = subplot(rowCount, colCount, event_i);
            end
            
            % Variables for plotting
            plotData = [psthBySubEvent{chan_i}{unit_i}(dataInd(event_i), :); psthBySubEventNull{chan_i}{unit_i}(dataInd(event_i), :)];
            plotErr = [psthErrBySubEvent{chan_i}{unit_i}(dataInd(event_i), :); psthErrBySubEventNull{chan_i}{unit_i}(dataInd(event_i), :)];
            
            % Plot
            [~, ~, ~, legendH] = plotPSTH(plotData, plotErr, [], subEventParams.psthParams, 'line', plotTitle, [plotLabels(event_i), {'Null'}]);
            
            % Label
            ylabel('Fr (Hz)');
            xlabel('Time from sub-event onset (ms)');
            legendH.Location = 'southeast';
            hold on
            
            if dataInd(event_i) ~= length(subEventNames) && ~figPerEvent
              axesH.XLabel.String = '';
            end
            
            % If one figure per event, save
            if figPerEvent
              saveFigure(figStruct.figDir, figTitleSave, [], figStruct, figStruct.figTag)
            end
            
          end
          
          if ~figPerEvent
            saveFigure(figStruct.figDir, figTitleSave, [], figStruct, [])
          end
          
        end
      end
    end
    
    % SubEvent Occurances for stimulus visual events, not saccades/blinks.
%     if ~isempty(eventsInData)
%       figTitle = sprintf('Sub event occurances in Stimuli');
%       monkeyInd = contains(eventsInData, 'monkey');
%       resortInd = [find(monkeyInd); find(~monkeyInd)];
%       eventMat2Figure(stimEventMat(resortInd, :, :), eventsInData(resortInd), strrep(eventsInEventData, '_', ' '), figTitle)      
%       saveFigure(figStruct.figDir, figTitle, [], figStruct, [])
%     end
    
  end
  
  % Package outputs
  subEventSigStruct = struct();
  subEventSigStruct.subEventPSTH = [psthBySubEvent, psthErrBySubEvent];
  subEventSigStruct.subEventNullPSTH = [psthBySubEventNull, psthErrBySubEventNull];
  subEventSigStruct.testResults = testResults;
  subEventSigStruct.testPeriod = testPeriod;
  subEventSigStruct.testPeriod = subEventParams.testPeriodPerEvent;
  subEventSigStruct.testPeriodName = subEventParams.possibleEvents;
  subEventSigStruct.saccadeWindow = saccadeWindow;  
  subEventSigStruct.psthWindow = [-subEventParams.psthParams.psthPre subEventParams.psthParams.psthImDur];
  subEventSigStruct.cohensD = cohensD;
  subEventSigStruct.sigInd = "{channel}{unit}{event}";
  subEventSigStruct.events = subEventNames;
  subEventSigStruct.noSubEvent = 0;
  
end

tableVarNames = strcat("subSel_", subEventNames)';
tableVarNames = strrep(tableVarNames, '-', '_');
% In cases where adequate examples don't exist (seems like blinks) then NaN
% is in the selTable, make them 1.
selArray(isnan(selArray)) = 1;

for ev_i = 1:length(tableVarNames)
  selTable.([tableVarNames{ev_i} '_pVal']) = selArray(:, ev_i);
  selTable.([tableVarNames{ev_i} '_cohensD']) = diffArray(:, ev_i);
end

if ~taskData.rewardParadigm
  selTable.subSel_rewardAbsent_pVal = ones(size(selArray(:,1), 1), 1);
  selTable.subSel_rewardAbsent_cohensD = zeros(size(selArray(:,1), 1), 1);
end

specSubEventStruct = struct();

if specSubEvent && exist('uniqueSubEvents', 'var')
  specSubEventStruct.subEventOrigins = specSubEventNamesSorted;
  specSubEventStruct.subEventNames = uniqueSubEvents;
  specSubEventStruct.subEventCounts = onsetsByEventSpecCounts;
  specSubEventStruct.subEventPSTH = specSubEventPSTH;
end

end

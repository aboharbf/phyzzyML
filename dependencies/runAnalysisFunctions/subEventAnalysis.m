function [subEventSigStruct, specSubEventStruct, selTable] = subEventAnalysis(eyeBehStatsByStim, spikesByChannel, taskData, ephysParams, subEventParams, selTable, figStruct)
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

% Unpack Variables, some preprocessing. Generate an index for individual 
% trials that take place during the run, matching them to the stimuli as 
% present in the taskData.eventIDs field.
frameMotionData = taskData.frameMotionData;
eventIDs = subEventParams.eventIDs;
taskEventIDs = taskData.taskEventIDs;
specSubEvent = subEventParams.specSubEvent;

taskEventIDInd = zeros(length(taskEventIDs),1);
for ev_i = 1:length(eventIDs)
  taskEventIDInd = taskEventIDInd + (strcmp(taskEventIDs, eventIDs{ev_i}) * ev_i);
end
assert(max(taskEventIDInd) == length(eventIDs))

% Initialize a list of all the events present in the code.
subEventNames = [];
onsetsByEvent = {};
onsetsByEventNull = {};

% Event Type 1 - Find events which take place in the stimuli presented,
% defined in the eventData file.
eventDataFile = dir([subEventParams.stimDir '/**/eventData.mat']);
% If we find a file, load it.
load(fullfile(eventDataFile(1).folder, eventDataFile(1).name), 'eventData');
eventsInData = intersect(eventIDs, eventData.Row);

% If all the events in the run are in the eventData
if length(eventsInData) == length(eventIDs)
  eventDataRun = eventData(eventIDs, :);
  emptyEntryInd = cellfun(@(x) size(x,1) ~= 0, table2cell(eventDataRun));
  presentEventInd = any(emptyEntryInd);
  emptyEntryInd = emptyEntryInd(:,presentEventInd);
  eventsInEventData = eventDataRun.Properties.VariableNames(presentEventInd);
  eventStimTable = cell2table(cell(0,4), 'VariableNames', {'eventName', 'stimName', 'startFrame', 'endFrame'});
  
  % If our run has eventData specific events
  %   if ~isempty(eventsInEventData)
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
  stimPSTHDur = subEventParams.stimPlotParams.psthPre + subEventParams.stimPlotParams.psthImDur + subEventParams.stimPlotParams.psthPost;
  stimEventMat = initNestedCellArray(length(eventsInEventData), 'zeros', [length(eventIDs) stimPSTHDur]);
  
  for event_i = 1:length(eventsInEventData)
    % Identify spaces with the events in the presented stimuli.
    eventData = eventStimTable(strcmp(eventStimTable.eventName, eventsInEventData{event_i}),2:end);
    stimWithEvent = unique(eventData.stimName);
    
    % Identify trials of stimuli without this event, for null sampling
    [~, stimWOEventInd] = setdiff(eventIDs, stimWithEvent);
    [~, stimWEventInd] = intersect(eventIDs, stimWithEvent);
    taskEventNullSampling = false(size(taskEventIDInd));
    for sti_i = stimWOEventInd'
      taskEventNullSampling(taskEventIDInd == sti_i) = true;
    end
    
    [eventStartTimes, eventEndTimes, eventStartTimesNull, eventEndTimesNull, eventStimSource] = deal([]);
    for stim_i = stimWEventInd'
      % Get the appropriate start and end frames, convert to
      stimTrialInd = (taskEventIDInd == stim_i);
      stimEventData = eventData(strcmp(eventData.stimName, eventIDs{stim_i}), :);
      
      % Cycle through and find event occurance times
      for tab_i = 1:size(stimEventData,1)
        startTime = round(stimEventData.startFrame(tab_i));
        endTime = round(stimEventData.endFrame(tab_i));
        eventTitle = {[eventsInEventData{event_i}, '|', extractBefore(eventIDs{stim_i}, strfind(eventIDs{stim_i}, '.')) '-' num2str(startTime) '-' num2str(endTime)]};
        
        % Label the stimEventMat to highlight sampled region.
        startEventInd = (startTime + subEventParams.stimPlotParams.psthPre);
        endEventInd = (startEventInd + subEventParams.psthParams.psthImDur);
        stimEventMat{event_i}(stim_i, startEventInd:endEventInd) = deal(1);
        
        % Label the beginning to make lines easier to see stim.
        stimEventMat{event_i}(stim_i, 1:100) = deal(1);
        
        % Find out when the stimulus happens, and add these structures to
        % the larger ones
        eventStartTimes = [eventStartTimes; taskData.taskEventStartTimes(stimTrialInd) + startTime];
        eventEndTimes = [eventEndTimes; taskData.taskEventStartTimes(stimTrialInd) + endTime];
        eventStimSource = [eventStimSource; repmat(eventTitle, [length(taskData.taskEventStartTimes(stimTrialInd)), 1])];
        
        % Generate a matrix which can be used to determine a null
        % distribution
        %stimEventMat(stimMatIndex, startTime:endTime) = deal(1);
        nullTimesStart = taskData.taskEventStartTimes(taskEventNullSampling) + startTime;
        nullTimesEnd = taskData.taskEventStartTimes(taskEventNullSampling) + endTime;
        
        eventStartTimesNull = [eventStartTimesNull; nullTimesStart];
        eventEndTimesNull = [eventEndTimesNull; nullTimesEnd];
        
      end
      %eventStimSource = [eventStimSource; repmat(stimWithEvent(stim_i), [length(eventStartTimes), 1])];
    end
    onsetsByEvent{event_i} = eventStartTimes;
    onsetsByEventNull{event_i} = eventStartTimesNull;
    stimSourceByEvent{event_i} = eventStimSource;
  end
  
  % If we want, we can parse the events into individual instances, incase
  % there is some effect taking place for 'head turns during idle'
  stimSourceByEventAll = vertcat(stimSourceByEvent{:});
  uniqueSubEvents = unique(stimSourceByEventAll);
  if specSubEvent
    onsetsByEventAll = vertcat(onsetsByEvent{:});
    % Grab times of the individual events.
    onsetsByEventSpec = cell(length(uniqueSubEvents),1);
    onsetsByEventSpecCounts = zeros(length(uniqueSubEvents),1);
    for ii = 1:length(uniqueSubEvents)
      onsetsByEventSpec{ii} = onsetsByEventAll(strcmp(stimSourceByEventAll, uniqueSubEvents{ii}));
      onsetsByEventSpecCounts(ii) = length(onsetsByEventSpec{ii});
    end
    
    specEventInds = [length(onsetsByEvent) + 1; length(onsetsByEvent) + length(onsetsByEventSpec)];
    onsetsByEvent = [onsetsByEvent; onsetsByEventSpec];
    subEventNames = [subEventNames; uniqueSubEvents];
  end
  
else
  % If there are no events in the data
  subEventSigStruct = struct();
  subEventSigStruct.noSubEvent = 1;
end

% Event Type 2 - Generate/Organize Event times for eye movements here, concatonate onto
% the 'onsetsByEvent' cue.
if ~isempty(eyeBehStatsByStim)
  if isempty(subEventNames)
    subEventNames = {'blinks'};
  else
    subEventNames = [subEventNames; 'blinks'];
  end
  [blinkTimes, nullBlinkTimes] = deal([]);
  for stim_i = 1:length(eyeBehStatsByStim)
    % Find all the start times for the stimulus
    stimuliStartTimes = subEventParams.onsetsByEvent{stim_i};
    for trial_i = 1:length(eyeBehStatsByStim{stim_i})
      % for every trial, extract blink times, adding the appropriate event
      % Onset to get absolute eye event start times.
            
      if ~isempty(eyeBehStatsByStim{stim_i}{trial_i}.blinktimes)
        blinkTimes = [blinkTimes; eyeBehStatsByStim{stim_i}{trial_i}.blinktimes(1,:)'] + stimuliStartTimes(trial_i);
      end
      
    end
  end
  
  % Sort these lists for the next processing steps
  blinkTimes = sort(blinkTimes);
  eyeEventTimes = {blinkTimes};
  eyeEventNullTimes = {nullBlinkTimes};
  
  % Generate a distribution of null times for blinks
  firstTime = blinkTimes(1);
  lastTime = blinkTimes(end);

  % Generate Null times, add a number (100 - 1000), shuffle them and check 
  % if its close to a number on the list. If not, its a null value.
  maxShift = 500;
  minShift = 100;
  minEventDist = 100;
  
  for eyeEvent_i = 1:length(eyeEventTimes)
    eyeEventsTiled = eyeEventTimes{eyeEvent_i}; %repmat(eyeEventTimes{eyeEvent_i}, [shuffleMult,1]);
    shiftCount = length(eyeEventsTiled);
    shifts = [randi([minShift, maxShift],[ceil(shiftCount/2),1]); randi([-maxShift, -minShift],[floor(shiftCount/2),1])];
    shifts = shifts(randperm(length(shifts)));
    nullTimesForEvent = eyeEventsTiled + shifts; %repmat(eyeEventTimes{eyeEvent_i}, [shuffleMult,1]) + shifts;
    
    % check if any of these null times are close to event times, if not,
    % remove.
    nullTimesForEventCheck = repmat(nullTimesForEvent, [1, length(nullTimesForEvent)]);
    
    nullTimesForEventCheck = (nullTimesForEventCheck' - eyeEventsTiled)';
    nullTimesForEventCheckRemoveInd = abs(nullTimesForEventCheck) > minEventDist;
    keepInd = ~any(~nullTimesForEventCheckRemoveInd,2);
    nullTimesForEvent = nullTimesForEvent(keepInd);
    
    % Remove ones that are after the last time stamp or before the first
    keepInd = ~logical((nullTimesForEvent < firstTime) + (nullTimesForEvent > lastTime));
    nullTimesForEvent = nullTimesForEvent(keepInd);
    
    % Store into appropriate cell array
    eyeEventNullTimes{eyeEvent_i} = nullTimesForEvent;
  end
  
  onsetsByEvent = [onsetsByEvent; eyeEventTimes];
  onsetsByEventNull = [onsetsByEventNull; eyeEventNullTimes];
end

% Event Type 3 - add the Reward as an event
if subEventParams.RewardEvent
  
  % If this is one of the 'rewardParadigm' runs (where some trials are
  % intentially unrewarded to see if reward effects persist)
  
  
  % Generate a distribution of null times 
  rewardTimes = taskData.juiceOnTimes(~isnan(taskData.juiceOnTimes));
  firstTime = min(rewardTimes);
  lastTime = max(rewardTimes);
  
  % Generate Null times Take the Reward times, add a number (100 - 1000), shuffle them and check
  % if its close to a number on the list. If not, its a null value.
  maxShift = 500;
  minShift = 100;
  minEventDist = 200;
  
  shiftCount = length(rewardTimes);
  shifts = [randi([minShift, maxShift],[ceil(shiftCount/2),1]); randi([-maxShift, -minShift],[floor(shiftCount/2),1])];
  shifts = shifts(randperm(length(shifts)));
  nullTimesForEvent = rewardTimes + shifts; %repmat(eyeEventTimes{eyeEvent_i}, [shuffleMult,1]) + shifts;
  
  % check if any of these null times are close to event times, if not,
  % remove.
  nullTimesForEventCheck = repmat(nullTimesForEvent, [1, length(nullTimesForEvent)]);
  
  nullTimesForEventCheck = (nullTimesForEventCheck' - rewardTimes)';
  nullTimesForEventCheckRemoveInd = abs(nullTimesForEventCheck) > minEventDist;
  keepInd = ~any(~nullTimesForEventCheckRemoveInd,2);
  nullTimesForEvent = nullTimesForEvent(keepInd);
  
  % Remove ones that are after the last time stamp or before the first
  keepInd = ~logical((nullTimesForEvent < firstTime) + (nullTimesForEvent > lastTime));
  nullTimesForEvent = nullTimesForEvent(keepInd);
  
  % Store into appropriate cell array
  rewardTimesNull = nullTimesForEvent;
  
  onsetsByEvent = [onsetsByEvent; rewardTimes];
  onsetsByEventNull = [onsetsByEventNull; rewardTimesNull];
  
  % Stick onto larger structures
  if isempty(subEventNames)
    subEventNames = {'reward'};
  else
    subEventNames = [subEventNames; 'reward'];
  end
  
  % Generate a 2nd reward comparison, contraining null samples to
  % unrewarded trials.
  if taskData.rewardParadigm
    unRewardedTrials = isnan(taskData.juiceOnTimes);
    unRewardedStimStart = taskData.taskEventStartTimes(unRewardedTrials);
    rewardTime = subEventParams.stimPlotParams.psthImDur + 200;
    
    % When rewards would have taken place, normally.
    unRewardedStimRewardTimes = unRewardedStimStart + rewardTime;
    unRewardedStimRewardTimes = repmat(unRewardedStimRewardTimes, [10, 1]);
    
    onsetsByEvent = [onsetsByEvent; rewardTimes];
    onsetsByEventNull = [onsetsByEventNull; unRewardedStimRewardTimes];
    
    % Stick onto larger structures
    subEventNames = [subEventNames; 'rewardAbsent'];
    addFalseRowRewardAbsent = 0;
  else
    addFalseRowRewardAbsent = 1;
  end
 

end

%Ideally this would happen in some form after alignSpikes, instead of
%before, but the nested cell structure is more difficult to repmat then
%this, and may not be read the same way.

% Sample the null times, expand each onsertsByEventNull cell by the
% parameter below (subEventParams.nullSampleMult).
% onsetsByEventNull = cellfun(@(x) repmat(x, [subEventParams.nullSampleMult,1]), onsetsByEventNull, 'UniformOutput', false);
if ~isempty(onsetsByEvent)
  subEventParams.refOffset = 0;
  [spikesBySubEvent, spikesEmptyBySubEvent] = alignSpikes(spikesByChannel, onsetsByEvent, ones(length(spikesByChannel),1), subEventParams);
  [spikesBySubEventNull, spikesEmptyBySubEventNull] = alignSpikes(spikesByChannel, onsetsByEventNull, ones(length(spikesByChannel),1), subEventParams);
  
  % follow with putting spikesBySubEvent into calcPSTH.
  if ~ephysParams.spikeTimes
    spikesBySubEventBinned = calcSpikeTimes(spikesBySubEvent, subEventParams.psthParams);
    spikesBySubEventNullBinned = calcSpikeTimes(spikesBySubEventNull, subEventParams.psthParams);
    
    [psthBySubEvent, psthErrBySubEvent] = calcStimPSTH(spikesBySubEventBinned, spikesEmptyBySubEvent, ephysParams.spikeTimes, subEventParams.psthParams, subEventParams);
    [psthBySubEventNull, psthErrBySubEventNull] = calcStimPSTH(spikesBySubEventNullBinned, spikesEmptyBySubEventNull, ephysParams.spikeTimes, subEventParams.psthParams, subEventParams);
  else
    [psthBySubEvent, psthErrBySubEvent] = calcStimPSTH(spikesBySubEvent, spikesEmptyBySubEvent, ephysParams.spikeTimes, subEventParams.psthParams, subEventParams);
    [psthBySubEventNull, psthErrBySubEventNull] = calcStimPSTH(spikesBySubEventNull, spikesEmptyBySubEventNull, ephysParams.spikeTimes, subEventParams.psthParams, subEventParams);
  end
  
  if specSubEvent
    % Remove the subEvents which are specific instances, they need to be processed differently.
    subEventInd = false(length(subEventNames),1);
    % In case there are only saccades and blinks.
    if exist('specEventInds', 'var')
      subEventInd(specEventInds(1):specEventInds(2)) = true;
    end
    
    % Copy the structure generated
    specSubEventPSTH = psthBySubEvent;
    for chan_i = 1:length(specSubEventPSTH)
      for unit_i = 1:length(specSubEventPSTH{chan_i})
        % Go to the base layer, segregate info for specific subEvents and general ones.
        specSubEventPSTH{chan_i}{unit_i} = specSubEventPSTH{chan_i}{unit_i}(subEventInd,:);
        psthBySubEvent{chan_i}{unit_i} = psthBySubEvent{chan_i}{unit_i}(~subEventInd,:);
      end
    end
    
    % Modify other structures where event is at the top layer
    spikesBySubEvent = spikesBySubEvent(~subEventInd);
    subEventNames = subEventNames(~subEventInd);
    % uniqueSubEvents for specificSubEvents
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
    testPeriod = subEventParams.testPeriodPerEvent(strcmp(subEventParams.possibleEvents, subEventNames{event_i}), :);
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
        if subEventParams.nonParametric
          [pVal, ~, ~] = ranksum(eventSpikes,  nullSpikes);
        else
          [~, pVal, ~] = ttest2(eventSpikes,  nullSpikes);
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
    
    % Plotting Below
    tabPerEvent = 1; % Plot line PSTHes for each event.
    
    % Generate an array for referencing subEvents correctly.
    eventsInEventDataPlot = strrep(subEventNames, '_',' ');
    if exist('uniqueSubEvents', 'var') && specSubEvent
      tmp = cellfun(@(x) strsplit(x, '|'), uniqueSubEvents, 'UniformOutput', 0);
      specSubEventNamesSorted = cellfun(@(x) x{1}, tmp, 'UniformOutput', 0);
      %   eventsInEventDataPlot = strrep(subEventNames, '_',' ');
    end
    
    for chan_i = 1:length(psthBySubEvent)
      unitCount = length(psthBySubEvent{chan_i});
      for unit_i = 1:unitCount
        
        if unit_i == 1 && unitCount == 2
          continue % Only make plot for MUA.
        end
        
        % Create per unit figure
        if ~isempty(ephysParams.channelUnitNames{chan_i})
          psthTitle = sprintf('SubEvent comparison - %s %s', ephysParams.channelNames{chan_i}, ephysParams.channelUnitNames{chan_i}{unit_i});
          figure('Name',psthTitle,'NumberTitle','off','units','normalized', 'position', figStruct.figPos);
          if ~tabPerEvent
            sgtitle(psthTitle)
          end
          plotLabels = [{'Event'}, {'Null'}];
          
          % Cycle through events
          for event_i = 1:length(subEventNames)
            if ~tabPerEvent
              axesH = subplot(length(subEventNames), 1, event_i);
            else
              objTab = uitab('Title', subEventNames{event_i});
              axesH = axes('parent', objTab);
            end
            
            plotData = [psthBySubEvent{chan_i}{unit_i}(event_i, :); psthBySubEventNull{chan_i}{unit_i}(event_i, :)];
            plotErr = [psthErrBySubEvent{chan_i}{unit_i}(event_i, :); psthErrBySubEventNull{chan_i}{unit_i}(event_i, :)];
            plotTitle = sprintf('%s (p = %s, %d - %d ms, cohensD = %s, N = %d)', eventsInEventDataPlot{event_i}, num2str(testResults{chan_i}{unit_i}{event_i}), testPeriod(1), testPeriod(2), num2str(cohensD{chan_i}{unit_i}{event_i}), length(spikeCounts{event_i}{chan_i}{unit_i}{1}.rates));
            [psthHand, ~, vertLineHands] = plotPSTH(plotData, plotErr, [], subEventParams.psthParams, 'line', plotTitle , plotLabels);
            hold on
            if event_i ~= length(subEventNames)
              axesH.XLabel.String = '';
            end
            
            % Plot individual traces
            if specSubEvent
              relevantEventInd = find(strcmp(specSubEventNamesSorted, subEventNames{event_i}));
              newLines = gobjects(length(relevantEventInd),1);
              psthRange = -subEventParams.psthParams.psthPre:subEventParams.psthParams.psthImDur + subEventParams.psthParams.psthPost;
              for ii = 1:length(relevantEventInd)
                newLines(ii) = plot(psthRange, specSubEventPSTH{chan_i}{unit_i}(relevantEventInd(ii),:));
              end
              
              % Expand vertical bars to capture new traces
              [vertLineHands(1).YData, vertLineHands(2).YData] = deal(ylim());
              
              % Update the legend
              handles = [psthHand.mainLine, newLines'];
              labels = [plotLabels, uniqueSubEvents(relevantEventInd)'];
              legend(handles, labels);
            end
            
          end
          saveFigure(figStruct.figDir, psthTitle, [], figStruct, [])
        end
      end
    end
    
    % SubEvent Occurances for stimulus visual events, not saccades/blinks.
    if ~isempty(eventsInData)
      plotName = 'SubEvent Occurances';
      h = figure('Name', plotName,'NumberTitle','off','units','normalized');
      for event_i = 1:length(stimEventMat)
        axesH = subplot(length(stimEventMat), 1, event_i);
        [~, colorbarH] = plotPSTH(stimEventMat{event_i}, [], [], subEventParams.stimPlotParams, 'color', eventsInEventDataPlot(event_i) , strrep(eventIDs, '_', ' '));
        delete(colorbarH)
        axesH.FontSize = 10;
        if event_i ~= length(eventsInEventData)
          axesH.XLabel.String = '';
        end
      end
      saveFigure(figStruct.figDir, plotName, [], figStruct, [])
    end
    
  end
  
  % Package outputs
  subEventSigStruct = struct();
  subEventSigStruct.subEventPSTH = [psthBySubEvent, psthErrBySubEvent];
  subEventSigStruct.subEventNullPSTH = [psthBySubEventNull, psthErrBySubEventNull];
  subEventSigStruct.testResults = testResults;
  subEventSigStruct.testPeriod = testPeriod;
  subEventSigStruct.psthWindow = [-subEventParams.psthParams.psthPre subEventParams.psthParams.psthImDur];
  subEventSigStruct.cohensD = cohensD;
  subEventSigStruct.sigInd = "{channel}{unit}{event}";
  subEventSigStruct.events = subEventNames;
  subEventSigStruct.noSubEvent = 0;
  
end

tableVarNames = strcat("subSel_", subEventNames)';

% In cases where adequate examples don't exist (seems like blinks) then NaN
% is in the selTable, make them 1.
selArray(isnan(selArray)) = 1;

for ev_i = 1:length(tableVarNames)
  selTable.([tableVarNames{ev_i} '_pVal']) = selArray(:, ev_i);
  selTable.([tableVarNames{ev_i} '_cohensD']) = diffArray(:, ev_i);
end

if addFalseRowRewardAbsent
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

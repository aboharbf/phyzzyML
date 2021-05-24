function [subEventData, dataType2] = compileSubEventPSTHArray(stimPSTH, dataType, eventData, eventLists, groupIterInd, psthPreData, params)
% Function which iterates across the stimulus specific stimPSTH array,
% extracts subEventData, and compiles it into an array
% Output: subEventData{event}{group_i}{dataType2}

stimEventCountarray = cellfun(@(x) size(x, 1), table2cell(eventData));
stimEventLogicArray = stimEventCountarray ~= 0;
stimEvent2Plot = any(stimEventLogicArray,2);        % Only worry about plotting stimuli below which have events.

traceLength = size(stimPSTH{1}, 2);
psthPre = params.subEventPSTHParams.psthPre;
psthImDur = params.subEventPSTHParams.psthImDur;
psthPost = params.subEventPSTHParams.psthPost;

dataType2 = {'subEventPSTH', 'Label', 'RunInd'};

% Generate subEvent Data structure
subEventData = cell(length(eventLists), 3, length(dataType2)); % I took out Dim 2, which is usually group. Everything will be set to MUA for now.
stimuli = eventData.Row;
eventList = eventData.Properties.VariableNames';
eventInstIndex = ones(size(eventList));

eventMaxDur = zeros(size(eventLists));
for ii = 1:length(eventLists)
  tmp = cellfun(@(x) max(x.endTime - x.startTime), eventData{:, ii}, 'UniformOutput', false);
  maxLen = max([tmp{:}]);
  if ~isempty(maxLen)
    eventMaxDur(ii) = ceil(max([tmp{:}]));
  end
end

for stim_i = 1:length(stimuli)
  if stimEvent2Plot(stim_i)
    for group_i = groupIterInd
      for event_i = 1:length(eventList)
        
        % Gather info for all the instances of a specific event in a
        % stimulus, and loop across them.
        evTable = eventData{stim_i, eventList{event_i}}{1};
        evStarts = round(evTable.startTime);
        evEnds = round(evTable.endTime);
        eventLengths = evEnds - evStarts;
        
        for inst_i = 1:length(eventLengths)
          % Expand the desired stop and start by the amount outside of
          % the event you'd like to see.
          startT = evStarts(inst_i) + psthPreData - psthPre;
          endT = min(evEnds(inst_i) + psthPreData + psthImDur + psthPost, traceLength);            % If endT exceeds the length of the trace, then adjust it.
          endDist = min(traceLength - evEnds(inst_i), psthPost);
          
          % Gather Relevant data
          subEventPSTH = stimPSTH{stim_i, group_i, strcmp(dataType, 'PSTH')}(:, startT:endT);
          subEventRunInds = stimPSTH{stim_i, group_i, strcmp(dataType, 'Run Ind')};
          subEventLabel = sprintf('%s, %d - %d', stimuli{stim_i}, evStarts(inst_i), evEnds(inst_i));
          
          % Retrieve and increment the storage index.
          storInd = eventInstIndex(event_i);
          eventInstIndex(event_i) = eventInstIndex(event_i) + 1;
          
          % subEventPSTH has to be padded to allow for stacking with
          % different events later. This Padding involves removing
          % post event activity.
          subEventPSTHPadded = nan(size(subEventPSTH,1), eventMaxDur(event_i));
          endIndForData = size(subEventPSTH, 2) - endDist;
          subEventPSTHPadded(:, 1:endIndForData) = subEventPSTH(:, 1:end - endDist); % Using endDist here to avoid possibly cutting out activity during event, when postEvent time can't be 200 ms.
          
          % Store relevant data
          subEventData{event_i, group_i, strcmp(dataType2, 'subEventPSTH')}{storInd} = subEventPSTHPadded;
          subEventData{event_i, group_i, strcmp(dataType2, 'Label')}{storInd} = subEventLabel;
          subEventData{event_i, group_i, strcmp(dataType2, 'RunInd')}{storInd} = subEventRunInds;
        end
      end
    end
  end
end
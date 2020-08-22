function eventMat = generateEventImage(eventData, stimDur)
% Generates an event image which can be shown next to an image to get a
% sense for correlations between activity and events.

% Get names
stimNames = eventData.Properties.RowNames;
eventNames = eventData.Properties.VariableNames;

% Initialize
eventMat = zeros(length(stimNames), stimDur, length(eventNames));

% Loop through and fill out matrix
for event_i = 1:length(eventNames)
  eventDataSingle = eventData(:,event_i);
  
  for stim_i = 1:length(stimNames)
    eventTable = eventDataSingle{stim_i,:}{1};
    if ~isempty(eventTable)
      starts = eventTable.startTime;
      stops = eventTable.endTime;
      
      % Fill in the eventMat for each instance of the event
      for ev_i = 1:length(starts)
        eventInd = round(starts(ev_i)):round(stops);
        eventMat(stim_i, eventInd, event_i) = 1;
      end
    end    
  end
  
end

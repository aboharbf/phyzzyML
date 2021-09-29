function dataNew = deEventSpikeData(ByEventData)
% A function which cycles through event based data and collapses across
% them.
% Inputs:
% - byEventData: some data with the structure {event_i}{chan_i}{unit_i}
% Outputs:
% - data: data with structure {chan_i}{unit_i}

totalTrials = sum(cellfun(@(x) size(x{1}{1}, 1), ByEventData));
binCount = size(ByEventData{1}{1}{1}, 2);
dataNew = initNestedCellArray(ByEventData{1}, 'cell', [0 0], 100);

for chan_i = 1:length(dataNew)
  for unit_i = 1:length(dataNew{chan_i})
    % Initialize
    unitDataAllTrials = nan(totalTrials, binCount);
    emptySlot = 1;
    
    % Cycle through event, concatonate data
    for event_i = 1:length(ByEventData)
      d2Add = ByEventData{event_i}{chan_i}{unit_i};
      endSpot = emptySlot+size(d2Add,1)-1;
      unitDataAllTrials(emptySlot:endSpot, :) = d2Add;
      emptySlot = endSpot+1;
    end
    
    dataNew{chan_i}{unit_i} = unitDataAllTrials;
    
  end
end

end
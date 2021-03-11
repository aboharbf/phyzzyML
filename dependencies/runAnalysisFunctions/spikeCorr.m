function spikeBinCorr = spikeCorr(spikesByEventBinned)
%functions performs a two-way ANOVA, comparing the firing rates for each
%Epoch for all members of "target" and all members of the rest of the
%groups. Returns the ANOVA table for each such comparison.

%spikesByEvent{event}{channel}{unit}(trial).times = spikes*1
%spikesByEventBinned{event}{channel}{unit} = trials*bins;

%Cycle through structure, concatonating the correct events.
for event_i = 1:length(spikesByEventBinned)
  for channel_i = 1:length(spikesByEventBinned{event_i})
    for unit_i = 1:length(spikesByEventBinned{event_i}{channel_i})
      spikeBins = spikesByEventBinned{event_i}{channel_i}{unit_i};
      spikeBinCorr{event_i}{channel_i}{unit_i} = corrcoef(spikeBins');
    end
  end
end
end

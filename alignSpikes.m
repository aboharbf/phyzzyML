function [ spikesByItem, psthEmptyByItem ] = alignSpikes( spikesByChannel, alignPointsByItem, spikeChannels, params )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if isempty(spikesByChannel)
  spikesByItem = {};
  psthEmptyByItem = {};
  return
end

preAlign = params.preAlign;
postAlign = params.postAlign;
refOffset = params.refOffset; % this value is subtracted from all times; so, to have all positive times, should be -msPreAlign
itemCount = length(alignPointsByItem);
chanCount = length(spikeChannels);

[spikesByItem, psthEmptyByItem] = deal(cell(itemCount, 1));
for item_i = 1:itemCount
  [itemSpikesByChannel, itemEmptyByChannel] = deal(cell(chanCount,1)); %channels x 1
  onsets = alignPointsByItem{item_i}; 
  for channel_i = 1:chanCount
    channelSpikes = spikesByChannel(channel_i);
    units = sort(unique(channelSpikes.units));
    [itemChannelSpikes, itemChannelEmpty] = deal(cell(length(units)+1,1));

    for unit_i = 1:length(units)
      tstamps = channelSpikes.times(channelSpikes.units == units(unit_i));
      itemUnitSpikes = repmat(struct('times',[]),length(onsets),1); 
      empty = 1;
      for trial_i = 1:length(onsets)
        itemUnitSpikes(trial_i) = struct('times', tstamps(tstamps>=onsets(trial_i) - preAlign & tstamps <= onsets(trial_i)+postAlign) - onsets(trial_i) - refOffset);
        if ~isempty(itemUnitSpikes(trial_i).times)
          empty = 0;
        end
      end
      itemChannelSpikes{unit_i} = itemUnitSpikes;
      itemChannelEmpty{unit_i} = empty;
    end
    
    % channel MUA
    tstamps = channelSpikes.times;
    itemChannelMUA = repmat(struct('times',[]),length(onsets),1);
    empty = 1;
    for trial_i = 1:length(onsets)
      itemChannelMUA(trial_i) = struct('times', tstamps(tstamps>=onsets(trial_i) - preAlign & tstamps <= onsets(trial_i)+postAlign) - onsets(trial_i) - refOffset);
      if ~isempty(itemChannelMUA(trial_i).times)
        empty = 0;
      end
    end
    
    itemChannelSpikes{end} = itemChannelMUA;
    itemChannelEmpty{end} = empty;
    itemSpikesByChannel{channel_i} = itemChannelSpikes;
    itemEmptyByChannel{channel_i} = itemChannelEmpty;
  end
  spikesByItem{item_i} = itemSpikesByChannel;
  psthEmptyByItem{item_i} = itemEmptyByChannel;
end
end


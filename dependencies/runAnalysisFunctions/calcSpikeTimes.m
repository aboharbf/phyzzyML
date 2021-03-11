function spikesByStimBinned = calcSpikeTimes(spikesByStim, psthParams)
%Wrapper which cycles through a particular Stimulus tier (Catagory or
%Event/Image) and feeds it into the chronux binspikes function with bin
%sizes of 1 ms. spikesByStim must be {event}{channel}{unit}

% Unpack Variables
psthPre = psthParams.psthPre;
psthImDur = psthParams.psthImDur;
psthPost = psthParams.psthPost;
movingWin = psthParams.movingWin;
smoothingWidth = psthParams.smoothingWidth;

assert((movingWin(1)/2 >= 3*smoothingWidth),'Error: current implementation assumes that movingWin/2 > 3*psthSmoothingWidth. Not true here');

spikesByStimBinned = cell(size(spikesByStim));
for image_i = 1:length(spikesByStim)
  spikesByStimBinned{image_i} = cell(length(spikesByStim{1}),1);
  for channel_i = 1:length(spikesByStim{1})
    spikesByStimBinned{image_i}{channel_i} = cell(length(spikesByStim{1}{channel_i}),1);
    for unit_i = 1:length(spikesByStim{1}{channel_i})
      % this tmp variable fixes a Chronux bug that leads the last entry in the binned spike array to be zero when it shouldn't be
      spikesTmp = binspikes(spikesByStim{image_i}{channel_i}{unit_i},1,[-1*(psthPre+movingWin(1)/2), psthImDur+1+psthPost+movingWin(1)/2])';
      spikesByStimBinned{image_i}{channel_i}{unit_i} = spikesTmp(:,1:end-1);
      %         for trial_i = 1:size(spikesByEventBinned{image_i}{channel_i}{unit_i},1)
      %           fieldName = sprintf('spikesBinned.%s_%s',channelNames{channel_i}(~isspace(channelNames{channel_i})),channelUnitNames{channel_i}{unit_i}(~isspace(channelUnitNames{channel_i}{unit_i})));
      %           trialDB = trialDatabaseSetField(fieldName,spikesByEventBinned{image_i}{channel_i}{unit_i}(trial_i,:),trialDB,trialIDsByEvent{image_i}(trial_i),'dateSubj',dateSubject,'runNum',runNum);
      %         end
    end
  end
end
end

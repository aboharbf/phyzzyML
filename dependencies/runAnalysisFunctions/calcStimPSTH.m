function [psthByStim, psthErrByStim] = calcStimPSTH(spikesByStim, psthEmptyByStim, spikeTimes, psthParams, spikeAlignParams)

%spikesByStim = some data structure input (spikesByEvent/Catagory,
%optionally binned).
preAlign = spikeAlignParams.preAlign;
postAlign = spikeAlignParams.postAlign;

smoothingWidth = psthParams.smoothingWidth;
movingWin = psthParams.movingWin;
times = -psthParams.psthPre:(psthParams.psthImDur+psthParams.psthPost);
psthErrorType = psthParams.errorType;
psthErrorRangeZ = psthParams.errorRangeZ;
psthBootstrapSamples = psthParams.bootstrapSamples;

spikesByItem = spikesByStim;
psthEmptyByItem = psthEmptyByStim;
numItems = length(spikesByStim);

if ~spikeTimes
  filterPoints = -3*smoothingWidth:3*smoothingWidth;
  smoothingFilter = exp(-1*filterPoints.^2/(2*smoothingWidth^2));
  smoothingFilter = smoothingFilter/sum(smoothingFilter);
end

% Initialize vectors for PSTH and accompanying error.
[psthByItem, psthErrByItem] = deal(cell(length(spikesByStim{1}),1));

% Cycle through channels, units, then items.
for channel_i = 1:length(spikesByStim{1})
  [channelItemPsth, channelItemPsthErr] = deal(cell(length(spikesByStim{1}{channel_i}),1));
  for unit_i = 1:length(spikesByStim{1}{channel_i})
    [unitItemPsth, unitItemPsthErr] = deal(zeros(numItems,length(times)));
    for item_i = 1:numItems
      spikeData = spikesByItem{item_i}{channel_i}{unit_i};
      if ~psthEmptyByItem{item_i}{channel_i}{unit_i}
        if spikeTimes
          [paddedPsth,~,paddedPsthErr] = psth(spikeData,smoothingWidth, 'n', [-preAlign postAlign], min(psthErrorType,2), -preAlign:postAlign);
          paddedPsth = 1000*paddedPsth;
          paddedPsthErr = 1000*paddedPsthErr;
          unitItemPsth(item_i,:) = paddedPsth(3*smoothingWidth+1:end-3*smoothingWidth);
          unitItemPsthErr(item_i,:) = paddedPsthErr(3*smoothingWidth+1:end-3*smoothingWidth)*psthErrorRangeZ/2; %note: psth returns +/- 2 stderr; thus the factor of 0.5
        else % use spike bins
          paddedPsth = 1000*conv(mean(spikeData, 1),smoothingFilter,'same');
          % chronux convention: 1 is poisson, 2 is trialwise bootstrap, 3 is across trial std for binned spikes, bootstrap for spike times 
          if psthErrorType == 1 || size(spikeData, 1) == 1  % if only one trial, can't use bootstrap
            %note: the factor of sqrt(1000) appears because we need to convert to spks/ms to get poisson error, then back to Hz
            paddedPsthErr = sqrt(paddedPsth)*(sqrt(1000)*psthErrorRangeZ/sqrt(size(spikeData, 1)));
          elseif psthErrorType == 2
            opt = statset('UseParallel',true);
            paddedPsthErr = std(bootstrp(psthBootstrapSamples, @(x) mean(x,1), convn(spikeData, smoothingFilter, 'same'), 'Options', opt), [], 1) * psthErrorRangeZ * 1000;
          elseif psthErrorType == 3
            paddedPsthErr = std(convn(spikeData, smoothingFilter,'same'),[],1) * psthErrorRangeZ * 1000/sqrt(size(spikeData,1));
          else
            error('psthErrorType parameter, when using spike bins, must be between 1 and 3')
          end
          unitItemPsth(item_i,:) = paddedPsth(movingWin(1)/2+1:end-movingWin(1)/2);
          unitItemPsthErr(item_i,:) = paddedPsthErr(movingWin(1)/2+1:end-movingWin(1)/2)*psthErrorRangeZ;
        end
      end
    end
    channelItemPsth{unit_i} = unitItemPsth;
    channelItemPsthErr{unit_i} = unitItemPsthErr;
  end
  psthByItem{channel_i} = channelItemPsth;
  psthErrByItem{channel_i} = channelItemPsthErr;
end

%Properly package outputs with the right names, save them.
psthByStim = psthByItem;
psthErrByStim = psthErrByItem;
end

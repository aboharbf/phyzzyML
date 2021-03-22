function [sigStruct, imageSortOrder, nullModelPvalues, nullTraceMeans, nullTraceSD] = nullModelStats(psthByImage, spikeCountsByImageByEpoch, firingRatesByImageByEpoch, firingRateErrsByImageByEpoch, trialCountsByImage, analysisGroups, epochLabels, eventIDs, ephysParams)

% Outputs
% sigStruct.

groups = analysisGroups.stimulusLabelGroups.groups{1};
channelNames = ephysParams.channelNames;

[imageSortOrder, imageSortedRates] = deal(initNestedCellArray(spikeCountsByImageByEpoch, 'cell', [0, 0], 3));
[nullModelPvalues, nullTraceMeans, nullTraceSD] = deal(initNestedCellArray(spikeCountsByImageByEpoch, 'cell', [length(groups), 1], 3));

% Method 1 - Using Serene's null model of spiking, binned across entire
% stimulus presentation.
for epoch_i = 1:length(spikeCountsByImageByEpoch)
  for channel_i = 1:length(spikeCountsByImageByEpoch{epoch_i})
    for unit_i = 1:size(spikeCountsByImageByEpoch{epoch_i}{channel_i},1)
      
      % Sort the firingRatesByImageByEpoch for the sake of the null model
      [imageSortedUnit, imageSortOrderUnit] = sort(firingRatesByImageByEpoch{epoch_i}{channel_i}(unit_i,:), 2, 'descend');  %todo: write the firing rates to file
      imageSortedRates{epoch_i}{channel_i}{unit_i} = imageSortedUnit;
      imageSortOrder{epoch_i}{channel_i}{unit_i} = imageSortOrderUnit;
%       spikeCountsUnitSorted = spikeCountsByImageByEpoch{epoch_i}{channel_i}{unit_i}(imageSortOrderUnit);
      imFrErrSorted = firingRateErrsByImageByEpoch{epoch_i}{channel_i}(unit_i,imageSortOrderUnit);
      %       sortedImageLabels = eventLabels(imageSortOrder{epoch_i}{channel_i}{unit_i});
      %       sortedEventIDs = eventIDs(imageSortOrder{epoch_i}{channel_i}{unit_i});
      trialCountsByImageSorted = trialCountsByImage(imageSortOrderUnit);
      
      %Calculate the Null model
      mu = mean(imageSortedUnit);
      if any(trialCountsByImageSorted > 3)
        rawSigmas = imFrErrSorted(trialCountsByImageSorted > 1).*sqrt(trialCountsByImageSorted(trialCountsByImageSorted > 1))';
        mixSEs = horzcat(imFrErrSorted(trialCountsByImageSorted > 1),randsample(rawSigmas,sum(trialCountsByImageSorted == 1),true));
      else
        mixSEs = zeros(1, size(trialCountsByImageSorted, 1));
      end
      
      nulltrials = 100;
      nullDistSortsMix = zeros(nulltrials,length(imageSortedUnit));
      for i = 1:length(imageSortedUnit)
        nullDistSortsMix(:,i) = normrnd(mu,mixSEs(i),nulltrials,1);
      end
      nullDistSortsMix = sort(nullDistSortsMix,2,'descend');
      nullDistSortsMix(nullDistSortsMix < 0) = 0;
      muNull = mean(nullDistSortsMix, 1);
      nullTraceSD{epoch_i}{channel_i}{unit_i} = std(nullDistSortsMix, 1);
      nullTraceMeans{epoch_i}{channel_i}{unit_i} = muNull;
      
      %Calculate significance
      zScore = (imageSortedUnit - muNull)./imFrErrSorted;
      nullModelPvalues{epoch_i}{channel_i}{unit_i} = exp(-0.717.*zScore-0.416*(zScore.^2));
      
      % Different Approach - generate null distribution, calculate t test
      % for each spot.
      % Check how error is calculated in spikeCounter - (std(s.rates)/sqrt(n))*sqrt(2/(n-1))*(gamma(n/2)/gamma((n-1)/2))
      
    end
  end
end

% After generate p values with the null model, generate the sigStruct.
% sigStruct.channels = ephysParams.spikeChannels;
groupingType = {'Unsorted','Unit','MUA'};
dataType = {'PSTH', 'Label','Stimuli'};
sigStruct.IndInfo = {epochLabels, groupingType, dataType};
sigStruct.channelNames = channelNames;
sigStruct.unitCount = zeros(length(channelNames),1);
sigStruct.sigInfo = zeros(length(channelNames), length(groupingType));

for channel_i = 1:length(psthByImage)
  data = cell([length(spikeCountsByImageByEpoch), length(groupingType), length(dataType)]);
  sigStruct.unitCount(channel_i) = length(psthByImage{channel_i}) - 2;
  for unit_i = 1:length(psthByImage{channel_i})
    unitPSTH = psthByImage{channel_i}{unit_i};
    for epoch_i = 1:size(data,1)
      % Find out which stimuli produced the significant activity
      runMask = (nullModelPvalues{epoch_i}{channel_i}{unit_i} < 0.05);
      
      if any(runMask)
        sortVec = imageSortOrder{epoch_i}{channel_i}{unit_i};
        stimOfInterest = sortVec(runMask);
        PSTHtoPlot = unitPSTH(stimOfInterest,:);       % Recover vectors of interest
        stimtoPlot = eventIDs(stimOfInterest);         % Sort the eventList
        
        % Identify the unit type
        if unit_i == 1
          groupingInd = 1;
          unitLabel = groupingType{groupingInd};
        elseif unit_i == length(psthByImage{channel_i})
          groupingInd = 3;
          unitLabel = groupingType{groupingInd};
        else
          groupingInd = 2;
          unitLabel = ['U' num2str(unit_i - 1)];
        end
        sigStruct.sigInfo(channel_i, groupingInd) =  sigStruct.sigInfo(channel_i, groupingInd) + 1;
        %{sprintf('Ch%d U%d',[channel_i, unit_i+1])}
        
        dataArray = cell(length(dataType),1);
        dataArray{1} = PSTHtoPlot;
        dataArray{2} = repmat({sprintf('Ch%d %s',[channel_i, unitLabel])}, [3,1]);
        dataArray{3} = stimtoPlot;
        
        for data_i = 1:length(dataArray)
          data{epoch_i, groupingInd, data_i} = [data{epoch_i, groupingInd, data_i}; dataArray{data_i}];
        end
        
      end
    end
  end
  sigStruct.data{channel_i} = data;
end

end

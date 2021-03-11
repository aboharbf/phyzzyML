function [sigStruct, imageSortOrder, nullModelPvalues, nullTraceMeans, nullTraceSD, frEpochsStats] = genStats(psthByImage, spikeCountsByImageByEpoch, firingRatesByImageByEpoch, firingRateErrsByImageByEpoch, trialCountsByImage, analysisGroups, epochLabels, eventIDs, ephysParams, genStatsParams)

groups = analysisGroups.stimulusLabelGroups.groups{1};
channelNames = ephysParams.channelNames;

[imageSortOrder, imageSortedRates] = deal(initNestedCellArray(spikeCountsByImageByEpoch, 'cell', [0, 0], 3));
[nullModelPvalues, nullTraceMeans, nullTraceSD] = deal(initNestedCellArray(spikeCountsByImageByEpoch, 'cell', [length(groups), 1], 3));

for epoch_i = 1:length(spikeCountsByImageByEpoch)
  for channel_i = 1:length(spikeCountsByImageByEpoch{epoch_i})
    for unit_i = 1:size(spikeCountsByImageByEpoch{epoch_i}{channel_i},1)
      % Sort the firingRatesByImageByEpoch
      [imageSortedUnit, imageSortOrderUnit] = sort(firingRatesByImageByEpoch{epoch_i}{channel_i}(unit_i,:),2,'descend');  %todo: write the firing rates to file
      imageSortedRates{epoch_i}{channel_i}{unit_i} = imageSortedUnit;
      imageSortOrder{epoch_i}{channel_i}{unit_i} = imageSortOrderUnit;
      spikeCountsUnitSorted = spikeCountsByImageByEpoch{epoch_i}{channel_i}{unit_i}(imageSortOrderUnit);
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
      runMask = (nullModelPvalues{epoch_i}{channel_i}{unit_i} < 0.05);   % Find out which stimuli produced the significant activity
      if sum(runMask) ~= 0
        sortVec = imageSortOrder{epoch_i}{channel_i}{unit_i};
        stimOfInterest = sortVec(runMask);
        PSTHtoPlot = unitPSTH(stimOfInterest,:);       % Recover vectors of interest
        stimtoPlot = eventIDs(stimOfInterest);         % Sort the eventList
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

% Linear Model Code
% ANOVA - functions performs a two-way ANOVA, comparing the firing rates for each
%Epoch for all members of "target" and all members of the rest of the
%groups. Returns the ANOVA table for each such comparison.
group = genStatsParams.ANOVAParams.group;
target = genStatsParams.ANOVAParams.target;
groupLabelsByImage = genStatsParams.ANOVAParams.groupLabelsByImage;

%a phase 2 trial w/ only social stuff, scrambles w/ only scrambles.
if length(strcmp(group,target)) == 1 || ~(length(unique(groupLabelsByImage)) > 1)
  target = 0;
end

%spikeCountsByImageByEpoch{epoch}{channel}{unit}{event}.rates = trials*1
%Will likely need changes in the future with respect to group variable
%will lead to errors on cases with more than one group.
%Cycle through structure, concatonating the correct events.
frEpochsStats = initNestedCellArray(spikeCountsByImageByEpoch{1}, 'struct',[0,0],2);

for channel_i = 1:length(spikeCountsByImageByEpoch{1})
  for unit_i = 1:length(spikeCountsByImageByEpoch{1}{channel_i})
    [trialSpikes, trialLabels, trialEpoch]  = deal([]);
    for epoch_i = 1:length(spikeCountsByImageByEpoch)
      unitResponsePerEvent = spikeCountsByImageByEpoch{epoch_i}{channel_i}{unit_i};
      %Target behaves as a switch. If there is a target, it becomes Target
      %v All ANOVA, if not, each eventID is considered.
      if target
        ANOVAvarName = ['Ch' num2str(channel_i) 'U' num2str(unit_i) ' - SocVsNonSoc Label'];
        %grab the relevant events
        targetInd = groupLabelsByImage == find(strcmp(group,target));
        targetSpikes = unitResponsePerEvent(targetInd);
        otherSpikes = unitResponsePerEvent(~targetInd);
        %Initialize relevant vecotrs
        spikeGroups = {targetSpikes otherSpikes};
        spikeGroupLabels ={(target) (['non-' target])};
        %Cluster and reshape the arrays properly
        for group_i = 1:length(spikeGroups)
          tmp = spikeGroups{group_i};
          tmp = [tmp{:}];
          dataVec = vertcat(tmp.rates);
          labelVec = repmat(spikeGroupLabels(group_i), length(dataVec),1);
          epochVec = repmat(epochLabels(epoch_i), length(dataVec),1);
          trialSpikes = vertcat(trialSpikes,dataVec);
          trialLabels = vertcat(trialLabels, labelVec);
          trialEpoch = vertcat(trialEpoch, epochVec);
        end
      else
        ANOVAvarName = ['Ch' num2str(channel_i) 'U' num2str(unit_i) ' - Event Labels'];
        spikes = unitResponsePerEvent;
        spikeLabels = eventIDs;
        for group_i = 1:length(spikes)
          trialSpikes = vertcat(trialSpikes,spikes{group_i}.rates);
          trialLabels = vertcat(trialLabels, repmat(spikeLabels(group_i),length(spikes{group_i}.rates),1));
          trialEpoch = vertcat(trialEpoch, repmat(epochLabels(epoch_i),length(spikes{group_i}.rates),1));
        end
      end
    end
    % Check for task modulation
    [frEpochsStats{channel_i}{unit_i}.taskModulatedP, ~, ~] = anova1(trialSpikes, trialEpoch, 'off');%,'model','interaction','varnames',{'Epoch'}, 'alpha', 0.05,'display','off');
    [pVal, cohensD] = deal(cell(length(epochLabels),1));
    % Perform t test on each epoch
    for epoch_i = 1:length(epochLabels)
      socSpikesPres = trialSpikes((strcmp(trialEpoch,epochLabels{epoch_i}) + strcmp(trialLabels, target)) == 2);
      nonSocSpikesPres = trialSpikes((strcmp(trialEpoch,epochLabels{epoch_i}) + ~strcmp(trialLabels, target)) == 2);
      %[Presp, ~, Pstats] = anovan(trialSpikesPres,{trialLabelsPres},'model','interaction','varnames',{ANOVAvarName}, 'alpha', 0.05,'display','off');
      %[Presp, ~, Pstats] = anova1(trialSpikesPres, trialLabelsPres, 'off'); %'model','interaction','varnames',{ANOVAvarName}, 'alpha', 0.05,'display','off');
      [~, pVal{epoch_i}, ~, tmpStats] = ttest2(socSpikesPres,nonSocSpikesPres);
      cohensD{epoch_i} = (mean(socSpikesPres) - mean(nonSocSpikesPres))/tmpStats.sd;
    end
    % Store in struct
    frEpochsStats{channel_i}{unit_i}.tTest.epochLabels = epochLabels;
    frEpochsStats{channel_i}{unit_i}.tTest.target = target;
    frEpochsStats{channel_i}{unit_i}.tTest.pVals = pVal;
    frEpochsStats{channel_i}{unit_i}.tTest.cohensD = cohensD;
  end
end

end

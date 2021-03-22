function [frEpochsStats] = anovaStats(spikesByEvent, epochLabels, eventIDs, ephysParams, ANOVAParams)
% This function performs a few statistical tests to determine whether there
% are meaningful difference between different periods of a task, and
% different conditions of the task.

% First, gather spike counts in the respective 
timeBins = ANOVAParams.times;
labels = ANOVAParams.labels;
for epoch_i = 1:size(timeBins,1)
  [spikeCounts, ~, ~] = spikeCounter(spikesByEvent, timeBins(epoch_i, 1), timeBins(epoch_i, 2));
  spikeCountsByImageByEpoch{epoch_i} = spikeCounts;
%   firingRatesByImageByEpoch{epoch_i} = fr;
%   firingRateErrsByImageByEpoch{epoch_i} = frErr;
end

% Linear Model Code
% ANOVA - functions performs a two-way ANOVA, comparing the firing rates for each
% Epoch for all members of "target" and all members of the rest of the
% groups. Returns the ANOVA table for each such comparison.
group = ANOVAParams.group;
target = ANOVAParams.target;
groupLabelsByImage = ANOVAParams.groupLabelsByImage;

%a phase 2 trial w/ only social stuff, scrambles w/ only scrambles.
if length(strcmp(group,target)) == 1 || ~(length(unique(groupLabelsByImage)) > 1)
  target = 0;
end

%spikeCountsByImageByEpoch{epoch}{channel}{unit}{event}.rates = trials*1

%Will likely need changes in the future with respect to group variable
%will lead to errors on cases with more than one group.
%Cycle through structure, concatonating the correct events.
% frEpochsStats = initNestedCellArray(spikeCountsByImageByEpoch{1}, 'struct',[0,0],2);

% NEW STUFF
% Use 'signrank' - https://www.mathworks.com/help/stats/signrank.html
% between rates for baseline and rates in period of interest.
% Use 'ranksum' - https://www.mathworks.com/help/stats/ranksum.html between
% rates for social and non-social conditions.

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

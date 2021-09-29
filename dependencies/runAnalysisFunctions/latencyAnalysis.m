function selTable = latencyAnalysis(selTable, spikesByEventBinned, spikesByEventBinnedFixAlign, psthParams, taskData)

% Determine the latency for each unit. Do the analysis based on the
% apperance of the fixation cue or the appearance of the stimulus.
% Inputs:
% - selTable: Where the outcome will be saved.
% - spikesByEventBinned: binned spikes aligned to the stimulus onset.
% - spikesByEventBinned: binned spikes aligned to the fixation point.
% Outputs:
% - selTable, row 'latencty' - an entry for each unit describing its
% latency in ms.
% prefixWin = spikeAlignParams.preAlign-500:spikeAlignParams.preAlign;
% fixWin = spikeAlignParams.preAlign:spikeAlignParams.preAlign+800;

preStimTime = psthParams.psthPre + psthParams.movingWin(1)/2;
preStimTime = psthParams.psthPre + psthParams.movingWin(1)/2;

% Determine the fixation duration for each trial
taskEventList = taskData.taskEventList;
latestReward = round(max(taskData.rewardTimePerTrial));
rewardEndInd = latestReward + preStimTime+ 50;


% fixDirSorted = [];
fixDirSorted = cell(size(taskEventList));
for ii = 1:length(taskEventList)
  fixDirSorted{ii} = taskData.fixTime(strcmp(taskData.taskEventIDs, taskEventList{ii}));
%   fixDirSorted = [fixDirSorted; fixDurStim];
end

spikesBinned = deEventArrayData(spikesByEventBinned);
spikesBinnedFixAlign = deEventArrayData(spikesByEventBinnedFixAlign);

for event_i = 1:length(spikesByEventBinned)
  fixTimesEvent = fixDirSorted{event_i};
  startInd = preStimTime - fixTimesEvent;
  for chan_i = 1:length(spikesByEventBinned{event_i})
    for unit_i = 1:length(spikesByEventBinned{event_i}{chan_i})
      
      % Find the mean firing rate for the entire trial, starting prior to the
      % fix cue, ending just after the reward.
      eventUnitData = spikesByEventBinned{event_i}{chan_i}{unit_i};
      
      % for each trial, extract the correct slice
      for trial_i = 1:size(eventUnitData,1)
        % Pull the trial raster from the appearance of the fix dot to just
        % after the reward
        trialRaster = eventUnitData(trial_i, startInd(trial_i):rewardEndInd);
        
        % Identify mean rate of unit on this trial
        meanRateTrial = sum(trialRaster)/length(trialRaster) * 1000;
        
        % Find spike times
        trialSpikeTimes = find(trialRaster);
        
        % count ISIs
        trialSpikeISIs = diff(trialSpikeTimes);
        
        % Find the spike rate for each pair of spikes.
        spikeRatePerPair = 2./trialSpikeISIs .* 1000;
        
        % Find the first one above the mean spike rate for the trials, make
        % this the initial value of T.
        firstPairAboveMean = find(spikeRatePerPair >= meanRateTrial, 1, 'first');
        initialT = trialSpikeISIs(firstPairAboveMean);
        
        % Add all subsequent spikes and their times to the train
        tVec = cumsum(trialSpikeISIs(firstPairAboveMean:end));
        
        % Calculate the surprisal index for these arrays. first by
        % calculating the chance of a poisson process producing it.
        spikesPerTrain = (1:length(tVec)) + 1;
        Pvec = spikeTrainProb(tVec, meanRateTrial, spikesPerTrain);

        
        
      end
      
      [binnedData, binStarts, binEnds] = binStuff(spikesBinnedFixAlign{chan_i}{unit_i}(:,prefixWin), 100, 25);
      
      % See at which point the sequence of bins rises above baseline (determined
      % by the fixAligned thing.
      
      % FixAligned spikes - use these to check for latency to fix dot.
      
      % stimAligned - use these to check for latency to stimulus.
      
    end
  end
end

end
function selTable = latencyAnalysisBaseline(selTable, psthByEventBinned, psthByEventBinnedFixAlign, psthParams, taskData)

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

preStimTime = psthParams.psthPre;% + psthParams.movingWin(1)/2;

% Determine the fixation duration for each trial
taskEventList = taskData.taskEventList;
latestReward = round(max(taskData.rewardTimePerTrial));
rewardEndInd = latestReward + preStimTime+ 50;

fixDirSorted = [];
% fixDirSorted = cell(size(taskEventList));
for ii = 1:length(taskEventList)
  fixDurStim = taskData.fixTime(strcmp(taskData.taskEventIDs, taskEventList{ii}));
  fixDirSorted = [fixDirSorted; fixDurStim];
end
fixDurMean = mode(fixDirSorted);

% psthBinned = deEventSpikeData(psthByEventBinned);
% psthBinnedFixAlign = deEventSpikeData(psthByEventBinnedFixAlign);

stimLatencyArray = cell(size(selTable,1),1);
fixLatencyArray = nan(size(stimLatencyArray));
unitInd = 1;

for chan_i = 1:length(psthByEventBinned)
  for unit_i = 1:length(psthByEventBinned{chan_i})
    
    % Identify the baseline for the unit. The data below is fix aligned,
    % and uses psthParams.
    eventUnitBaselineData = psthByEventBinnedFixAlign{chan_i}{unit_i}(:, psthParams.psthPre-500:psthParams.psthPre);
    baselineDist = eventUnitBaselineData(:);
    baselineThres = prctile(baselineDist, 95);
    
    % Establish latency with respect to fixation point
    fixationData = mean(psthByEventBinned{chan_i}{unit_i}(:,psthParams.psthPre-fixDurMean:psthParams.psthPre));
    fixLatencyUnit = find(fixationData > baselineThres, 1);
    if ~isempty(fixLatencyUnit)
      fixLatencyArray(unitInd) = fixLatencyUnit;
    end
    
    stimPresData = psthByEventBinned{chan_i}{unit_i}(:,psthParams.psthPre:psthParams.psthPre+psthParams.psthImDur);
    stimPresLatency = nan(size(taskEventList));
    for stim_i = 1:size(stimPresData,1)
      aboveBaseline = stimPresData(stim_i,:) > baselineThres;
      if any(aboveBaseline)
        % Find the threshold crossing
        baselineThresCross = find(aboveBaseline, 1);
        peaksInData = findpeaks(stimPresData(stim_i,:));
        
        % Find the first peak post crossing
        firstPeakPostCross = peaksInData.loc(find(peaksInData.loc > baselineThresCross, 1));
        
        % Store
        if ~isempty(firstPeakPostCross)
          stimPresLatency(stim_i) = firstPeakPostCross;
        end
        
      end
    end
    
    stimLatencyArray{unitInd} = stimPresLatency;
    unitInd = unitInd + 1;
  end
end

stimLatencyMinArray = cellfun(@(x) min(x), stimLatencyArray);

% Plotting - generate a per stim histogram
if 0
figure()
hold on
for stim_i = 1:length(taskEventList)
  histogram(cellfun(@(x) x(stim_i), stimLatencyArray), 28)
end
legend(taskEventList)

end

% Store in the outputs
selTable.latency_min = stimLatencyMinArray;      % The first response to any stimulus
selTable.latency_Fix = fixLatencyArray;          % The response latency during the fixation period.
selTable.latency_perStim = stimLatencyArray;     % The latency on a per stimulus basis

end
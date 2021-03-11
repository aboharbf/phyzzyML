function spikePupilCorrStruct = spikePupilCorr(spikesByEvent, eyeDataStruct, taskData, spikeTimes, psthParams, eventIDs, spikeAlignParams, figStruct)
% Plots the correlation between the PSTH and eye dilation.
spikePupilCorrStruct = struct();
% Step 1 - turn the spikesByEvent into a similarly structured
% psthByEventTrial w/ structure {stim}(trial, time)

onlyStim = 0; % For pupil and saccade analysis, remove the pre and post stimulus time.

pupilByStim = eyeDataStruct.pupilByStim;
padding = psthParams.movingWin(1)/2 + 1;

if onlyStim
  % find stim start and end times
  stimStart = abs(psthParams.psthPre);
  stimEnd = stimStart + psthParams.psthImDur;
  
  % Modify structures contain pupil and spikeBins to accomadate.
  pupilByStim = cellfun(@(x) x(:, stimStart:stimEnd), pupilByStim, 'UniformOutput',0);
  
  % Change stimStart and end due to padding on ephys.
  stimStart = abs(psthParams.psthPre)+padding;
  stimEnd = stimStart + psthParams.psthImDur;
  for stim_i = 1:length(spikesByEvent)
    for chan_i = 1:length(spikesByEvent{stim_i})
      spikesByEvent{stim_i}{chan_i} = cellfun(@(x) x(:, stimStart:stimEnd), spikesByEvent{stim_i}{chan_i}, 'UniformOutput',0);
    end
  end
else
  % Remove padding from spikesByEvent
  for stim_i = 1:length(spikesByEvent)
    for chan_i = 1:length(spikesByEvent{stim_i})
      spikesByEvent{stim_i}{chan_i} = cellfun(@(x) x(:, padding:end-padding), spikesByEvent{stim_i}{chan_i}, 'UniformOutput',0);
    end
  end
end

smoothingWidth = psthParams.smoothingWidth;
lfpPaddedBy = psthParams.movingWin(1)/2 + 1;

if ~spikeTimes
  filterPoints = -3*smoothingWidth:3*smoothingWidth;
  smoothingFilter = exp(-1*filterPoints.^2/(2*smoothingWidth^2));
  smoothingFilter = smoothingFilter/sum(smoothingFilter);
end

% Method - make histograms of pupil values where spikes happen vs don't. t
% test between them

% Plot 1 - for each stimulus, check if unit's response was coupled to
% pupil.
figTitle = 'Pupil Spike vs No Spike Histogram';
h = figure('Name',figTitle,'NumberTitle','off','units', 'normalized', 'Position', [0 0 1 0.9]);
objTbGrp = uitabgroup('Parent', h);
spikePupilSigArray = cell(length(spikesByEvent{1}), 2);

for chan_i = 1:length(spikesByEvent{1})
  [spikePupilSigArray{chan_i, :}] = deal(zeros(length(spikesByEvent), length(spikesByEvent{1}{chan_i})));
end

for stim_i = 1:length(spikesByEvent)
  objTab = uitab('Parent', objTbGrp, 'Title', eventIDs{stim_i});
  axesH = axes('parent', objTab);
  pupilStimData = pupilByStim{stim_i};
  %plotInd = 1;
  
  % Plot stuff
  histHandle = histogram(pupilStimData);
  hold on
  pupilDataMean = nanmean(nanmean(pupilStimData));
  maxBin = max(histHandle.Values);
  legendHandles = plot([pupilDataMean pupilDataMean], [0 maxBin*1.05], 'color','k', 'Linewidth', 3);
  legendLabels = {sprintf('%s = Grand Mean', num2str(pupilDataMean,3))};
  for chan_i = 1:length(spikesByEvent{stim_i})
    meanLineHandles = gobjects(length(spikesByEvent{stim_i}{chan_i}),1);
    meanLabels = cell(length(spikesByEvent{stim_i}{chan_i}),1);
    for unit_i = 1:length(spikesByEvent{stim_i}{chan_i})
      pupilStimCopy = pupilStimData;
      spikeBins = spikesByEvent{stim_i}{chan_i}{unit_i}(:,lfpPaddedBy:end-lfpPaddedBy);
      pupilSpikeVals = [];
      
      % For every trial, pull pupil values where spikes occur, replacing
      % them with NaNs.
      for trial_i = 1:size(spikeBins,1)
        spikeInds = find(spikeBins(trial_i,:));
        pupilSpikeVals = [pupilSpikeVals, pupilStimData(trial_i, spikeInds)];
        pupilStimCopy(trial_i,spikeInds) = deal(NaN);
      end
      
      % Plot Stuff
      grandMeanSubset = nanmean(pupilStimCopy(:));
      unitMean = nanmean(pupilSpikeVals);
      [~, B, ~, stats] = ttest2(pupilSpikeVals, pupilStimCopy(:));
      histogram(pupilSpikeVals);
      meanLineHandles(unit_i) = plot([unitMean unitMean], [0 maxBin], 'Linewidth', 3);
      meanLabels{unit_i} = sprintf('%s = Ch%dU%d (%s)', num2str(unitMean,3), chan_i, unit_i, num2str(B,3));
      cohensD = (unitMean - grandMeanSubset)/stats.sd;
      
      if B < 0.05
        text(nanmean(pupilSpikeVals), maxBin*1.01, '*', 'Fontsize', 14)
      end
      
      % Store significance test things
      spikePupilSigArray{chan_i, 1}(stim_i, unit_i) = B;
      spikePupilSigArray{chan_i, 2}(stim_i, unit_i) = cohensD;

    end
    
    % Add Legends
    legendHandles = [legendHandles; meanLineHandles];
    legendLabels = [legendLabels; meanLabels];
  end
  
  % Title, Legends
  title(sprintf('Pupil Dilations for %s', eventIDs{stim_i}))
  legend(legendHandles, legendLabels);
end

% Save figure
saveFigure(figStruct.figDir, figTitle, [], figStruct, figStruct.figTag);


% Plot 2 - for all stimulus together, check Unit-Pupil coupling
spikePupilAllStim = cell(length(spikesByEvent{1}), 2);
for chan_i = 1:length(spikesByEvent{1})
  [pVec, dVec] = deal(zeros(1, length(spikesByEvent{1}{chan_i})));
  for unit_i = 1:length(spikesByEvent{1}{chan_i})
    % Generate a figure which describes the unit's coupling to pupil
    % diameter.
    figTitle = sprintf('Pupil Spike vs No Spike Histogram, %s', figStruct.channelUnitNames{chan_i}{unit_i});
    h = figure('Name', figTitle, 'NumberTitle', 'off', 'units', 'normalized', 'Position', [0 0 1 0.9]);
    
    [spikingBins, nonSpikingBins] = deal([]);
    for stim_i = 1:length(spikesByEvent)
      spikeEventData = spikesByEvent{stim_i}{chan_i}{unit_i};
      pupilEventData = pupilByStim{stim_i};
      for trial_i = 1:size(spikeEventData, 1)
        spikeEventTrial = spikeEventData(trial_i, :);
        nonSpikingBins = [nonSpikingBins, pupilEventData(trial_i, ~logical(spikeEventTrial))];
        loopSub = 0;
        while any(spikeEventTrial)
          % Keep looping until spikeEventTrial is empty.
          spikeEventTrial = max(spikeEventTrial-loopSub,0);
          if any(spikeEventTrial)
            spikingBins = [spikingBins, pupilEventData(trial_i, logical(spikeEventTrial))];
            loopSub = 1;
          end
        end
      end
      
    end
    
    if ~any([isempty(spikingBins), isempty(nonSpikingBins)])
      % Make a histogram of the results
      axesH = axes(h);
      nSHist = histogram(nonSpikingBins);
      hold on
      sHist = histogram(spikingBins, nSHist.NumBins);
      
      % Add lines at the means
      nonSpikeMean = nanmean(nonSpikingBins);
      spikeMean = nanmean(spikingBins);
      
      maxBin = max(nSHist.Values);
      grandMeanHand = plot([nonSpikeMean nonSpikeMean], [0 maxBin*1.05], 'color','k', 'Linewidth', 3);
      grandMeanLabel = {sprintf('%s = Grand Mean', num2str(nonSpikeMean, 3))};
      
      spikeMeanHand = plot([spikeMean spikeMean], [0 maxBin*1.05], 'color','b', 'Linewidth', 3);
      spikeMeanLabel = {sprintf('%s = Unit Spike Mean', num2str(spikeMean, 3))};
      
      % T test, Add Star if Significant
      [A, pVal, ~, stats] = ttest2(nonSpikingBins, spikingBins);
      cohensD = (spikeMean - nonSpikeMean)/stats.sd;
      if A
        text(spikeMean, spikeMeanHand.YData(2)*1.00, '*', 'Fontsize', 14)
      end
      
      % Title, Legends
      title(sprintf('Pupil Spike Dilations for %s activity, p = %s, Cohens D = %s, ', figTitle, num2str(pVal, 3), num2str(cohensD, 3)))
      legend([grandMeanHand, spikeMeanHand], [grandMeanLabel, spikeMeanLabel]);
      
      % Store outputs
      pVec(unit_i) = pVal;
      dVec(unit_i) = cohensD;
      
      % Save figure
      saveFigure(figStruct.figDir, figTitle, [], figStruct, figStruct.figTag);
    else
      [pVec(unit_i), dVec(unit_i)] = deal(NaN);
    end
  end
  
  % Store outputs
  spikePupilAllStim{chan_i, 1} = pVec;
  spikePupilAllStim{chan_i, 2} = dVec;
  
end

spikePupilCorrStruct.spikePupilUnitStim.pVals = spikePupilSigArray(:,1);
spikePupilCorrStruct.spikePupilUnitStim.dVals = spikePupilSigArray(:,2);
spikePupilCorrStruct.spikePupilUnitAllStim.pVals = spikePupilAllStim(:,1);
spikePupilCorrStruct.spikePupilUnitAllStim.dVals = spikePupilAllStim(:,2);

end

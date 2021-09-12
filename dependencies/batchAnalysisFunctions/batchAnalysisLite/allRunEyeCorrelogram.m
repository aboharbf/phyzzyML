function allRunEyeCorrelogram(spikePathBank, params, figStruct)
% A function which gauges the amount of correlation across all trials for
% each stimulus.
% Inputs:
% - spikePathBank: Typical structure
% - params: spikePaathLoadParams
% - figStruct: typical
% Outputs:
% - Figures
%     - Per stimulus image, showing correlation across all trials of a
%     particular stimulus.
%     - Rows/Columns labeled with runs.

paradigmList = unique(spikePathBank.paradigmName);
monkeyArray = {'Mo', 'Sam', 'Combo'};

% make outputDir
outputDir = fullfile(params.outputDir, 'eyeCorr');

for m_i = 1:length(monkeyArray)
  for par_i = 1
    
    % Focus on the specific subset of data to process
    pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
    if ~strcmp(monkeyArray{m_i}, 'Combo')
      mInd = contains(spikePathBank.Properties.RowNames, monkeyArray{m_i});
    else
      mInd = true(size(pInd));
    end
    spikePathBankParadigm = spikePathBank(pInd & mInd, :);
    
    % Extract the relevant data from each run
    [eyeInByEventPerRun, psthParamsPerRun, eventIDsPerRun] = spikePathLoad(spikePathBankParadigm, {'eyeInByEvent', 'psthParams', 'eventIDs'}, params.spikePathLoadParams);
    stimStartInd = psthParamsPerRun{1}.psthPre;
    stimEndInd = stimStartInd + psthParamsPerRun{1}.psthImDur;
    
    % Create a structure which gathers all the trials for a particular
    % stimulus into the same cell
    
    uniqueIDs = unique([eventIDsPerRun{:}]);
    uniqueIDsPlot = strrep(uniqueIDs, '_', ' ');
    [allTrialeyeByEvent, runNameByEvent] = deal(cell(size(uniqueIDs)));
    for uni_i = 1:length(uniqueIDs)
      % Identify which data belongs to the stimulus per run.
      event2Find = uniqueIDs{uni_i};
      eventIndex = cellfun(@(x) find(strcmp(x, event2Find)), eventIDsPerRun, 'UniformOutput', false);
      run2Check = find(~cellfun('isempty', eventIndex));
      
      % Pull and stack the relevant vectors, take the relevant slice of
      % data.
      eyeInByStimPerRun = cell(size(run2Check));
      for run_i = 1:length(run2Check)
        eyeInByStimPerRun{run_i} = eyeInByEventPerRun{run2Check(run_i)}{eventIndex{run2Check(run_i)}}(:, :, stimStartInd:stimEndInd);
      end
      
      % Store into the larger structure
      allTrialeyeByEvent{uni_i} = eyeInByStimPerRun;
      runNameByEvent{uni_i} = spikePathBankParadigm.Row(run2Check, :);
      
    end
    
    %sampRate = 1/eyeCalParams.samplingRate;
    [eventCorr, trialPerEvent] = deal(nan(size(runNameByEvent)));
    
    %Cycle through "by Event" eye data"
    for event_i = 1:length(allTrialeyeByEvent)
      
      % Pull relevant data
      eyeInByStim = allTrialeyeByEvent{event_i};
      trialPerRun = cellfun(@(x) size(x, 2), eyeInByStim);
      trialPerEvent(event_i) = sum(trialPerRun);
      eyeInByStim = cat(2, eyeInByStim{:});
      
      eyeInByStim(isnan(eyeInByStim)) = 0;
      eyeXMat = squeeze(eyeInByStim(1, :, :));
      eyeYMat = squeeze(eyeInByStim(2, :, :));
      
      % eventPower(event_i) = (mean(bandpower(eyeXMat)) + mean(bandpower(eyeYMat)))/2;
      
      %Get rid of blinks
      eyeXMat(eyeXMat > 12) = nan;
      eyeYMat(eyeYMat > 12) = nan;
      
      % Round values since we don't actually have that level of percition
      eyeXMat = round(eyeXMat,2)';
      eyeYMat = round(eyeYMat,2)';
      
      %Calculate intertrial corr for this stimulus
      eyeXCorr = corr(eyeXMat, 'rows', 'pairwise'); % NaNs are excluded on a column by column basis.
      eyeYCorr = corr(eyeYMat, 'rows', 'pairwise');
      
      %Store correlations across trials and variance in each signal.
      eyeCorrAll = cat(3, eyeXCorr, eyeYCorr);
      eyeMatCorr = mean(eyeCorrAll, 3);
      eventCorr(event_i) = mean(eyeCorrAll(:));
      
      %Plotting Time
      eyeTitle = sprintf('Eye signal Intertrial correlation - %s', uniqueIDsPlot{event_i});
      figure('Name', eyeTitle, 'NumberTitle', 'off', 'units', 'normalized', 'position', [0.1 0.1 0.54 0.69]);
      subplot(4,1, [1:3])
      imageH = imagesc(eyeMatCorr);
      axesH = imageH.Parent;
      title(eyeTitle)
      colorbar
      
      %Chop it up
      tLen = length(eyeMatCorr);
      lineSpots = cumsum(trialPerRun);
      for line_i = 1:length(lineSpots)
        lineInd = lineSpots(line_i) + 0.5;
        %Draw vertical line
        line([lineInd lineInd], [0 tLen], 'Linewidth',.5, 'color', 'k')
        %Draw horzontal line
        line([0 tLen], [lineInd lineInd], 'Linewidth',.5, 'color', 'k')
      end
      
      [axesH.YTick, axesH.XTick] = deal(lineSpots - (trialPerRun/2));
      axesH.YTickLabel = arrayfun(@(x) sprintf('%d:%s', x, runNameByEvent{event_i}{x}(6:end)), 1:length(runNameByEvent{event_i}), 'UniformOutput', false);
      axesH.XTickLabel = 1:length(runNameByEvent{event_i});
      
      axH3 = subplot(4,1,4);
      binSize = 100;
      stepSize = 25;
      [corrDataX, binStarts, binEnds] = slidingWindowCorr(eyeXMat', binSize, stepSize);
      [corrDataY, ~, ~] = slidingWindowCorr(eyeYMat', binSize, stepSize);
      binLabels = mean(cat(1, binStarts, binEnds), 1);
      plot(mean(cat(2, corrDataX, corrDataY), 2))
      axH3.XTick = 1:10:length(corrDataY);
      axH3.XTickLabel = binLabels(1:10:length(corrDataY));
      axH3.XLim = [1 length(corrDataY)];
      
      % Save
      saveFigure(outputDir, sprintf('EyeCorr_%s', uniqueIDs{event_i}), [], figStruct, []);
      
      % Simpiler map - take averages across all the runs for this stim
      eyeMatCorrSimp = zeros(length(runNameByEvent{event_i}));
      sampRang = [0; lineSpots];
      for row_i = 1:size(eyeMatCorrSimp,1)
        for col_i = 1:size(eyeMatCorrSimp,2)
          eyeMatCorrSimp(row_i, col_i) = mean(mean(eyeMatCorr(sampRang(row_i)+1:sampRang(row_i+1), sampRang(col_i)+1:sampRang(col_i+1))));
        end
      end
      
      % Plot
      eyeAvgTitle = sprintf('Eye signal inter-trial correlation, Averaged across Runs - %s', uniqueIDsPlot{event_i});
      figure('Name', eyeAvgTitle, 'NumberTitle', 'off', 'units', 'normalized', 'position', [0.1 0.1 0.54 0.69]);
      imgH = imagesc(eyeMatCorrSimp);
      axeH2 = imgH.Parent;
      title(eyeAvgTitle)
      colorbar
      
      [axeH2.YTick, axeH2.XTick] = deal(1:length(runNameByEvent{event_i}));
      axeH2.YTickLabel = arrayfun(@(x) sprintf('%d:%s', x, runNameByEvent{event_i}{x}(6:end)), 1:length(runNameByEvent{event_i}), 'UniformOutput', false);
      axeH2.XTickLabel = 1:length(runNameByEvent{event_i});
      
      saveFigure(outputDir, sprintf('EyeCorrAvg_%s', uniqueIDs{event_i}), [], figStruct, []);
      
    end
  end
end

end
function eyeDataStruct = pupilDilation(analogInByEvent, psthParams, eventLabels, eyeDataStruct, catIndStruct, figStruct)

%blinkVals = taskData.eyeCal.blinkVals;
lfpPaddedBy = (psthParams.movingWin(1)/2)+1;
blinkPad = 15;      % May need to be slightly widdened. Consider simply tiling 1st value.
pupilImg = cell(length(analogInByEvent),1);

returnBlinks = 0;
fltord = 60;
lowpasfrq = 30;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);     % low pass filter internal, determined by analysisParam.
maxEndTime = size(analogInByEvent{1}, 4);

for stim_i = 1:length(analogInByEvent)
%   figure()
%   hold on
  % Blink processing
  stimData = squeeze(analogInByEvent{stim_i}(:,3,:,:));
  blinkThres = mean(mean(stimData)) - (std(std(stimData)) * 2.5);
  stimThresBlinks = stimData < blinkThres;
  stimDataDiff = diff(stimData, 1, 2);
  
  for trial_i = 1:size(stimData,1)
    stimDataOrig = stimData(trial_i,:);
    % Plot original trace
%     subplot(2,1,1)
%     hold on
%     plot(stimDataOrig)
    
    % Identify blink locations
    blinkDiff = stimDataDiff(trial_i, :);
    blinkEnds = find(blinkDiff > 0.3) + 1;
    blinkStarts = find(blinkDiff < -0.3) + 1;
    blinkMid = find(stimThresBlinks(trial_i,:));
    blinkInds = unique(sort([blinkMid, blinkStarts, blinkEnds]));

    if ~isempty(blinkInds)
      blinkTimes = BehavioralIndexPlus(blinkInds, maxEndTime);      
      %blinkTimes = BehavioralIndex(blinkInds);
      for blink_i = 1:size(blinkTimes, 2)
        % Slide end further if needed
        blinkEnd = blinkTimes(2, blink_i);
        if (blinkEnd ~= maxEndTime)
          blinkEnd = blinkEnd + find(blinkDiff(blinkEnd:end) < 0, 1);
          if ~isempty(blinkEnd)
            blinkTimes(2, blink_i) = blinkEnd;
          end
        end
        % Use this loop to remove blinks from the data.
        preBlinkInd = max(blinkTimes(1, blink_i) - blinkPad, 1);
        postBlinkInd = min(blinkTimes(2, blink_i) + blinkPad, maxEndTime);
        % Interpolate between the two points found, fill in the blink with
        % this value. If close to the edge, fill with the non-edge point.
        if (preBlinkInd == 1) || (postBlinkInd == maxEndTime)
          if preBlinkInd == 1
            newVal = stimData(trial_i, postBlinkInd);
          else
            newVal = stimData(trial_i, preBlinkInd);
          end
          stimData(trial_i, preBlinkInd:postBlinkInd) = deal(newVal);
        else
          periBlinkVals = [stimData(trial_i, preBlinkInd), stimDataOrig(postBlinkInd)];
          %Overwrite the blink period with data from the sides. Interp b/t the
          %points for now.
          newInd = linspace(1,2, length(stimData(trial_i, preBlinkInd:postBlinkInd)));
          stimData(trial_i, preBlinkInd:postBlinkInd) = interp1([1 2], periBlinkVals, newInd);
        end
        
      end
    end
    
    % Smooth data
    stimDataSmooth = filtfilt(flt, 1, stimData(trial_i, :));
    
    % Blink Management post smoothing
    if ~isempty(blinkInds) && returnBlinks
      % Return blinks to their place
      for blink_i = 1:size(blinkTimes, 2)
        stimDataSmooth(blinkTimes(1,:): blinkTimes(2,:)) = stimDataOrig(blinkTimes(1,:): blinkTimes(2,:));
      end
    elseif ~isempty(blinkInds) && ~returnBlinks
      % NaN out the regions of blink, since pupil is not seen.
      for blink_i = 1:size(blinkTimes, 2)
        stimDataSmooth(blinkTimes(1,:): blinkTimes(2,:)) = deal(nan());
      end
    end
    
    % Overwrite original
    stimData(trial_i,:)=  stimDataSmooth;
    
    % Plot
%     subplot(2,1,2)
%     plot(stimDataSmooth)
%     hold on
  end
%   title(eventIDs{stim_i})
  % Chop off the buffer and store into the output structure.
  pupilImg{stim_i} = stimData(:, lfpPaddedBy:end-lfpPaddedBy);
end

% Plot distributions Catagory vs non-Catagory
catagoryPlot = {'socialInteraction'};
catStats = cell(length(catagoryPlot),1);
for cat_i = 1:length(catagoryPlot)
  % see if the stimulus is present
  stimInd = catIndStruct.catIndMat(:, strcmp(catagoryPlot{cat_i}, catIndStruct.categoryList));
  if length(unique(stimInd)) > 1
    % Prepare figure
    figTitle = sprintf('Catagory based Pupil dilations, %s vs Non-%s', catagoryPlot{cat_i}, catagoryPlot{cat_i});
    h = figure('Name',figTitle,'NumberTitle','off','units','normalized','outerposition',figStruct.figPos);
    
    % Split pupil counts;
    catPupil = pupilImg(stimInd);
    nonCatPupil = pupilImg(~stimInd);
    
    % Compile, take mean, and plot
    catPupil = vertcat(catPupil{:});
    catPupil = catPupil(:);
    nonCatPupil = vertcat(nonCatPupil{:});
    nonCatPupil = nonCatPupil(:);
    
    catPupilMean = nanmean(catPupil);
    nonCatPupilMean = nanmean(nonCatPupil);
    
    hist2Hand = histogram(nonCatPupil);
    hold on
    hist1Hand = histogram(catPupil);
    hist1Mean = plot([catPupilMean catPupilMean], [0 h.Children.YLim(2)], 'lineWidth', 3);
    hist2Mean = plot([nonCatPupilMean nonCatPupilMean], [0 h.Children.YLim(2)], 'lineWidth', 3, 'color', 'k');
    maxBin = max([hist1Hand.Values, hist2Hand.Values]);
    h.Children.YLim(2) = maxBin * 1.05;
    
    % Significance test
    [A, pVal, ~, stats] = ttest2(catPupil, nonCatPupil);
    if A
      text(catPupilMean, maxBin*1.01, '*', 'Fontsize', 14)
    end
    cohensD = (nanmean(catPupil) - nanmean(nonCatPupil))/stats.sd;
    
    % Legends, Titles
    title(sprintf('%s (p = %s, Cohens D = %s)', figTitle, num2str(pVal,3), num2str(cohensD, 3)))
    h.Children.XLabel.String = 'Raw Pupil values';
    h.Children.YLabel.String = 'Frequency';
    
    % Save the figure
    if figStruct.saveFig
      % saveFigure( outDir, filename, figData, saveFig, exportFig, saveData, varargin )
      saveFigure(figStruct.figDir, figTitle, [], figStruct, figStruct.figTag );
    end
    
    % Save to Output Struct
    catStats{cat_i} = {pVal, cohensD};
  end
end

% Plot adaptation curves
stimStart = psthParams.psthPre;
stimEnd = psthParams.psthPre + psthParams.psthImDur;
[meanVec, SDVec] = deal(cell(length(eventLabels),1));

figTitle = 'Pupil adaptation';
h = figure('Name',figTitle,'NumberTitle','off','units','normalized','outerposition',figStruct.figPos);
objTbGrp = uitabgroup('Parent', h);

stimPupilAdaptLine = cell(length(pupilImg),1);
for stim_i = 1:length(pupilImg)
  objTab = uitab('Parent', objTbGrp, 'Title', eventLabels{stim_i});
  axesH = axes('parent', objTab);
  pupilStimData = pupilImg{stim_i};
  
  % Plot 1 - all trials and their collective means
  subplot(2,1,1)
  hold on
  plot(1:length(pupilStimData), pupilStimData)
  plot(1:length(pupilStimData), nanmean(pupilStimData, 1), 'LineWidth', 3, 'color', 'r')
  xlim([0 length(pupilStimData)])
  plot([stimStart stimStart], ylim(), 'color', 'k', 'LineWidth', 4);
  plot([stimEnd stimEnd], ylim(), 'color', 'k', 'LineWidth', 4);
  
  title(sprintf('Individual Trials - %s', eventLabels{stim_i}))
  legendLabels = [arrayfun(@(x) ['Trial ' num2str(x)], 1:size(pupilStimData,1), 'UniformOutput', 0), 'Mean']; 
  legend(legendLabels, 'location', 'northeastoutside', 'AutoUpdate', 'off');  
  
  % Plot 2 - Is there a change over trial period?
  subplot(2,1,2)
  hold on
  [meanV, SDV] = deal(zeros(size(pupilStimData,1), 1));
  [ally, allx] = deal([]);
  for trial_i = 1:size(pupilStimData,1)
    meanV(trial_i) = nanmean(pupilStimData(trial_i, :));
    SDV(trial_i) = nanstd(pupilStimData(trial_i, :));
    
    ydata = pupilStimData(trial_i, :);
    xdata = ones(size(ydata)) * trial_i;
    
    ally = [ally, ydata];
    allx = [allx, xdata];

  end
  plot(xdata, ydata, 'r*', 'MarkerSize', 1, 'LineWidth', 1);
  hold on
  meanVec{stim_i} = meanV;
  SDVec{stim_i} = SDV;
  
  % Use this data to fit a linear model
  linMod = fitlm(allx, ally);
  plot(linMod);
  [stimPupilAdaptLine{stim_i}, lineVals] = deal(linMod.Coefficients.Estimate);
  
  title(sprintf('Pupilary Adaptation : x1 = %s, int = %s', num2str(lineVals(2), 2), num2str(lineVals(1), 2)))
end

% 
% % All Stim adaptation plot
% objTab = uitab('Parent', objTbGrp, 'Title', 'All Stim');
% axesH = axes('parent', objTab);
% trialCounts = cellfun('length', meanVec);
% [allMeanVec, allSDVec] = deal(nan(length(meanVec), max(trialCounts)));
% for stim_i = 1:size(allMeanVec)
%   allMeanVec(stim_i, 1:trialCounts(stim_i)) = meanVec{stim_i};
%   allSDVec(stim_i, 1:trialCounts(stim_i)) = SDVec{stim_i};
% end
% mseb(1:size(allMeanVec,2), nanmean(allMeanVec), nanmean(allSDVec));

% Once done, Save the figure
if figStruct.saveFig
  % saveFigure( outDir, filename, figData, saveFig, exportFig, saveData, varargin )
  saveFigure(figStruct.figDir, figTitle, [], figStruct, figStruct.figTag);
end

eyeDataStruct.pupilByStim = pupilImg;
eyeDataStruct.pupilStats.catComp = catagoryPlot;
eyeDataStruct.pupilStats.catStats = catStats;
eyeDataStruct.pupilStats.pupilAdapt = stimPupilAdaptLine;

end

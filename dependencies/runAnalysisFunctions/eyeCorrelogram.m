function eyeCorrelogram(eyeInByEvent, psthParams, eventLabels, figStruct)

outDir = figStruct.figDir;
stimStartInd = psthParams.psthPre;
stimEndInd = stimStartInd + psthParams.psthImDur;

%sampRate = 1/eyeCalParams.samplingRate;
[eventCorr, trialPerEvent, eventPower] = deal(nan(length(eyeInByEvent), 1));

%Cycle through "by Event" eye data"
for event_i = 1:length(eyeInByEvent)
  assert(size(eyeInByEvent{event_i}, 3) > 1, 'Not enough trials to do correlation')
  eyeInByStim = eyeInByEvent{event_i}(:,:,stimStartInd:stimEndInd);
  trialPerEvent(event_i) = size(eyeInByStim, 2);
  eyeInByStim(isnan(eyeInByStim)) = 0;
  eyedat = cell(1, size(eyeInByStim,2));
  [eyeXMat, eyeYMat] = deal(nan(size(eyeInByStim,3), size(eyeInByStim,2)));
  
  for trial_i = 1:size(eyeInByStim, 2)
    eyedat{trial_i} = [squeeze(eyeInByStim(1,trial_i,:))' ; squeeze(eyeInByStim(2,trial_i,:))'];
    eyeXMat(:,trial_i) = squeeze(eyeInByStim(1,trial_i,:));
    eyeYMat(:,trial_i) = squeeze(eyeInByStim(2,trial_i,:));
  end
  
  eventPower(event_i) = (mean(bandpower(eyeXMat)) + mean(bandpower(eyeYMat)))/2;
  
  %Get rid of blinks
  eyeXMat(eyeXMat > 12) = nan;
  eyeYMat(eyeYMat > 12) = nan;
  
  %Calculate intertrial corr for this stimulus
  eyeXCorr = corr(eyeXMat,'rows', 'complete');
  eyeYCorr = corr(eyeYMat,'rows', 'complete');
  
  %Store correlations across trials and variance in each signal.
  eventCorr(event_i) = (mean(mean(eyeXCorr)) + mean(mean(eyeYCorr)))/2;
  
  
  %Store for large correlation map.
  if event_i == 1
    eyeXMatTot = eyeXMat;
    eyeYMatTot = eyeYMat;
  else
    eyeXMatTot = [eyeXMatTot eyeXMat];
    eyeYMatTot = [eyeYMatTot eyeYMat];
  end
end

eyeXMatCorr = corr(eyeXMatTot,'rows', 'complete');
eyeYMatCorr = corr(eyeYMatTot,'rows', 'complete');

eyeMatCorr =(eyeXMatCorr + eyeYMatCorr)/2;

%Plotting Time
eyeTitle = sprintf('Eye signal intertrial correlation');
eyeSig = figure('Name',eyeTitle,'NumberTitle','off');
imagesc(eyeMatCorr)
title(eyeTitle)
colorbar

%Chop it up
tLen = length(eyeMatCorr);
for event_i = 1:length(trialPerEvent)
  lineInd = sum(trialPerEvent(1:event_i)) + .5;
  %Draw vertical line
  line([lineInd lineInd], [0 tLen], 'Linewidth',.5, 'color', 'k')
  %Draw horzontal line
  line([0 tLen], [lineInd lineInd], 'Linewidth',.5, 'color', 'k')
end
axis off

% Label each row
for event_i = 1:length(trialPerEvent)
  labelYInd = sum(trialPerEvent(1:event_i)) - (trialPerEvent(event_i)/2);
  labelXInd = -trialPerEvent(1)*1.5;
  str = [num2str(event_i) ':' eventLabels{event_i}];
  %Draw vertical line
  text(labelXInd,labelYInd,str)
end

%Label the row of power signals.
labelYInd = sum(trialPerEvent(1:event_i)) + sum(trialPerEvent(event_i))/2;
text(labelXInd*1.5,labelYInd,'Power in Eye Signal','FontSize',14)

% Label the Columns
for event_i = 1:length(trialPerEvent)
  labelYInd = sum(trialPerEvent)+3;
  labelXInd = sum(trialPerEvent(1:event_i)) - (trialPerEvent(event_i)/2);
  str = num2str(event_i);
  pwr = num2str(eventPower(event_i));
  text(labelXInd,labelYInd,str,'FontSize',12)
  text(labelXInd,labelYInd+5,pwr,'FontSize',12)
end

figData = eyeMatCorr;
saveFigure(outDir, sprintf('EyeCorr_%s',figStruct.figTag), figData, figStruct, figStruct.figTag);


%Simpiler map
corrVec = nan(length(trialPerEvent));
for event_i = 1:length(trialPerEvent)
  startXInd = 1 + sum(trialPerEvent(1:event_i-1));
  endXInd = sum(trialPerEvent(1:event_i));
  for event_j = 1:length(trialPerEvent)
    startYInd = 1 + sum(trialPerEvent(1:event_j-1));
    endYInd = sum(trialPerEvent(1:event_j));
    corrVec(event_i,event_j) = mean(mean(eyeMatCorr(startXInd:endXInd,startYInd:endYInd)));
  end
end

eyeAvgTitle = sprintf('Eye signal inter-trial correlation, Averaged across stimuli');
eyeSigAvg = figure('Name',eyeAvgTitle,'NumberTitle','off');
imagesc(corrVec)
title(eyeAvgTitle)
colorbar
axis off

% Label each row
for event_i = 1:length(trialPerEvent)
  labelYInd = event_i;
  labelXInd = -3;
  str = [num2str(event_i) ':' eventLabels{event_i}];
  text(labelXInd,labelYInd,str,'FontSize',12)
end

% Label the Columns
for event_i = 1:length(trialPerEvent)
  labelYInd = length(trialPerEvent)+1;
  labelXInd = event_i;
  str = num2str(event_i);
  pwr = num2str(eventPower(event_i));
  text(labelXInd,labelYInd,str,'FontSize',12)
  text(labelXInd,labelYInd+3,pwr,'FontSize',12)
end

figData = corrVec;
saveFigure(outDir, sprintf('EyeCorrAvg_%s',figStruct.figTag), figData, figStruct, figStruct.figTag);

end

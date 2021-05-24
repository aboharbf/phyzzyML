function saccadeByStimLabeledImage(saccadeByStim, psthParams)

load('subEventAnalysisTesting_naturalSocial.mat', 'eyeDataStruct', 'psthParams')
saccadeByStim = eyeDataStruct.saccadeByStim;

% Variables - 
trialPerStim = cellfun(@(x) size(x, 1), saccadeByStim);
ISI = psthParams.ITI;
pre = psthParams.psthPre;
stimDur = psthParams.psthImDur;
post = psthParams.psthPost;

fixDotOn = ISI;
stimOn = pre;
stimOff = pre+stimDur;
stimEvents = stimOn:500:stimOff;
reward = pre+stimDur+200;

% X spots to label
xticks2Label = [fixDotOn, stimEvents, reward];
xLines = [fixDotOn, stimOn, stimOff, reward];
% 
saccadeByStim = vertcat(eyeDataStruct.saccadeByStim{:});

% Plot
figH = figure();
axesH = axes(figH);
imagesc(saccadeByStim);

% Fix Axes
axesH.XTick = xticks2Label;
axesH.XTickLabel = xticks2Label - stimOn;
hold on

for line_i = 1:length(xLines)
  plot([xLines(line_i), xLines(line_i)], ylim(), 'lineWidth', 3);
end

xticklabels 
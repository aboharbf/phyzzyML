function saccadeByStimLabeledImage(saccadeByStim, psthParams)

load('subEventAnalysisTesting_naturalSocial.mat', 'eyeDataStruct', 'psthParams', 'eventIDs')
saccadeByStim = eyeDataStruct.saccadeByStim;

% Variables - 
trialPerStim = cellfun(@(x) size(x, 1), saccadeByStim);
ISI = 400;
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
saccImg = imagesc(saccadeByStim);
% saccImg.AlphaData = (saccadeByStim ~= 1);

% Fix Axes
axesH.XTick = xticks2Label;
axesH.XTickLabel = xticks2Label - stimOn;
hold on

for line_i = 1:length(xLines)
  plot([xLines(line_i), xLines(line_i)], ylim(), 'lineWidth', 3);
end

legH = legend('Fixation Onset*', 'Stimulus Onset', 'Stimulus Offset', 'Mean Reward Delivery');
legH.AutoUpdate = 'off';
legH.Location = 'northeastoutside';

% Add Horizontal line
stimuliBreaks = cumsum(trialPerStim) + 0.5;
xlimSize = xlim();
xlim(xlimSize);
for ii = 1:length(stimuliBreaks)
  plot(xlimSize, [stimuliBreaks(ii), stimuliBreaks(ii)], 'linewidth', 2, 'color', 'k');
end

% Add Stim Labels
stimLabels = round(stimuliBreaks - trialPerStim/2);
yticks(stimLabels);
yticklabels(strrep(extractBefore(eventIDs, '.avi'), '_', ''));

title('Eye Behavior Plot (Purple = Fix, Blue = Pre Saccade, Green = Saccade, Yellow = Blink)')

disp('y');
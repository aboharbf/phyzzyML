function [ ] = rasterSimple(spikesByItem, pictureLabels, rasterAxes, psthParams, ISI)
%RASTER makes a raster plot in the current figure
%    Note: if no spikes on any trial of an image, that image will not appear in legend

preAlign = psthParams.psthPre;
postAlign = psthParams.psthPost;
imDur = psthParams.psthImDur;

if isempty(rasterAxes)
  rasterAxes = gca;  
end

spikesByItemFixed = cell(size(spikesByItem));
for ii = 1:length(spikesByItem)
  raw = spikesByItem{ii};
  raw = raw(:, psthParams.spikeBuffer:end-(psthParams.spikeBuffer+1));
  spikesByItemFixed{ii} = raw;
end
spikesByItem = spikesByItemFixed;

tmpSCHand = imsc(~vertcat(spikesByItem{:}), 'hicontrast');
tmpSCHand.Parent = rasterAxes;
delete(gcf)

xlim(rasterAxes, [-preAlign,imDur+postAlign]); 
hold(rasterAxes, 'on');
axis(rasterAxes, 'ij');
yLevel = -0.5; % accumulated trial index, sets height in raster
legendHandles = [];

for item_i = 1:length(spikesByItem)
  yLevelStart = yLevel;
  for trial_i = 1:size(spikesByItem{item_i}, 1)
    yLevel = yLevel + 1;
%     trialSpikes = find(spikesByItem{item_i}(trial_i,:)) - psthParams.psthPre;
%     yLvl1 = repmat(yLevel-0.4, size(trialSpikes));
%     yLvl2 = repmat(yLevel+0.4, size(trialSpikes));
%     plot(rasterAxes, [trialSpikes; trialSpikes],[yLvl1; yLvl2],'color', 'k');
  end
  
  %Plot relevant stimulus onset and offset marker
  h = plot(rasterAxes, [0 0],[yLevelStart+0.5 yLevel+0.5],'color', 'k');
  legendHandles = vertcat(legendHandles,h);
  plot(rasterAxes, [imDur imDur],[yLevelStart+0.5 yLevel+0.5],'color', 'k');
  if item_i ~= length(spikesByItem)
    plot(rasterAxes, [xlim(rasterAxes)],[yLevel+0.5 yLevel+0.5],'color','black', 'LineWidth', 2, 'LineStyle', '--')
  end
end
xlimits = xlim(rasterAxes);
ylim(rasterAxes, [0,yLevel+0.5]);
if xlimits(2) < imDur + ISI
  plot(rasterAxes, [imDur+ISI imDur+ISI],[0 yLevel],'b-');
end
xlabel(rasterAxes, 'Time after stimulus onset (ms)');
ylabel(rasterAxes, 'single trials');
[rasterAxes.YTick, rasterAxes.YTickLabel] = deal([]);
legend(rasterAxes, legendHandles, pictureLabels, 'Location','northeastoutside');
hold(rasterAxes, 'off');

tmpSCHand.YData = rasterAxes.YLim;
tmpSCHand.XData = rasterAxes.XLim;
ylim(rasterAxes, rasterAxes.YLim)
xlim(rasterAxes, rasterAxes.XLim)
end




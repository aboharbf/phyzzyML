function [ ] = rasterColorCoded(inputHandle, spikesByItem, pictureLabels, psthParams, ISI, channel_i, unit_i, eyeDataStruct, colorSwitch, spikePattern)
%RASTER makes a raster plot in the current figure of the specified unit and
%channel data. Uses information from the attendedObjData structure to color
%code the spikes based on what object the subject was looking at. 
%    Note: if no spikes on any trial of an image, that image will not appear in legend
% colorSwitch (in AnalysisParam): 1 = Spikes, 2 = Shaded Regions for attended Obj, 3 =
% Shaded region for Saccades, 4 = shaded region for pupil dil, no saccades.

% Color Spikes vs Color Shaded Region
%Set the figure handle as the current handle

if ~exist('spikePattern', 'var')
  spikePattern = 1;
end

switch class(inputHandle)
  case 'matlab.ui.Figure'
    set(0, 'CurrentFigure', inputHandle)
    rasterAxes = axes();
  case 'matlab.graphics.axis.Axes'
    rasterAxes = inputHandle;
end

fieldsUnpack = fields(eyeDataStruct);
for ii = 1:length(fieldsUnpack)
  eval(sprintf('%s = eyeDataStruct.%s;', fieldsUnpack{ii}, fieldsUnpack{ii}))
end

%Unpack Variables
preAlign = psthParams.psthPre;
postAlign = psthParams.psthPost;
imDur = psthParams.psthImDur;

xlim([-preAlign,imDur+postAlign]);
hold on;
axis ij

%Reshape attendedObjData to be correctly indexed.
if exist('attendedObjData', 'var')
  try
    
    attendObjInd = zeros(length(pictureLabels),1);
    for obj_ind = 1:length(pictureLabels)
      attendObjInd(obj_ind) = find(strcmp(attendedObjData.eventIDs, pictureLabels{obj_ind}));
    end
    
    attendedObjData.eventIDs = attendedObjData.eventIDs(attendObjInd);
    attendedObjData.attendedObjVect = attendedObjData.attendedObjVect(attendObjInd);
    attendedObjData.frameStartInd = attendedObjData.frameStartInd(attendObjInd);
    attendedObjData.frameEndInd = attendedObjData.frameEndInd(attendObjInd);
    attendedObjData.tracePlotData = attendedObjData.tracePlotData(attendObjInd);
    
    %attendedObjData.colorCode{end+1} = [0, 0, 0]; %Black for non-stim spikes.
    uniqueColorTrials = zeros(length(attendedObjData.colorCode),1);
    uniqueColorInds = zeros(length(attendedObjData.colorCode),3);
    
  catch
    attendObjInd = 1:length(pictureLabels);
  end
  
end

if colorSwitch == 1
  %Create the color index for all the spikes
  for item_i = 1:length(spikesByItem)
    for trial_i = 1:length(spikesByItem{item_i}{channel_i}{unit_i})
      %Initialize a vector for colors
      spikeTimes = spikesByItem{item_i}{channel_i}{unit_i}(trial_i).times;
      spikeColorTmp = zeros(length(spikeTimes),3);
      spikeTimes = spikeTimes - psthParams.latency; %Shift by latency.
      
      %Assign remanining spikes colors
      for spike_i = 1:length(spikeTimes)
        if ~(spikeTimes(spike_i) < 0 || spikeTimes(spike_i) > imDur)     %Pre and Post stimuli spikes are black.
          frameInd = find((spikeTimes(spike_i)-attendedObjData.frameStartInd{item_i}) > 0,1,'last');
          objAtt = attendedObjData.attendedObjVect{item_i}(trial_i, frameInd);
          colorInd = strcmp(attendedObjData.objList, objAtt);
          spikeColorTmp(spike_i,:) = attendedObjData.colorCode{colorInd};
          if sum(logical(uniqueColorTrials)) ~= sum(logical(uniqueColorTrials+colorInd)) %If adding the index makes a new logical vecotr, this is a new color.
            uniqueColorInds(colorInd,:) = [item_i trial_i spike_i]; %Used later to get handles for legend.
            uniqueColorTrials = uniqueColorTrials + colorInd;
          end
        end
      end
      spikeColor{item_i}{trial_i}.color = spikeColorTmp;
    end
  end
end

%Plot the Shaded area, if appropriate
try
  if colorSwitch == 2
    attObjData = cat(1, attendedObjData.tracePlotData{1:end});
    im = image(attObjData);
    imOffset = 0;
    im.XData = [0+imOffset psthParams.psthImDur-imOffset];
    im.AlphaData = 0.5;
    im.YData = [im.YData(1)-0.5 im.YData(2)-0.5];
    inputHandle.UserData.shadedAreaHandle = im;
  elseif colorSwitch == 3 ||  colorSwitch == 4
    if colorSwitch == 3
      saccadeData = vertcat(saccadeByStim{attendObjInd});
      im = imagesc(saccadeData);
    else
      pupilByStim = vertcat(pupilByStim{attendObjInd});
      im = imagesc(pupilByStim);
      colorbar('location','westoutside')
    end
    im.XData = [-preAlign,imDur+postAlign];
    im.YData = [im.YData(1)-0.5 im.YData(2)-0.5];
    im.AlphaData = 0.5;
    inputHandle.UserData.shadedAreaHandle = im;
  end
catch
  warning('colorSwitch failed')
end

%Plot spikes
yLevel = -0.5; % accumulated trial index, sets height in raster
legendHandles = [];
% spikeHandles = [];
trialLabels = [];
switch spikePattern
  case 1
    % Traditional representation of spikes as lines.
    for item_i = 1:length(spikesByItem)
      yLevelStart = yLevel;
      trialLabels = [trialLabels 1:length(spikesByItem{item_i}{channel_i}{unit_i})];
      for trial_i = 1:length(spikesByItem{item_i}{channel_i}{unit_i})
        yLevel = yLevel + 1;
        trialSpikes = spikesByItem{item_i}{channel_i}{unit_i}(trial_i);
        for spike_i = 1:length(trialSpikes.times)
          if colorSwitch == 1
            trialColors = spikeColor{item_i}{trial_i}.color;
            plot([trialSpikes.times(spike_i) trialSpikes.times(spike_i)],[yLevel-0.4 yLevel+0.4],'color', trialColors(spike_i,:));
          else
            plot([trialSpikes.times(spike_i) trialSpikes.times(spike_i)],[yLevel-0.4 yLevel+0.4],'color', 'black');
          end
        end
      end
      %Plot Stim start, stim end, and Bottom dividing bar
      h = plot([0 0],[yLevelStart+0.5 yLevel+0.5],'color', 'black', 'LineWidth', 2); % Stimulus Onset
      plot([imDur imDur],[yLevelStart+0.5 yLevel+0.5],'color', 'black', 'LineWidth', 2); % Stimulus Offset
      if item_i ~= length(spikesByItem)
        plot([xlim()],[yLevel+0.5 yLevel+0.5],'color','black', 'LineWidth', 2, 'LineStyle', '--'); %Stimuli dividing line
      end
      legendHandles = vertcat(legendHandles,h);
    end
    
  case 2
    % Spikes represented as single image. - Seems to not work atm.
    spikeImg = [];
    traceTime = preAlign + postAlign + imDur;
    spikeImgIndOffset = preAlign * 10;
    spiDi = 1;
    yLevel = -0.5; % accumulated trial index, sets height in raster
    
    for item_i = 1:length(spikesByItem)
      trialCount = length(spikesByItem{item_i}{channel_i}{unit_i});
      tmp = ones(trialCount, traceTime*10, 3);  %Grid with Trial rows and 1/10th ms columns
      yLevelStart = yLevel;
      for trial_i = 1:trialCount
        yLevel = yLevel + 1;
        trialSpikes = spikesByItem{item_i}{channel_i}{unit_i}(trial_i);
        for spike_i = 1:length(trialSpikes.times)
          spikeImgInd = round(trialSpikes.times(spike_i) * 10) + spikeImgIndOffset;
          if spikeImgInd > 0 && spikeImgInd < (traceTime * 10)
            tmp(trial_i, spikeImgInd-spiDi:spikeImgInd+spiDi, :) = deal(0);
          end
        end
      end
      spikeImg = [spikeImg; tmp];
    end
    spikeImgH = image(spikeImg);
    spikeImgH.AlphaData = ~squeeze(spikeImg(:,:,1));
    spikeImgH.XData = [-preAlign postAlign + imDur];
    ylim([0,size(spikeImg,1)+0.5]);
    
  case 3
    [spikeX, spikeY, spikeC] = deal([]);
    % Traditional representation of spikes as lines.
    for item_i = 1:length(spikesByItem)
      yLevelStart = yLevel;
      trialLabels = [trialLabels 1:length(spikesByItem{item_i}{channel_i}{unit_i})];
      for trial_i = 1:length(spikesByItem{item_i}{channel_i}{unit_i})
        yLevel = yLevel + 1;
        trialSpikes = spikesByItem{item_i}{channel_i}{unit_i}(trial_i);
        spikeX = [spikeX; trialSpikes.times];
        spikeY = [spikeY; repmat(yLevel, [length(trialSpikes.times), 1])];
        if colorSwitch == 1
          spikeC = [spikeC; repmat(spikeColor{item_i}{trial_i}.color, [length(trialSpikes.times), 1])];
        else
          spikeC = [spikeC; repmat([0 0 0], [length(trialSpikes.times), 1])];
        end
      end
      
      %Plot Stim start, stim end, and Bottom dividing bar
      h = plot([0 0],[yLevelStart+0.5 yLevel+0.5],'color', 'black', 'LineWidth', 2); % Stimulus Onset
      plot([imDur imDur],[yLevelStart+0.5 yLevel+0.5],'color', 'black', 'LineWidth', 2); % Stimulus Offset
      if item_i ~= length(spikesByItem)
        plot([xlim()],[yLevel+0.5 yLevel+0.5],'color','black', 'LineWidth', 2, 'LineStyle', '--'); %Stimuli dividing line
      end
      legendHandles = vertcat(legendHandles,h);
      
    end
    
    x = scatter(rasterAxes, spikeX, spikeY, 8, spikeC, 'filled', '|', 'MarkerEdgeColor', 'black');

    
end

if isa(inputHandle, 'matlab.ui.Figure')
  title(inputHandle.Name)
end

xlimits = xlim();
ylim([0,yLevel+0.5]);
if xlimits(2) < imDur + ISI
  plot([imDur+ISI imDur+ISI],[0 yLevel],'b-');
end
trialLabelsTmp = cell(1, length(trialLabels));
for tick_i = 1:length(trialLabels)
  trialLabelsTmp{tick_i} = num2str(trialLabels(tick_i));
end

set(gca,'YTickLabels',trialLabelsTmp);
set(gca,'TickLength',[0 0]);
yticks((1:length(trialLabelsTmp))-0.5)
xlabel('Time after stimulus onset (ms)');
ylabel('single trials');
legend(legendHandles, pictureLabels, 'Location', 'northeastoutside');
hold off;

%Create a legend which gives the color code for the shaded region or the
%individual spikes. Invisible axes with bar plots for each color.
if colorSwitch == 2
  invisAx = axes('Color','none','XColor','none');
  
  handleToThisBarSeries = gobjects(length(attendedObjData.objList),1);
  for b = 1 : length(attendedObjData.objList)
    handleToThisBarSeries(b) = bar(nan, nan, 'BarWidth', 0.5);
    set(handleToThisBarSeries(b), 'FaceColor', attendedObjData.colorCode{b});
    hold on;
  end
  legend(attendedObjData.objList, 'Location', 'southeastoutside');
  invisAx.Visible = 'off';
  linkprop([rasterAxes invisAx],'Position');
end

end


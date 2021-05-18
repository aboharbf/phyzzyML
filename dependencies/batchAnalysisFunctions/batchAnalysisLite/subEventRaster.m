function subEventRaster()
% Function which creates a raster for every stimulus presented to a unit,
% and shades the region in the back w/ subEvent data.

% Unit to use - 20201206Mo 001 Ch43, U2

% Collect raster info
eventDataFile = 'C:\Users\aboha\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories\eventData.mat';
rasterDir = 'D:\DataAnalysis\20201206Mo\Basic\001';
rasterFile = fullfile(rasterDir, 'preprocessedData.mat');
paramFile = fullfile(rasterDir, 'AnalysisParams.mat');
load(eventDataFile); % Puts eventData into the workspace.

% Create the image suitable to these rasters.
load(rasterFile, 'spikesByEvent', 'eventIDs');
load(paramFile, 'psthParams');
chan_i = 43;
unit_i = 3; % Unit 2 = unit index 3.
ISI = 500;

ITI = psthParams.ITI - 100;
pre = psthParams.psthPre - ITI;
ImDur = psthParams.psthImDur;
post = psthParams.psthPost;
periodWindows = [pre ImDur post];
trialCount = cellfun(@(x) length(x{1}{1}), spikesByEvent);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the data from the table, use generateEventImage to make into
% matrix.
eventData = eventData(eventIDs,:);
eventList = eventData.Properties.VariableNames(1:4);
eventDataStim = eventData(:, eventList);
eventList = strrep(eventList, '_', ' ');
eventMat = generateEventImage(eventDataStim, 2800);

% Resort data and labels
stimInTable = eventData.Properties.RowNames;
stimInTable = extractBefore(stimInTable, '.avi');

nonMonkeyInd = ~contains(stimInTable, 'monkey');
goalDirInd = contains(stimInTable, 'GoalDir');
idleInd = contains(stimInTable, 'Idle');
socInd = ~goalDirInd & ~idleInd & ~nonMonkeyInd;

resortInd = [find(socInd); find(goalDirInd); find(idleInd); find(nonMonkeyInd)];
eventIDs = eventIDs(resortInd);
spikesByEvent = spikesByEvent(resortInd);
stimInTable = stimInTable(resortInd);
eventMat = eventMat(resortInd, :, :);
trialCount = trialCount(resortInd);

% Process the names
monkeyInd = contains(stimInTable, 'monkey');
stimInTablePlot = strcat(extractBetween(stimInTable(monkeyInd), 'monkey', '_'), cellfun(@(x) x(end), stimInTable(monkeyInd)));
otherNames = strcat(extractBefore(stimInTable(~monkeyInd), '_'), cellfun(@(x) x(end), stimInTable(~monkeyInd)));
stimInTablePlot = [stimInTablePlot; otherNames];

% Expand the eventMat to match the trials per stimulus on this data.
eventMatSize = size(eventMat);
eventMatStack = nan(eventMatSize(1) * length(eventList), eventMatSize(2));
stackingIndex = reshape(1:eventMatSize(1) * length(eventList), [eventMatSize(3) eventMatSize(1)]);

% Populate the larger matrix correctly
for ii = 1:size(eventMat,1)
  for jj = 1:4
    eventMatStack(stackingIndex(jj,ii), :) = squeeze(eventMat(ii, :, jj)) * jj;
  end
end

% Expand the eventMat so that layering can be done effectively. For every 1
% real row, make it into 10.
eventMatSize = size(eventMat);
newEventMat = nan(sum(trialCount), eventMatSize(2), eventMatSize(3));

for event_i = 1:eventMatSize(3)
  % initialize new layer
  eventMatNew = nan(sum(trialCount), eventMatSize(2));
  for row_i = 1:size(eventMat,1)
    startInd = find(isnan(eventMatNew(:,1)), 1);
    endInd = startInd + trialCount(row_i) - 1;
        
    eventMatNew(startInd:endInd,:) = repmat(eventMat(row_i,:,event_i), [trialCount(row_i), 1]);
    
    % Change the start ind for later events to allow overlap
    eventMatNew(startInd:startInd + event_i, :) = deal(0);
    
  end  
  newEventMat(:,:,event_i) = eventMatNew;
end

% Pad the event mat to fit on the raster
eventMatSize = size(newEventMat);
frontPad = nan(eventMatSize(1), ITI + pre, eventMatSize(3));
endPad = nan(eventMatSize(1), post, eventMatSize(3));

eventMatTmp = cat(2, frontPad, newEventMat);
newEventMat = cat(2, eventMatTmp, endPad);

% Plot the raster
% attendedObjData.tracePlotData - {1:end} is stacked.
% Each cell needs to be the total time long.
figH = figure('units', 'normalized', 'position', [0.3026 0.0380 0.52 0.50]);

eventColor = [[0.9 0 0]; [0.6 0 0]; [0 0 0.6]; [0 0.6 0]];
ytickPos = 5:10:length(stimInTable)*10 - 5;

eyeDataStruct = struct();
rasterColorCoded(figH, spikesByEvent, eventIDs, psthParams, ISI, chan_i, unit_i, eyeDataStruct, 0)

axesRaster = findobj(figH.Children, 'Type', 'Axes');
% axesRaster.YTick = 1;
% axesRaster.YTickLabel = {''};

legendH = findobj(figH.Children, 'Type', 'Legend');
delete(legendH);

axesArray = gobjects(size(newEventMat,3), 1);
for ii = 1:size(newEventMat,3)
  % make image
  eventImg = binary2RGB(newEventMat(: , :, ii), eventColor(ii,:));
  axesArray(ii) = axes(figH, 'Position', axesRaster.Position);
  imgH = image(eventImg);
  
  % Change the alpha data so non-stimulus stuff is 100% transparent
  %   imgH.CData = imgH.CData * ii;
  imgH.AlphaData = newEventMat(: , :, ii);
  
%   if ii == 1
    axesArray(ii).Visible = 'off';
%   else
%     hold on
%     imgH.Parent.FontSize = 12;
%     %   xlim([0 2800]);
%     yticks(ytickPos);
%     yticklabels(strrep(stimInTable, '_', ' '));
%     xticks(1);
%     xticklabels('');
%     ylabel('Stimulus');
%   end
  
end

chi = figH.Children;
figH.Children = flipud(chi);
axesRaster.Visible = 'off';

% Create the right labels


% Fix the axes
figH.Children(end).Visible = 'on';
figH.Children(end).YTick = cumsum(trialCount) - 5;
figH.Children(end).YTickLabel = stimInTablePlot;
ylabel('Stimulus Name');
figH.Children(end).XTick = axesRaster.XTick + psthParams.psthPre;
figH.Children(end).XTickLabel = axesRaster.XTickLabel;
xlabel('Time (ms)');
title('Raster');
figH.Children(end).FontSize = 15;

% Add Text
for event_i = 1:length(eventList)
  textH2 = text(1, 1, eventList{event_i}, 'color', eventColor(event_i, :), 'FontWeight', 'bold', 'FontSize', 14, 'Units', 'normalized');
  yPoint = 0.22 - (0.06 * (event_i - 1));
  textH2.Position = [0.7 yPoint 0];
end

outputDir = 'C:\Users\aboha\OneDrive\Lab\Documents\Figures April 2020';
figStruct = struct();
figStruct.saveFig = 1;                
figStruct.closeFig = 0;               
figStruct.exportFig = 1;             
figStruct.saveFigData = 0;  

saveFigure(outputDir, sprintf('ExampleRasterSubEventOverLay'), [], figStruct, []);

end
function subEventPlotPerStim(eventDataFile)
% Function which turns the eventData table into an image, where each row is
% a stimulus, each column is a time bin, and each row contains 4 different
% lines demarking the 4 different events in the eventData.

eventDataFile = 'C:\Users\aboha\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories\eventData.mat';
outputDir = 'C:\Users\aboha\OneDrive\Lab\Documents\Figures April 2020';

figStruct = struct();
figStruct.saveFig = 1;                
figStruct.closeFig = 0;               
figStruct.exportFig = 1;             
figStruct.saveFigData = 0;            

load(eventDataFile); % Puts eventData into the workspace.

% Stim in table - Make the table just stimuli in the 'NaturalSocial'
% paradigm.
stimInTable = eventData.Properties.RowNames;
stimInTable = extractBefore(stimInTable, '.avi');
stimInTable = stimInTable(1:68);
stimInTable = stimInTable(~contains(stimInTable, 'Following') & ~contains(stimInTable, 'Foraging') & ~contains(stimInTable, 'Holding') & ~contains(stimInTable, 'human') & ~contains(stimInTable, 'landscape') & ~contains(stimInTable, '5'));
stimInTable = stimInTable(~contains(stimInTable, 'human') & ~contains(stimInTable, 'landscape') & ~contains(stimInTable, '5'));
stimInTable = stimInTable(~contains(stimInTable, 'monkeyIdle_1304'));
stimInTable = [stimInTable; 'monkeyIdle_1305'];

% Extract the data from the table, use generateEventImage to make into
% matrix.
eventDataStim = eventData(strcat(stimInTable,'.avi'), :);
eventList = eventDataStim.Properties.VariableNames(1:4);
eventDataStim = eventDataStim(:, eventList);
eventList = strrep(eventList, '_', ' ');
eventMat = generateEventImage(eventDataStim, 2800);

% Make plot Labels.
stimInTablePlot = strcat(extractBetween(stimInTable, 'monkey', '_'), cellfun(@(x) x(end), stimInTable));
stimInTable = strrep(stimInTable, '_', '');

% Resort data and labels
goalDirInd = contains(stimInTablePlot, 'GoalDir');
idleInd = contains(stimInTablePlot, 'Idle');
socInd = ~goalDirInd & ~idleInd;

resortInd = [find(socInd); find(goalDirInd); find(idleInd)];
stimInTablePlot = stimInTablePlot(resortInd);
eventMat = eventMat(resortInd, :, :);

eventMatSize = size(eventMat);

eventMatStack = nan(eventMatSize(1)*4, eventMatSize(2));
stackingIndex = reshape(1:eventMatSize(1)*4, [eventMatSize(3) eventMatSize(1)]);

% Populate the larger matrix correctly
for ii = 1:size(eventMat,1)
  for jj = 1:4
    eventMatStack(stackingIndex(jj,ii), :) = squeeze(eventMat(ii, :, jj)) * jj;
  end
end

% Expand the eventMat so that layering can be done effectively. For every 1
% real row, make it into 10.
eventMatSize = size(eventMat);
newEventMat = nan(eventMatSize(1)*10, eventMatSize(2), eventMatSize(3));
for event_i = 1:eventMatSize(3)
  % initialize new layer
  eventMatNew = nan(eventMatSize(1) * 10, eventMatSize(2));
  for row_i = 1:size(eventMat,1)
    rowEnd = row_i * 10;
    rowStart = rowEnd-(9-event_i);
    rowCount = length(rowStart:rowEnd);
    eventMatNew(rowStart:rowEnd,:) = repmat(eventMat(row_i,:,event_i), [rowCount, 1]);
  end  
  newEventMat(:,:,event_i) = eventMatNew;
end
eventMat = newEventMat;

figH = figure('units', 'normalized', 'position', [0.3026 0.0380 0.6427 0.5027]);

eventColor = [[0.9 0 0]; [0.6 0 0]; [0 0 0.6]; [0 0.6 0]];
ytickPos = 5:10:length(stimInTable)*10 - 5;

for ii = 1:size(eventMat,3)
  % make image
  eventImg = binary2RGB(eventMat(: , :, ii), eventColor(ii,:));
  axesH = axes(figH);
  imgH = image(eventImg);
  
  % Change the alpha data so non-stimulus stuff is 100% transparent
  %   imgH.CData = imgH.CData * ii;
  imgH.AlphaData = eventMat(: , :, ii);
  
  if ii == 1
    hold on
    imgH.Parent.FontSize = 12;
    xlim([0 2800]);
    yticks(ytickPos);
    yticklabels(stimInTablePlot);
    xlabel('Time (ms)');
    ylabel('Stimulus');
  else
    axesH.Visible = 'off';
    
  end
  
end


% Add Text
for event_i = 1:length(eventList)
  textH2 = text(1, 1, eventList{event_i}, 'color', eventColor(event_i, :), 'FontWeight', 'bold', 'FontSize', 14, 'Units', 'normalized');
  yPoint = 0.75 - (0.05 * (event_i - 1));
  textH2.Position = [0.77 yPoint 0];
end

saveFigure(outputDir, sprintf('Sub event occurances in Stimuli'), [], figStruct, []);

end
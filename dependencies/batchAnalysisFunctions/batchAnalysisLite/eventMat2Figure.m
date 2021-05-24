function eventMat2Figure(eventMat, stimInTablePlot, eventList, figTitle)
% Turns a 3 dimensional matrix into a figure
% - Input:
% eventMat: a stim*bin*event matrix
% stimInTablePlot: the labels for the stim dimension of eventMat.

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
ytickPos = 5:10:length(stimInTablePlot)*10 - 5;

for ii = 1:size(eventMat,3)
  % make image
  eventImg = binary2RGB(eventMat(: , :, ii), eventColor(ii,:));
  if ii ~= 1
    axesH = axes(figH, 'position', refAx.Position);
  else
    axesH = axes(figH);
  end
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
    title(figTitle);
    refAx = axesH;
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


end
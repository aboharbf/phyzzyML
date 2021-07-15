function plotCountsOnBrain(selTable, params)
% a function which plots a counts grid onto an image of the brain, with
% labeled axes.

[selIndGridMat, labelArray] = selInd2GridMat(selTable);

% For every label presented, Create an image with the axes relabeled
% according to the proper dimensions.

MLlabel = labelArray{1};
APlabel = labelArray{2};
eventLabels = strrep(labelArray{3}, '_', ' ');
unitType = labelArray{4};

positionC = find(strcmp(MLlabel, 'C'));
position8 = find(strcmp(APlabel, '8'));

if all(contains(selTable.dateSubj, 'Sam'))
  MLcoords = [4 5 6];
  APcoords = [25.5 26.5 27.5 28.5 29.5];
  [X, Y] = meshgrid(1:length(APcoords), 1:length(MLcoords));
else
  MLcoords = [4 5 6];
  APcoords = [25.5 26.5 27.5 28.5];
  [X, Y] = meshgrid(1:length(APcoords), 1:length(MLcoords));
end

cmap = colormap();
xCoord = reshape(X, [], 1);
yCoord = reshape(Y, [], 1);

for unit_i = 1:length(unitType)
  for event_i = 1:length(eventLabels)
    
    eventData = squeeze(selIndGridMat(:,:, event_i, unit_i));
    eventName = eventLabels{event_i};
    
    % Collect the data
    eventData = reshape(eventData, [], 1);
    
    if ~all(eventData == 0)
      % See which points to keep
      keepInd = eventData ~= 0;
      xCoordunit = xCoord(keepInd);
      yCoordunit = yCoord(keepInd);
      eventData = eventData(keepInd);
      
      % Scale to color line, and create dummy objects to create the desired
      % colorbar.
      figure();
      dummyImageHandle = imagesc(1:max(eventData));
      dummyImageHandle.Parent.Visible = 'off';
      colorbarH = colorbar();
      colorbarH.Label.String = 'Count';
      
      % Put Axes on top.
      axesH = axes();
      linkprop([dummyImageHandle.Parent, axesH], 'Position');
      
      % Create an index for referencing the color map.
      scaleFactor = max(eventData)/size(cmap,1);      % Scale the highest value in the data to the index for the highest value in the color map.
      sizeDataScaled = round(eventData/scaleFactor);
            
      % Plot them
      scatter(xCoordunit, yCoordunit, sizeDataScaled, cmap(sizeDataScaled, :), 'filled');
      
      % Labels, colors
      ylabel('ML Position (mm)');
      axesH.YLim = [0.5 length(MLcoords)+0.5];
      axesH.YTick = 1:length(MLcoords);
      axesH.YTickLabel = MLcoords;
      set(axesH, 'YDir','reverse')
      
      xlabel('AP Position (mm)');
      axesH.XLim = [0.5 length(APcoords)+0.5];
      axesH.XTick = 1:length(APcoords);
      axesH.XTickLabel = APcoords;
      
      axesH.Color = [0.5 0.5 0.5];
      axesH.FontSize = 11;
      
      figTitle = sprintf('%s - %s', eventName, unitType{unit_i});
      title(figTitle);
      saveFigure(fullfile(params.outputDir, 'plotsOnBrains'), figTitle, [], params.figStruct, [])
      
    end
    
  end
end


end
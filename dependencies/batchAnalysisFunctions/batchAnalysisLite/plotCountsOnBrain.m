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

MLcoords = [4 5 6 7];
APcoords = [25.5 26.5 27.5 28.5];
[X, Y] = meshgrid(1:length(MLcoords), 1:length(APcoords));

xCoord = reshape(X', [length(X(:)), 1]);
yCoord = reshape(Y', [length(X(:)), 1]);

cmap = colormap();

for unit_i = 1:length(unitType)
  for event_i = 1:length(eventLabels)
    
    eventData = squeeze(selIndGridMat(:,:, event_i, unit_i));
    eventName = eventLabels{event_i};
    figTitle = sprintf('%s - %s', eventName, unitType{unit_i});
    
    % Collect the data
    eventData = reshape(eventData', [length(xCoord),1]);
    
    if ~all(eventData == 0)
      
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
      scaleFactor = max(eventData)/256;
      sizeDataScaled = round(eventData/scaleFactor);
      
      % See which points to keep
      keepInd = sizeDataScaled ~= 0;
      xCoordunit = xCoord(keepInd);
      yCoordunit = yCoord(keepInd);
      sizeDataScaled = sizeDataScaled(keepInd);
      
      % Plot them
      scatter(xCoordunit, yCoordunit, sizeDataScaled, cmap(sizeDataScaled), 'filled');
      
      % Labels, colors
      ylabel('ML Position (mm)');
      axesH.YLim = [0.5 4.5];
      axesH.YTick = 1:4;
      axesH.YTickLabel = MLcoords;
      set(axesH, 'YDir','reverse')
      
      xlabel('AP Position (mm)');
      axesH.XLim = [0.5 4.5];
      axesH.XTick = 1:4;
      axesH.XTickLabel = APcoords;
      
      axesH.Color = [0.5 0.5 0.5];
      axesH.FontSize = 11;
      
      title(figTitle);
      saveFigure(fullfile(params.outputDir, 'plotsOnBrains'), figTitle, [], params.figStruct, [])
      
    end
    
  end
end


end
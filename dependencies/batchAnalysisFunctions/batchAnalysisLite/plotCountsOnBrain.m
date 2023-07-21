function plotCountsOnBrain(selTable, params)
% a function which plots a counts grid onto an image of the brain, with
% labeled axes.

selTable = selTable(contains(selTable.unitType, digitsPattern), :);

% For every label presented, Create an image with the axes relabeled
% according to the proper dimensions.

monkeyVec = {'Mo', 'Sam'};
[~, labelArray] = selInd2GridMat(selTable, params);
eventLabels = strrep(labelArray{3}, '_', ' ');
[eventDataAll, shapeDataAll, xCoordunitAll, yCoordunitAll, MLcoordsAll, APcoordsAll] = deal(cell(size(eventLabels)));

for m_i = 1:length(monkeyVec)
  
  % filter the selTable
  selTableMonkey = selTable(contains(selTable.dateSubj, monkeyVec{m_i}),:);
  [selIndGridMat, labelArray] = selInd2GridMat(selTableMonkey, params);
  
  MLlabel = labelArray{1};
  APlabel = labelArray{2};
  eventLabels = strrep(labelArray{3}, '_', ' ');
  unitType = labelArray{4};
  
  %if strcmp(params.monkeyTag, 'Sam')
  switch monkeyVec{m_i}
    case 'Sam'
    positionA = find(strcmp(MLlabel, 'A'));
    position3 = find(strcmp(APlabel, '3'));
    
    coordsA = 5;
    coords3 = 30;
    MLcoords = coordsA - (positionA - 1): coordsA - (positionA - 1) + length(MLlabel) - 1;
    APcoords = coords3 - (position3 - 1): coords3 - (position3 - 1) + length(APlabel) - 1;
    
  %elseif strcmp(params.monkeyTag, 'Mo')
    case 'Mo'
    
    positionC = find(strcmp(MLlabel, 'C'));
    position8 = find(strcmp(APlabel, '8'));
    
    coordsC = 5;
    coords8 = 28;
    MLcoords = coordsC - (positionC - 1): coordsC - (positionC - 1) + length(MLlabel) - 1;
    APcoords = coords8 - (position8 - 1): coords8 - (position8 - 1) + length(APlabel) - 1;
    
  end
      
  [X, Y] = meshgrid(APcoords, MLcoords);
  xCoord = reshape(X, [], 1);
  yCoord = reshape(Y, [], 1);
  
  for unit_i = 1:length(unitType)
    for event_i = 1:length(eventLabels)
      
      MLcoordsAll{event_i} = [MLcoordsAll{event_i}, MLcoords];
      APcoordsAll{event_i} = [APcoordsAll{event_i}, APcoords];
      
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
      end
      
      switch monkeyVec{m_i}
        case 'Sam'
          xCoordunit = xCoordunit + 0.05;
          yCoordunit = yCoordunit + 0.05;
          shapeDataAll{event_i} = [shapeDataAll{event_i}; repmat("s", size(xCoordunit))];
        case 'Mo'
          xCoordunit = xCoordunit - 0.05;
          yCoordunit = yCoordunit - 0.05;
          shapeDataAll{event_i} = [shapeDataAll{event_i}; repmat("o", size(xCoordunit))];
      end
      
      % Store in the largerMat
      xCoordunitAll{event_i} = [xCoordunitAll{event_i}; xCoordunit];
      yCoordunitAll{event_i} = [yCoordunitAll{event_i}; yCoordunit];
      eventDataAll{event_i} = [eventDataAll{event_i}; eventData];
      
    end
  end
end

% Scale to color line, and create dummy objects to create the desired
% colorbar.
for event_i = length(eventLabels) %[1:length(eventLabels)-2 length(eventLabels)]
  figTitle = sprintf('Units Per Recording Site - %s - %s', eventLabels{event_i}, unitType{unit_i});
  figH = figure('Name', figTitle, 'units', 'normalized', 'position', [0.5521 0.3972 0.3646 0.4861]);
  colormap hot
  cmap = colormap();

  colormap hot
  dummyImageHandle = imagesc(1:max(eventDataAll{event_i}));
  dummyImageHandle.Parent.Visible = 'off';
  colorbarH = colorbar();
  colorbarH.Label.String = 'Count';
  
  % Put Axes on top.
  axesH = axes();
  linkprop([dummyImageHandle.Parent, axesH], 'Position');
  
  % Create an index for referencing the color map.
  scaleFactor = max(eventDataAll{event_i})/size(cmap,1);      % Scale the highest value in the data to the index for the highest value in the color map.
  sizeDataScaled = round(eventDataAll{event_i}/scaleFactor);
  
  % If you want to simplify the plot, make all the shapes the same size
  colorData = sizeDataScaled;
  if 1
    sizeDataScaled = ones(size(sizeDataScaled)) * max(sizeDataScaled);
  end
  
  shapeDataUnique = unique(shapeDataAll{event_i});
  % Plot them
  hold on
  colormap hot
  for shape_i = 1:length(shapeDataUnique)
    plotInd = strcmp(shapeDataAll{event_i}, shapeDataUnique(shape_i));
    if any(sizeDataScaled(plotInd))
      scatterH = scatter(xCoordunitAll{event_i}(plotInd), yCoordunitAll{event_i}(plotInd), sizeDataScaled(plotInd), cmap(colorData(plotInd), :), 'filled');
      scatterH.Marker = shapeDataUnique(shape_i);
    end
  end
  
  % Labels, colors
  ylabel('ML Position (mm)');
  set(axesH, 'YDir','reverse')
  axesH.YLim = [3.5, 7.5];
  axesH.YTick = 4:7;
  
  xlabel('AP Position (mm)');
  set(axesH, 'XDir','reverse')
  axesH.XLim = [26.5 32.5];
  axesH.XTick = 27:32;
  
  axesH.Color = [1 0 0];
  axesH.FontSize = 12;
  colorbarH.FontSize = 12;
  legend({'Monkey 1', 'Monkey 2'}, 'location', 'northwest', 'FontSize', 20)
  title(figTitle);

  saveFigure(fullfile(params.outputDir, unitType{unit_i}, 'plotsOnBrains'), figTitle, [], params.figStruct, [])
end

end
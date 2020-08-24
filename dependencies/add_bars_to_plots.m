function [newBarAxes, barDummyHands] = add_bars_to_plots(targAx, targFig, dataMat, colorMat, show_ind, yLims)
% Function intended to add a bar to plot, typically with the focus on a
% segmented bar highlighting the time of events in a continous data trace.
%
% Inputs:
% - targAx = target axes for the new bar. (default if [] = gca)
% - dataMat = bar * bins.
% - colorMat = bar * 1 cell array with RGB values in each cell.
% - show_ind = bar * bin logical array, showing where the bar should be
% visible. (default if [] = all visible).
% - yVals = the y values [top, bottom] in the current axes which correspond
% to edges of bar. (default = bottom 1/10th of the plot).
% Output:
% - barImgHand = a handle to the newly created bar.
% - barDummyHands = dummy handles for lines of the same color for adding to
% the legend of a plot.

assert(size(dataMat,1) <= length(colorMat), 'Insufficient colors defined in colorMat')

if isempty(targAx)
  targAx = gca;
end

if isempty(targFig)
  targFig = gcf;
end

if isempty(show_ind)
  show_ind = true(size(dataMat));
end

if isempty(yLims)
  % Default is bottom 10% of Axes. Need Bottom and Height/10.
  yLims = [targAx.Position(2), targAx.Position(4)/10];
end

axes2Add = cellfun(@(x) any(any(x)), show_ind);

% use yLims to determine where new axes will go
heightPerBar = yLims(1)/sum(axes2Add);
bottomPerBar = yLims(1):heightPerBar:(yLims(1) + yLims(2));
bottomPerBar = fliplr(bottomPerBar);
% yDataForImg = [yLims(1) + yRangeSize/4, yLims(2) - yRangeSize/4];

% for each line, generate an axis to holds the new image, and gives it a
% distinct color map.
[newBarAxes, newBarImg] = deal(gobjects(sum(axes2Add)));
colorBarLocations = {'westoutside', 'west', 'eastoutside', 'east'};

for img_i = 1:length(colorMat)
  % Check if there are any events in the catagory to represent.
  if axes2Add(img_i)
    newBarAxes(img_i) = axes();   % Generate an axis
    newBarColorMap = colorGradient([1 1 1], colorMat{img_i}, 256);
    
    newBarImg(img_i) = imagesc(dataMat{img_i});     % Add the image data
    newBarImg(img_i).AlphaData = show_ind{img_i};   % Blank out empty regions.
    colormap(newBarAxes(img_i), newBarColorMap);    % Assign the correct colormap
    
    if length(unique(unique(dataMat{img_i}))) > 2
      colorbar(colorBarLocations{img_i});             % Add colorbar, if needed
    end
    
    newBarAxes(img_i).Position = targAx.Position;         % Match the axes of the larger figure
    newBarAxes(img_i).Position(2) = bottomPerBar(img_i);  % Position it correctly, bottom
    newBarAxes(img_i).Position(4) = heightPerBar;         % Position it correctly, height
    newBarAxes(img_i).Visible = 'off';                    % Turn off the axes
  end
end

% Add dummy lines matching the desired colors and return handles, for the
% sake of complete legends
dummyLineYVals = linspace(yLims(1), yLims(2), 2+size(dataMat{1},1));
dummyLineYVals = dummyLineYVals(2:end-1);

if length(dummyLineYVals) == 1
  % Plot a single point at the first instance of each line, assuming it has
  % any entries, and add it to the handles
  barDummyHands = plot(targAx, 0, dummyLineYVals, 'color', colorMat{1}, 'linewidth', 15);
else
  barDummyHands = plot(targAx, 0, dummyLineYVals, 'color', 'k', 'linewidth', 1);
end

set(gcf, 'CurrentAxes', targAx);

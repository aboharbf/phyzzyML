function [barImgHand, barDummyHands] = add_bars_to_plots(targAx, dataMat, colorMat, show_ind, yLims)
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

if isempty(show_ind)
  show_ind = true(size(dataMat));
end

if isempty(yLims)
  yLims = [targAx.YLim(1), targAx.YLim(2)/10];
end

% use yVals to generate plot
yRangeSize = yLims(2) - yLims(1);
yDataForImg = [yLims(1) + yRangeSize/4, yLims(2) - yRangeSize/4];

% Create CData matrix of correct size
CDataForPlot = zeros([size(dataMat), 3]);
for ii = 1:size(dataMat,1)
  colorSlice(1, 1, :) = colorMat{ii};
  CDataForPlot(ii, :, :) = repmat(colorSlice, [1, size(dataMat,2), 1]);
end

% Flip the data so the data/colors defined first are on top
dataMat = flipud(dataMat);
CDataForPlot = flipud(CDataForPlot);
show_ind = flipud(show_ind);

barImgHand = imagesc([0 size(dataMat,2)], yDataForImg, CDataForPlot);
barImgHand.AlphaData = show_ind;

% Add dummy lines matching the desired colors and return handles, for the
% sake of complete legends
dummyLineYVals = linspace(yDataForImg(1), yDataForImg(2), 2+size(dataMat,1));
dummyLineYVals = dummyLineYVals(2:end-1);
barDummyHands = gobjects(size(dummyLineYVals));
for ii = 1:length(gobjects)
  % Plot a single point at the first instance of each line, assuming it has
  % any entries, and add it to the handles
  firstBin = find(dataMat(ii,:),1);
  barDummyHands(ii) = plot(firstBin, dummyLineYVals, 'color', colorMat{ii}, 'linewidth', 15);
end
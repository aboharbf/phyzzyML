function [pointsForLabels, binsForLabels, binLabelsNew, colorPoints] = bin2Label(binLabelInfo, binsForLabels, binLabels, segmentColors)
% A function which takes in desired bins to be labeled, along with the
% array of the present data, and returns the correct bins to label and
% where.
% Inputs:
% - binsPresent: an array of n*1 points showing what the labels are for the
% bins present
% - binsForLabels: an array showing the value of the bins which need to be
% labeled, or where those labels should be.
% - binLabels: the labels of the bins. Fed in incase the corresponding bin
% isn't present and the label must be deleted.

binsPresent = binLabelInfo.the_bin_start_times_shift;

startSegColor = [0.4 0.4 0.4];
endSegColor = [0.7 0.7 0.7];

% Make sure binsForLabels doesn't exceed binsPresent on either end
bins2Pick = (binsPresent(1) <= binsForLabels) & (binsPresent(end) >= binsForLabels);
binLabels = binLabels(bins2Pick);
binsForLabels = binsForLabels(bins2Pick);
if bins2Pick(end) == 0
  segmentColors = segmentColors(1:end-1, :);
elseif bins2Pick(1) == 0
  segmentColors = segmentColors(2:end, :);
end

% If the ends are present or not, add them on
if ~ismember(binsPresent(end), binsForLabels)
  binsForLabels = [binsForLabels binsPresent(end)];
  binLabels = [binLabels 'End'];
  segmentColors = [segmentColors; endSegColor];
end

if ~ismember(binsPresent(1), binsForLabels)
  binsForLabels = [binsPresent(1) binsForLabels];
  binLabels = ['Start' binLabels];
  segmentColors = [startSegColor; segmentColors];
end

pointsForLabels = round(interp1(binsPresent, 1:length(binsPresent), binsForLabels));
allPoints = round(interp1(binsPresent, 1:length(binsPresent), binsForLabels));
binLabelsNew = binLabels;

% Create the segment colors
colorPoints = nan(length(binsPresent), 3);
for ii = 1:size(segmentColors, 1)
  inds2Add = pointsForLabels(ii):pointsForLabels(ii+1);
  colorPoints(inds2Add, :) = repmat(segmentColors(ii,:), [length(inds2Add), 1]);
end

% Package outputs correctly
xAxisLabelStruct.pointsForLabels = pointsForLabels;
xAxisLabelStruct.binsForLabels = binsForLabels;
xAxisLabelStruct.labels = binLabelsNew;
xAxisLabelStruct.colorPerBin = colorPoints;

end
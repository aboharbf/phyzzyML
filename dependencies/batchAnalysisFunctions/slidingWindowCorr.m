function [corrDataMean, binStarts, binEnds] = slidingWindowCorr(eyeMat, binSize, stepSize)
% Performs sliding window correlation across the data present, starting at
% the beginning, going in steps.

binCount = size(eyeMat,2);

binStarts = 1:stepSize:binCount - binSize;
binEnds = binStarts + binSize;

corrData = nan([size(eyeMat,1), size(eyeMat,1),length(binStarts)]);
for ii = 1:length(binStarts)
  corrData(:, :, ii) = corr(eyeMat(:, binStarts(ii):binEnds(ii))', 'rows', 'pairwise');
end

corrDataMean = squeeze(mean(corrData, [1, 2]));
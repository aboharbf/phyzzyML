function [binnedData, binStarts, binEnds] = binStuff(data2Bin, binSizes, binStep)
% Function which bins previously binned events, by summing across the bins
% of a particular index.

% 7.22.2021 - probably going to run into an error, i changed the arguments,
% not sure where else it is called.

if isempty(binStep)
  binStep = binSizes;
end

totalBinsCurrent = size(data2Bin,2);
binStarts = 1:binStep:totalBinsCurrent;
binEnds = binStarts + binSizes;
binEdges = [binStarts; binEnds];

% Make sure the bins fit
removeBin = any(binEdges > totalBinsCurrent);
binEdges = binEdges(:, ~removeBin);

binnedData = nan(size(data2Bin,1), length(binEdges) - 1);
for ii = 1:length(binEdges)
  binnedData(:, ii) = sum(data2Bin(:, binEdges(1,ii):binEdges(2, ii)), 2);
end

binStarts = binEdges(1,:);
binEnds = binEdges(2,:) - 1;

end
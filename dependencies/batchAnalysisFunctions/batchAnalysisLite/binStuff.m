function binnedData = binStuff(data2Bin, binSizes)
% Function which bins previously binned events, by summing across the bins
% of a particular index.

totalBinsCurrent = size(data2Bin,2);
binEdges = [1:binSizes:totalBinsCurrent totalBinsCurrent];

binnedData = nan(size(data2Bin,1), length(binEdges) - 1);
for ii = 1:length(binEdges) - 1
  binnedData(:, ii) = sum(data2Bin(:, binEdges(ii):binEdges(ii+1)-1), 2);
end

end
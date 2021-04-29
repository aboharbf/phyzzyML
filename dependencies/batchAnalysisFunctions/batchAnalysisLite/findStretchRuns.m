function stretchesPerRow = findStretchRuns(pValInd2Keep)

zPad = zeros(size(pValInd2Keep,1),1);

for layer_i = 1:size(pValInd2Keep,3)
  % Pad the vector, diff
  pValStackDiff = diff([zPad, pValInd2Keep(:,:,layer_i), zPad], 1, 2); %  pad of zeros, diff
  
  % for each row, find the stretch length
  stretchLength = zeros(size(pValStackDiff,1),1);
  for row_i = 1:length(stretchLength)
    starts = find(pValStackDiff(row_i, :) >= 1);
    stops =  find(pValStackDiff(row_i, :) <= -1);
    lengths = stops - starts;
    
    if ~isempty(lengths)
      stretchLength(row_i) = max(lengths);
    end
    
  end
  
  stretchesPerRow(:, :, layer_i) = stretchLength;

end
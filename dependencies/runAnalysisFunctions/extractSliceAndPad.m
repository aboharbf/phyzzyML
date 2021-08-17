function slicedData = extractSliceAndPad(data, ind, pre, post, pad)
% A function which returns a slice of a larger array defined by the window
% Pre:Post in length, centered on Ind. Pad is used if the Pre or Post
% extend past the edge of the data provided.
% Inputs
% - Data: an n*bins Array of data, which Ind goes into. if n == 1 in data,
%   and not in Ind, the same data is used for all inds.
% - Ind: a 1*n set of indices.
% - Pre: The amount before each ind to be extracted
% - Post: The amount after each ind to be extracted.
% - Pad: What is added in the occassion that pre or post extend beyond
% edges of data.

switch pad
  case 'NaN'
    slicedData = nan(length(ind), length(-pre:post));
end

 for ind_i = 1:length(ind)
  indWin = ind(ind_i)-pre:ind(ind_i)+pre;
  
  % Check front
  if any(indWin <= 0)
    % Adjust
    startInd = indWin > 0;
    indWin = indWin(startInd);
    startInd = find(startInd, 1, 'first');
  else
    startInd = 1;
  end
  
 % Check Back
 if any(indWin > length(data))
      % Adjust
    endInd = indWin <= length(data);
    indWin = indWin(endInd);
    endInd = find(endInd, 1, 'last');
 else
   endInd = length(-pre:post);
 end
 
 % Insert adjusted slice
 slicedData(ind_i, startInd:endInd) = data(indWin);
 
end

end
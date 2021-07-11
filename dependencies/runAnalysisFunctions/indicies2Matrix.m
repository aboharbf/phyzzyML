function indexMat = indicies2Matrix(indexArray, shift, pre, dur, post)
% A function which converts an input of a cell array with indicies in each
% into a matrix
% Input:
% - indexArray: trial*1 cell array consisting of indicies for a cell array
% - shift: trial*1 vector denoting the shift added to the indicies (converts times to actual indices)
% - pre: the negative threshold for inclusion
% - post: the positive threshold for inclusion.

% Output:
% indexMat: a trial*(pre+post) zero array with a 1 where indices denote.

if length(shift) == 1
  shift = repmat(shift, [length(indexArray), 1]);
end

% Initialize the vector
indexMat = zeros(length(indexArray), pre+dur+post);

% Cycle through trials
for trial_i = 1:length(indexArray)
  if ~isempty(indexArray{trial_i})
    
    % Add the shift and the pre.
    saccTimes = indexArray{trial_i}(1,:) + shift(trial_i) + pre;
    
    % Filter based on 0 > and < pre + post
    saccKeep = saccTimes > 0 & saccTimes < (pre+dur+post);
    if any(saccKeep)
      % Add them to the matrix
      saccTimes = saccTimes(saccKeep);
      indexMat(trial_i, saccTimes) = deal(1);
    end
  end
end

end
function corrDist = splitHalvesCorrelation(eventTimes, dataTrace, loopCount, gaussFilt)
% A function which takes halves of the data iteratively and generates a
% distribution of correlations.

corrDist = nan(loopCount, 1);

% If odd number of events, remove 1 randomly.
if mod(length(eventTimes),2) ~= 0
  eventTimes = eventTimes(randperm(length(eventTimes), length(eventTimes)-1));
end
    
for split_i = 1:loopCount
  % Pick half
  trainInd = randperm(length(eventTimes), length(eventTimes)/2);
  trainIndMat = false(size(eventTimes));
  trainIndMat(trainInd) = true;
  
  % Check Correlations
  eyePosSpikeTrain = dataTrace(:,eventTimes(trainIndMat));
  eyePosSpikeTest = dataTrace(:,eventTimes(~trainIndMat));
  
  corrSets = cat(3, eyePosSpikeTrain, eyePosSpikeTest);
  activityGrid = zeros(55, 55, size(corrSets,3));
  
  for ii = 1:size(corrSets,3)
    eyePosSpikeInd = round(corrSets(:,:,ii) + 30);
    for jj = 1:size(corrSets, 2)
      activityGrid(eyePosSpikeInd(1,jj), eyePosSpikeInd(2,jj), ii) = activityGrid(eyePosSpikeInd(1,jj), eyePosSpikeInd(2,jj), ii) + 1;
    end
    activityGrid(:,:,ii) = convn(activityGrid(:,:,ii), gaussFilt, 'same');
  end
  
  % Vectorize and see correlation b/t split halves
  half1 = activityGrid(:,:,1);
  half2 = activityGrid(:,:,2);
  corrDist(split_i) = corr(half1(:), half2(:));
  
end

end
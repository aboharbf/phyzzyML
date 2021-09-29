function selTable = detectingRampingActivity(selTable, spikesByEventBinned, epochTargParams, psthParams)
% a function which detects ramping activity across epochs of trial.
% Does this by taking half the trials, finding the ramps in the epoch, and
% reporting whether they exist in the other half. Groups things
spikesBinned = deEventSpikeData(spikesByEventBinned);
trialCount = size(spikesBinned{1}{1},1);

% Windows to test
epochKeepInd = ~strcmp(epochTargParams.labels, 'preFix');
epochTimes = epochTargParams.times(epochKeepInd, :);
epochLabels = epochTargParams.labels(epochKeepInd);

binInc = 25;
minimumRampSize = 200;

% Remove the pre-fixation period contents, and post reward
startInd = epochTimes(strcmp(epochLabels, 'Fix'), 1) + psthParams.psthPre + psthParams.movingWin(1)/2;
endInd = epochTimes(strcmp(epochLabels, 'reward'), 2) + psthParams.psthPre + psthParams.movingWin(1)/2;
for c_i = 1:length(spikesBinned)
  for u_i = 1:length(spikesBinned{c_i})
    spikesBinned{c_i}{u_i} = spikesBinned{c_i}{u_i}(:, startInd:endInd);
  end
end

% Convert binStarts to times - Temporary: Don't record epoch of the ramp.
% Leave that for later.
[~, binStarts, binEnds] = binStuff(spikesBinned{1}{1}, binInc, binInc);
binStartsTimes = binStarts + epochTimes(strcmp(epochLabels, 'Fix'), 1);
binEndsTimes = binEnds + epochTimes(strcmp(epochLabels, 'Fix'), 1);

% Create a grid for sampling
binDist = -(binStarts' - binStarts);

% Create a vector to determine when the ramp begins, epoch wise.
epochBins = zeros(size(binStartsTimes));
for epoch_i = 1:size(epochTimes,1)
  epochBins(epochTimes(epoch_i,1) < binStartsTimes & epochTimes(epoch_i,2) > binStartsTimes) = epoch_i;
end

% Modify binEnds2Use to account for epoch differences we want to
% encorperate
epochMinDist = [400, 250, 1000, 200];
binEnds2Use = false(size(binDist));
for epoch_i = 1:length(epochMinDist)
  epochIndex = epochBins == epoch_i;
  binEnds2Use(epochIndex,:) = binDist(epochIndex,:) >= epochMinDist(epoch_i);
  
  % Remove bins which extend past the edge of the current epoch
  lastEpochInd = find(epochIndex,1,'last');
  if lastEpochInd ~= length(epochIndex)
    binEnds2Use(epochIndex,lastEpochInd+1:end) = deal(false);
  end
  
  % Make sure after the clipping at the end that all the bins still meet
  % the minimum distance requirement.
%   binLength = sum(binEnds2Use(epochIndex, :), 2) * binInc;
  
end

% figure(); imagesc(binEnds2Use); 
% xticks(1:length(binStartsTimes)); 
% xticklabels(binStartsTimes(xticks)); 
% yticks(1:length(binStartsTimes)); 
% yticklabels(binStartsTimes(yticks));

% Find bin start and ends to retry test
[X, Y] = meshgrid(1:length(binStarts), 1:length(binStarts));

% Visualize windows being sampled
if 0
  binIndsTmp = find(binEnds2Use);
  arrayTmp2 = false(sum(binEnds2Use(:)), size(binEnds2Use,1));
  for ii = 1:length(binIndsTmp)
    arrayTmp2(ii, Y(binIndsTmp(ii)):X(binIndsTmp(ii))) = true;
  end
  figure(); imagesc(arrayTmp2);
end

rampsPresent = cell(size(selTable,1),1);
unitInd = 1;
for chan_i = 1:length(spikesBinned)
  for unit_i = 1:length(spikesBinned{chan_i}) 
    
    % Skip unsorted activity.
    if unit_i == 1
      unitInd = unitInd + 1;
      continue
    end
        
    % Extract data, and increase bin size, non-overlapping.
    unitRasters = binStuff(spikesBinned{chan_i}{unit_i}, binInc, binInc);
    
    testingHalfInd = randperm(trialCount, round(trialCount/2));
    trainInd = false(trialCount,1);
    trainInd(testingHalfInd) = true;
    
    % take means of each half
    meanPSTHtrain = mean(unitRasters(trainInd,:), 1);
    meanPSTHtest = mean(unitRasters(~trainInd,:), 1);
    
    % Test Ramping
    SSTwholeTrace = sum((meanPSTHtrain-mean(meanPSTHtrain)).^2);
    rampingCheckResult = nan([size(binEnds2Use), 2]);

    for bin_s = 1:length(binStarts)
      binEnds2UseTmp = find(binEnds2Use(bin_s,:));
      for bin_e = binEnds2UseTmp
        % Run the regression
        dataSlice = meanPSTHtrain(bin_s:bin_e);
        lmFit = fitlm(1:length(dataSlice), dataSlice);
        
        % Save params
        if lmFit.Coefficients.pValue(2) <= 0.01
          rampingCheckResult(bin_s, bin_e, 1) = lmFit.SSR/SSTwholeTrace;
          rampingCheckResult(bin_s, bin_e, 2) = lmFit.Coefficients.pValue(2);
        end
       
      end
    end
    
    % Pick the best candidates based on R^2, sort them.
    priorityValues = rampingCheckResult(:, :, 1);
    [sortedVals, sortedValsInd] = sort(priorityValues(:), 'descend');
    keepInd = ~isnan(sortedVals);
    sortedValsInd = sortedValsInd(keepInd); % Remove NaNs
    sortedVals = sortedVals(keepInd);
           
    % Retrieve their bin indicies;
    binStartsTest = Y(sortedValsInd);
    binEndsTests = X(sortedValsInd);
    
    % Filter based on epoch/length
    binEdges = [binStartsTest, binEndsTests];
        
    % Make sure none to be investigated are overlapping
    elements2Keep = nonOverlappingWindows(binEdges, sortedVals);
    
    % Find the top 10
    elements2Keep = find(elements2Keep,10);
    binStartsTest = binStartsTest(elements2Keep);
    binEndsTests = binEndsTests(elements2Keep);
    sortedVals = sortedVals(elements2Keep);

    % Run the regression on the left out testing data
    rampingTestResults = nan(length(binStartsTest), 8);
    SSTotalTrace = sum((meanPSTHtest - mean(meanPSTHtest)).^2);
    for bin_i = 1:length(binStartsTest)
      dataSlice = meanPSTHtest(binStartsTest(bin_i):binEndsTests(bin_i));
      lmFit = fitlm(1:length(dataSlice), dataSlice);
      
      % Save params
      rampingTestResults(bin_i, 1) = epochBins(binStartsTest(bin_i));
      rampingTestResults(bin_i, 2) = binStartsTest(bin_i);
      rampingTestResults(bin_i, 3) = binEndsTests(bin_i);
      rampingTestResults(bin_i, 4) = lmFit.Coefficients.Estimate(2);
      rampingTestResults(bin_i, 5) = lmFit.Coefficients.pValue(2);
      rampingTestResults(bin_i, 6) = lmFit.SSR;                 % Model explained variance
      rampingTestResults(bin_i, 7) = lmFit.Rsquared.Adjusted;   % Explained fractional variance (Adjusted)
      rampingTestResults(bin_i, 8) = lmFit.SSR/SSTotalTrace;    % Explained total variance
      
    end
    
    % filter based on pVal
    rampingTestResults = rampingTestResults(rampingTestResults(:,5) <= 0.01, :);
    
    % If a ramp is detected, save the result
    if any(rampingTestResults)
      
      % Convert bins to times
      rampingTestResults(:,2:3) =  [binStartsTimes(:, rampingTestResults(:,2))', binEndsTimes(:, rampingTestResults(:,3))'];
      
      % store vector w/ epochInd, StartTime, Endtime, Slope, pVal, SSR, R^2
      % adjusted, and whole fractional variance explained.
      rampsPresent{unitInd} = rampingTestResults;
            
    end
    unitInd = unitInd + 1;
    
  end
end

selTable.rampStats = rampsPresent;

end
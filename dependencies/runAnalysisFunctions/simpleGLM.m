function simpleGLM(spikesByEventBinned, eyeDataStruct, movingWin, taskData, selTable)

% A function which dices and bins activity traces into 100 ms blocks, and
% applies labels to each bin according to
% - stimulusType: social and nonSocial.
% - stimulusCategory: chasing, fighting, etc
% - stimulusLabel: 'chasing1'
% - saccade: 0 for none, 1 for pre, 2 for duration.
% - stimulusEpoch: 'isi', 'fix', 'stimEarly', 'stimLate', 'Reward'
% - blink: 0, or 1 if this period is during a supposed blink.
% - Juice delivered: 0, or 1 if juice was just delivered.
% - rewardLastTrial = 0, or 1 if juice was delivered on the previous
% trial.

padLength = movingWin(1)/2;
trialPerStim = cellfun(@(x) size(x, 1), eyeDataStruct.saccadeByStim);

% Take spikesByEventBinned, find out how many 100 ms bins exist in the
% whole stimulus run, and reduce the bins to this number
totalBins = size(spikesByEventBinned{1}{1}{1}, 2) - movingWin(1);

newBinCount = floor(totalBins/100);
binStarts = 1:100:totalBins-100;
binEnds = [binStarts(2:end)-1 newBinCount*100];
timeBins = [binStarts', binEnds'];

% Bin spikes
spikesByEventBinnedLarge = cell(size(timeBins,1),1);
for stim_i = 1:length(spikesByEventBinned)
  for chan_i = 1:length(spikesByEventBinned{stim_i})
    for unit_i = 1:length(spikesByEventBinned{stim_i}{chan_i})
      unitAct = nan(trialPerStim(stim_i), newBinCount);
      for epoch_i = 1:size(timeBins,1)
        % Get Rid of Padding
        unitActivity = spikesByEventBinned{stim_i}{chan_i}{unit_i}(:, padLength + 1:end-padLength);
        
        unitAct(:, epoch_i) = sum(unitActivity(:, timeBins(epoch_i,1): timeBins(epoch_i,2)), 2);
        %spikeCountsByImageByEpoch{epoch}{channel}{unit}{event}.rates = trials*1
      end
      
      % Store
      spikesByEventBinnedLarge{stim_i}{unit_i}{unit_i} = unitAct;
      
    end
  end
end

% Identify bins based on eye behavior
eyeBehByEventBinnedLarge = cell(size(spikesByEventBinned, 1), 1);
for stim_i = 1:length(eyeDataStruct.saccadeByStim)
  % Initialize eye activity
  unitEye = nan(trialPerStim(stim_i), newBinCount);
  
  % Loop across epochs, binning eye activity by taking the mean.
  for epoch_i = 1:size(timeBins,1)
    unitEye(:, epoch_i) = round(mean(eyeDataStruct.saccadeByStim{stim_i}(:, timeBins(epoch_i,1): timeBins(epoch_i,2)), 2));
  end
  
  % Save to larger structure
  eyeBehByEventBinnedLarge{stim_i} = unitEye;
  
end

% Create the structures to feed into
rewardedTrials = ~isnan(taskData.juiceOnTimes);                                       % Reward of the previous trial
[taskEventIDsUnique, ~, trialsByEventNumber] = unique(taskData.taskEventIDsMerged);   % stimulusLabel
                                                                                      % stimulusCategory: chasing, fighting, etc
                                                                                      
% Run glm
% n = number of observations (a row of 0s and 1s representing which labels
% this bin has) - probably 46.
% p = predictor variables (a column of 0s and 1s representing a specific variable)

% y is an n*1 vector, each entry an observation (a spike count).


% X - Predictor variables, specified as an n-by-p numeric matrix, where n is the number of observations and p is the number of predictor variables. 
%     Each column of X represents one variable, and each row represents one observation.
% y - Response variable, specified as a vector or matrix. y must be an n-by-1 vector, where n is the number of observations. Each entry in y is 
%     the response for the corresponding row of X. The data type must be single or double.



end
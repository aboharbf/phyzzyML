function [selTable] = epochCompareStats(spikesByEvent, spikesByEventFixAlign, selTable, epochTargParams)
% This function performs a few statistical tests to determine whether there
% are meaningful difference between different periods of a task.

% Bin spikes
timeBins = epochTargParams.times;
timeLabels = epochTargParams.labels;

spikeCountsByImageByEpoch = cell(size(timeBins,1),1);
% First bin is found using fixAligned data
[spikeCountsByImageByEpoch{1}, ~, ~] = spikeCounter(spikesByEventFixAlign, timeBins(1, 1), timeBins(1, 2));

for epoch_i = 2:size(timeBins,1)
  [spikeCountsByImageByEpoch{epoch_i}, ~, ~] = spikeCounter(spikesByEvent, timeBins(epoch_i, 1), timeBins(epoch_i, 2));
end

%% Determine Epoch selectivity, compared to baseline
allUnitInd = 1;
pullRates = @(x) x.rates;
chanCount = length(spikesByEvent{1});
unitCounts = cellfun('length', spikesByEvent{1});
labelVsBaseline = timeLabels(2:end);

% Initialize
labelVsBasePVal = ones(sum(unitCounts), length(labelVsBaseline));
labelVsBaseMeans = zeros(sum(unitCounts), length(labelVsBaseline));

epochSelName = cell(sum(unitCounts), 1);
epochSelInd = zeros(sum(unitCounts), 1);

for chan_i = 1:chanCount
  for unit_i = 1:unitCounts(chan_i)
    
    % Check for fixation period selectivity -
    allBaseline = cellfun(pullRates, spikeCountsByImageByEpoch{strcmp(timeLabels, 'preFix')}{chan_i}{unit_i}, 'UniformOutput', false);
    allBaseline = vertcat(allBaseline{:});
    
    % Cycle through epochs which aren't baseline.
    epochValues = cell(length(labelVsBaseline),1);
    for epoch_i = 1:length(labelVsBaseline)
      
      allEpoch = cellfun(pullRates, spikeCountsByImageByEpoch{strcmp(timeLabels, labelVsBaseline{epoch_i})}{chan_i}{unit_i}, 'UniformOutput', false);
      allEpoch = vertcat(allEpoch{:});
      epochValues{epoch_i} = allEpoch;
      
      % Run the Test
      if epochTargParams.nonParametric
        [labelVsBasePVal(allUnitInd, epoch_i), ~, ~] = signrank(allEpoch, allBaseline);
      else
        [~, labelVsBasePVal(allUnitInd, epoch_i), ~] = ttest(allEpoch, allBaseline);
      end
      
      % save the difference
      labelVsBaseMeans(allUnitInd, epoch_i) = mean(allEpoch) - mean(allBaseline);
      
    end
    
    % Epoch selectivity index - Single comparison across epochs.
    tmpVals = labelVsBaseMeans(allUnitInd, :) + mean(allBaseline);
    [maxEp, maxEpI] = max(tmpVals);
    epochSelName(allUnitInd) = labelVsBaseline(maxEpI);
    tmpVals(maxEpI) = NaN;
    epochSelInd(allUnitInd) = (maxEp - nanmean(tmpVals))/(maxEp + nanmean(tmpVals));
    
    allUnitInd = allUnitInd + 1;
    
  end
end

% Add Epoch Selectivity Label and index
selTable.epochPref_Name = epochSelName;
selTable.epochPref_Ind = epochSelInd;

% Add Baseline comparisons to Sel table
for col_i = 1:length(labelVsBaseline)
  % Store p values
  tableVarNames = strcat('baseV_', labelVsBaseline{col_i}, '_pVal');
  selTable.(tableVarNames) = labelVsBasePVal(:, col_i);
  
  % Store mean differences b/t groups
  tableVarNames = strcat('baseV_', labelVsBaseline{col_i}, '_diff');
  selTable.(tableVarNames) = labelVsBaseMeans(:, col_i);
end

end

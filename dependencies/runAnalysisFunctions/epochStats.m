function [selTable] = epochStats(spikesByEvent, selTable, eventIDs, paradigm, epochStatsParams)
% This function performs a few statistical tests to determine whether there
% are meaningful difference between different periods of a task, and
% different conditions of the task.

% Use the simpler cell array of names to see which groups have been removed
% during preprocessing due to being absent.

% Bin spikes
timeBins = epochStatsParams.times;
timeLabels = epochStatsParams.labels;
spikeCountsByImageByEpoch = cell(size(timeBins,1),1);
for epoch_i = 1:size(timeBins,1)
  [spikeCounts, ~, ~] = spikeCounter(spikesByEvent, timeBins(epoch_i, 1), timeBins(epoch_i, 2));
  spikeCountsByImageByEpoch{epoch_i} = spikeCounts;
  %spikeCountsByImageByEpoch{epoch}{channel}{unit}{event}.rates = trials*1
end

%% Determine Epoch selectivity, compared to baseline
allUnitInd = 1;
pullRates = @(x) x.rates;
chanCount = length(spikesByEvent{1});
unitCounts = cellfun('length', spikesByEvent{1});
labelVsBaseline = timeLabels(2:end);
[labelVsBaseMeans, labelVsBasePVal] = deal(zeros(sum(unitCounts), length(labelVsBaseline)));
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
      if epochStatsParams.nonParametric
        [labelVsBasePVal(allUnitInd, epoch_i), ~, ~] = signrank(allEpoch, allBaseline);
      else
        [~, labelVsBasePVal(allUnitInd, epoch_i), ~] = ttest(allEpoch, allBaseline);
      end
      
      % If sig, save the difference
      % if ~isnan(pVal0) && pVal0 < alpha
      labelVsBaseMeans(allUnitInd, epoch_i) = mean(allEpoch) - mean(allBaseline);
      % end
      
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
selTable.epochPrefName = epochSelName;
selTable.epochPrefInd = epochSelInd;

% Add to Sel table
for col_i = 1:length(labelVsBaseline)
  % Store p values
  tableVarNames = strcat('BaseV', labelVsBaseline{col_i}, '_PVal');
  selTable.(tableVarNames) = labelVsBasePVal(:, col_i);
  
  % Store mean differences b/t groups
  tableVarNames = strcat('BaseV', labelVsBaseline{col_i}, '_Mean');
  selTable.(tableVarNames) = labelVsBaseMeans(:, col_i);
end

%% Determine Selectivity per Epoch selectivity, comparing target and non-target

targNames = epochStatsParams.(paradigm).targNames;
targLabelList = epochStatsParams.(paradigm).targ;
targetEpochsParadigm = epochStatsParams.(paradigm).targetEpochs;

labels = epochStatsParams.labels;
alpha = epochStatsParams.alpha;

% For plotIndex
plotParams.stimParamsFilename = epochStatsParams.stimParamsFilename;
plotParams.plotLabels = targLabelList;
plotMat = plotIndex(eventIDs, plotParams);

for group_i = 1:length(targLabelList)
  
  % Extract relevant variables for this specific group
  targName = targNames{group_i};
  targetEpochs = targetEpochsParadigm(group_i, :);
  groupLabelsByEvent = plotMat(:, group_i);
  
  % ID Target vs Non-target events
  targInd = max(groupLabelsByEvent);
  targEventInd = groupLabelsByEvent == targInd;                                     % Find the target index
  nonTargEventInd = ~(groupLabelsByEvent == targInd) & groupLabelsByEvent ~= 0;     % 0 is not part of the analysis group.
  
  % Initialize matrix of entries for results.
  allUnitInd = 1;
  epochCompare = find(targetEpochs);
  epochCompareCount = sum(targetEpochs);
  [targVNonTargMat, baselineMat] = deal(zeros(sum(unitCounts), epochCompareCount));
  
  for chan_i = 1:chanCount
    for unit_i = 1:unitCounts(chan_i)
            
      % Collect baseline for the target of the test.
      targBaseline = cellfun(pullRates, spikeCountsByImageByEpoch{1}{chan_i}{unit_i}(targEventInd), 'UniformOutput', false);
      targBaseline = vertcat(targBaseline{:});
      
      for ep_i = 1:epochCompareCount
        % Identify the epoch data to collect
        epInd = epochCompare(ep_i);
        
        % Rates for tests
        targRates = cellfun(pullRates, spikeCountsByImageByEpoch{epInd}{chan_i}{unit_i}(targEventInd), 'UniformOutput', false);
        nonTargRates = cellfun(pullRates, spikeCountsByImageByEpoch{epInd}{chan_i}{unit_i}(nonTargEventInd), 'UniformOutput', false);
        targRates = vertcat(targRates{:});
        nonTargRates = vertcat(nonTargRates{:});
        
        % Compare target epoch vs Baseline.
        if epochStatsParams.nonParametric
          [pVal1, sigSwitch1, ~] = signrank(targRates, targBaseline);
        else
          [sigSwitch1, pVal1, ~] = ttest(targRates, targBaseline);
        end
        
        if isnan(sigSwitch1)
          baselineMat(allUnitInd, ep_i) = nan;
        elseif ~isnan(pVal1) && pVal1 < alpha
          baselineMat(allUnitInd, ep_i) = mean(targRates) - mean(targBaseline);
        end
        
        % If there are non-target events, compare target vs non-target
        % epochs. Otherwise save a NaN.
        if any(nonTargEventInd)
          
          if epochStatsParams.nonParametric
            [pVal, sigSwitch, ~] = ranksum(targRates, nonTargRates);
          else
            [sigSwitch, pVal, ~] = ttest2(targRates, nonTargRates);
          end
          
          if isnan(sigSwitch)
            targVNonTargMat(allUnitInd, ep_i) = nan;
          elseif ~isnan(pVal) && pVal < alpha
            targVNonTargMat(allUnitInd, ep_i) = mean(targRates) - mean(nonTargRates);
          end
          
        else
          targVNonTargMat(allUnitInd, ep_i) = nan;
        end
                
      end
      
      % Increment index
      allUnitInd = allUnitInd + 1;
      
    end
  end
  
  % Save things to the selTable
  for ep_i = 1:epochCompareCount
    
    % baseline difference for target stimulus
    fieldName =  sprintf('%s_baseDiff_%s', targName, labels{epochCompare(ep_i)});
    selTable.(fieldName) = baselineMat(:,ep_i);

    % soc v non-soc differences
    fieldName =  sprintf('%sSel_%s', targName, labels{epochCompare(ep_i)});
    selTable.(fieldName) = targVNonTargMat(:,ep_i);
    
  end
  
end

end

function [selTable] = epochStats(spikesByEvent, selTable, eventIDs, stimulusLabelGroups, epochStatsParams)
% This function performs a few statistical tests to determine whether there
% are meaningful difference between different periods of a task, and
% different conditions of the task.

% Use the simpler cell array of names to see which groups have been removed
% during preprocessing due to being absent.
analysisGroups = stimulusLabelGroups.groups;
[~, group2Analyze] = intersect(epochStatsParams.names, stimulusLabelGroups.names);
group2Analyze = sort(group2Analyze);
epochStatsParams.names = epochStatsParams.names(group2Analyze);
epochStatsParams.targ = epochStatsParams.targ(group2Analyze);
epochStatsParams.targNames = epochStatsParams.targNames(group2Analyze);

timeBins = epochStatsParams.times;
labels = epochStatsParams.labels;
chanCount = length(spikesByEvent{1});
unitCounts = cellfun('length', spikesByEvent{1});

% Bin spikes
spikeCountsByImageByEpoch = cell(size(timeBins,1),1);
for epoch_i = 1:size(timeBins,1)
  [spikeCounts, ~, ~] = spikeCounter(spikesByEvent, timeBins(epoch_i, 1), timeBins(epoch_i, 2));
  spikeCountsByImageByEpoch{epoch_i} = spikeCounts;
  %spikeCountsByImageByEpoch{epoch}{channel}{unit}{event}.rates = trials*1
end

% Determine fixation selectivity
allUnitInd = 1;
fixSel = zeros(sum(unitCounts), 1);

for chan_i = 1:chanCount
  for unit_i = 1:unitCounts(chan_i)
    
    % Check for fixation period selectivity -
    allBaseline = cellfun(@(x) x.rates, spikeCountsByImageByEpoch{1}{chan_i}{unit_i}, 'UniformOutput', false);
    allBaseline = vertcat(allBaseline{:});
    allFix = cellfun(@(x) x.rates, spikeCountsByImageByEpoch{2}{chan_i}{unit_i}, 'UniformOutput', false);
    allFix = vertcat(allFix{:});
    
    % Run the Test
    if epochStatsParams.nonParametric
      [~, sigSwitch0, ~] = signrank(allFix, allBaseline);
    else
      [sigSwitch0, ~, ~] = ttest(allFix, allBaseline);
    end
    
    % If sig, save the difference
    if ~isnan(sigSwitch0) && sigSwitch0
      fixSel(allUnitInd) = mean(allFix) - mean(allBaseline);
    end
  end
end

% Add Results to table
selTable.fixationSel = fixSel;

for group_i = 1:length(epochStatsParams.names)
  
  % Extract relevant variables for this specific group
  target = epochStatsParams.targ{group_i};                        % When performing tests, the label from groups which is used.
  targName = epochStatsParams.targNames{group_i};
  groups = analysisGroups{group_i};
  targetEpochs = epochStatsParams.targetEpochs(group_i, :);
  groupLabelsByEvent = epochStatsParams.groupLabelsByImage(:, group_i);
  
  % ID Target vs Non-target events
  targInd = find(strcmp(groups, target));
  targEventInd = groupLabelsByEvent == targInd;
  if any(contains(eventIDs, 'monkey')) % hack, since SocVNonSoc comparisons should be only for agent containing stim.
    agentInd = ~(contains(eventIDs, 'landscape') | contains(eventIDs, 'objects'));
    nonTargEventInd = ~(groupLabelsByEvent == targInd) & agentInd;
  else
    nonTargEventInd = ~(groupLabelsByEvent == targInd) & groupLabelsByEvent ~= 0;    % 0 is not part of the analysis group.
  end
  
  % Initialize matrix of entries for results.
  allUnitInd = 1;
  [targVNonTargMat, baselineMat] = deal(zeros(sum(unitCounts), sum(targetEpochs)));
  epochCompare = find(targetEpochs);
  epochCompareCount = length(epochCompare);
  
  for chan_i = 1:chanCount
    for unit_i = 1:unitCounts(chan_i)
            
      % Collect baseline for the target of the test.
      targBaseline = cellfun(@(x) x.rates, spikeCountsByImageByEpoch{1}{chan_i}{unit_i}(targEventInd), 'UniformOutput', false);
      targBaseline = vertcat(targBaseline{:});
      
      for ep_i = 1:epochCompareCount
        % Identify the epoch data to collect
        epInd = epochCompare(ep_i);
        
        % Rates for tests
        targRates = cellfun(@(x) x.rates, spikeCountsByImageByEpoch{epInd}{chan_i}{unit_i}(targEventInd), 'UniformOutput', false);
        nonTargRates = cellfun(@(x) x.rates, spikeCountsByImageByEpoch{epInd}{chan_i}{unit_i}(nonTargEventInd), 'UniformOutput', false);
        targRates = vertcat(targRates{:});
        nonTargRates = vertcat(nonTargRates{:});
        
        % Compare target epoch vs Baseline.
        if epochStatsParams.nonParametric
          [~, sigSwitch2, ~] = signrank(targRates, targBaseline);
        else
          [sigSwitch2, ~, ~] = ttest(targRates, targBaseline);
        end
        
        if isnan(sigSwitch2)
          baselineMat(allUnitInd, ep_i) = nan;
        elseif ~isnan(sigSwitch2) && sigSwitch2
          baselineMat(allUnitInd, ep_i) = mean(targRates) - mean(targBaseline);
        end
        
        % If there are non-target events, compare target vs non-target
        % epochs. Otherwise save a NaN.
        if any(nonTargEventInd)
          
          if epochStatsParams.nonParametric
            [~, sigSwitch, ~] = ranksum(targRates, nonTargRates);
          else
            [sigSwitch, ~, ~] = ttest2(targRates, nonTargRates);
          end
          
          if isnan(sigSwitch)
            targVNonTargMat(allUnitInd, ep_i) = nan;
          elseif ~isnan(sigSwitch) && sigSwitch
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
  
  % Add to the output table

  % baseline difference
  for ep_i = 1:epochCompareCount
    fieldName =  sprintf('%s_baseDiff_%s', targName, labels{epochCompare(ep_i)});
    selTable.(fieldName) = baselineMat(:,ep_i);
  end
  
  % soc v non-soc differences
  for ep_i = 1:epochCompareCount
    fieldName =  sprintf('%sSel_%s', targName, labels{epochCompare(ep_i)});
    selTable.(fieldName) = targVNonTargMat(:,ep_i);
  end
  
end

end

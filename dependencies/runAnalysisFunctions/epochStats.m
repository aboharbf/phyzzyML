function [selTable] = epochStats(spikesByEvent, selTable, eventIDs, stimulusLabelGroups, epochStatsParams)
% This function performs a few statistical tests to determine whether there
% are meaningful difference between different periods of a task, and
% different conditions of the task.

% Non-parametric alternatives
% Use 'signrank' - https://www.mathworks.com/help/stats/signrank.html
% between rates for baseline and rates in period of interest.
% Use 'ranksum' - https://www.mathworks.com/help/stats/ranksum.html between
% rates for social and non-social conditions.

epochStatsParams.stimParamsFilename = epochStatsParams.stimParamsFilename;
epochStatsParams.plotLabels = stimulusLabelGroups.groups;
epochStatsParams.outLogic = 0;
epochStatsParams.removeEmpty = 0;

eventCatMat = plotIndex(eventIDs, epochStatsParams);

[~, group2Analyze] = intersect(epochStatsParams.names, stimulusLabelGroups.names);
epochStatsParams.names = epochStatsParams.names(group2Analyze);

for group_i = 1:length(stimulusLabelGroups.groups)
  
  % Extract relevant variables for this specific group
  target = stimulusLabelGroups.groups{group_i}{1};                        % When performing a t test, the label from groups which is used. the rest are 'non-' label.
  groups = stimulusLabelGroups.groups{group_i};
  timeBins = epochStatsParams.times;
  labels = epochStatsParams.labels;
  targetEpochs = epochStatsParams.targetEpochs(group_i, :);
  groupLabelsByEvent = epochStatsParams.groupLabelsByImage(:, group_i);
  
  % Bin spikes
  spikeCountsByImageByEpoch = cell(size(timeBins,1),1);
  for epoch_i = 1:size(timeBins,1)
    [spikeCounts, ~, ~] = spikeCounter(spikesByEvent, timeBins(epoch_i, 1), timeBins(epoch_i, 2));
    spikeCountsByImageByEpoch{epoch_i} = spikeCounts;
  end
  
  %spikeCountsByImageByEpoch{epoch}{channel}{unit}{event}.rates = trials*1
  
  % ID Target vs Non-target events
  targInd = find(strcmp( groups, target));
  targEventInd = groupLabelsByEvent == targInd;
  if any(contains(eventIDs, 'monkey')) % Hack, since non-headTurn isn't labeled here.
    agentInd = ~(contains(eventIDs, 'landscape') | contains(eventIDs, 'objects'));
    nonTargEventInd = ~(groupLabelsByEvent == targInd) & agentInd;
  else
    nonTargEventInd = ~(groupLabelsByEvent == targInd) & groupLabelsByEvent ~= 0;    % 0 is not part of the analysis group.
  end
  
  chanCount = length(spikeCountsByImageByEpoch{1});
  unitCounts = cellfun('length', spikeCountsByImageByEpoch{1});
  allUnitInd = 1;
  [targVNonTargMat, baselineMat] = deal(zeros(sum(unitCounts), sum(targetEpochs)));
  fixSel = targVNonTargMat(:,1);
  epochCompare = find(targetEpochs);
  epochCompareCount = length(epochCompare);
  
  for chan_i = 1:chanCount
    for unit_i = 1:unitCounts(chan_i)
      
      % Check for fixation period selectivity -
      allBaseline = cellfun(@(x) x.rates, spikeCountsByImageByEpoch{1}{chan_i}{unit_i}, 'UniformOutput', false);
      allBaseline = vertcat(allBaseline{:});
      allFix = cellfun(@(x) x.rates, spikeCountsByImageByEpoch{2}{chan_i}{unit_i}, 'UniformOutput', false);
      allFix = vertcat(allFix{:});
      [sigSwitch0, ~, ci] = ttest(allFix, allBaseline);
      if ~isnan(sigSwitch0) && sigSwitch0
        fixSel(allUnitInd) = mean(ci);
      end
      
      % Collect baseline for the target of the test.
      targBaseline = cellfun(@(x) x.rates, spikeCountsByImageByEpoch{1}{chan_i}{unit_i}(targEventInd), 'UniformOutput', false);
      targBaseline = vertcat(targBaseline{:});
      
      for ep_i = 1:epochCompareCount
        % identify the epoch data to collect
        epInd = epochCompare(ep_i);
        
        % Rates for tests
        targRates = cellfun(@(x) x.rates, spikeCountsByImageByEpoch{epInd}{chan_i}{unit_i}(targEventInd), 'UniformOutput', false);
        nonTargRates = cellfun(@(x) x.rates, spikeCountsByImageByEpoch{epInd}{chan_i}{unit_i}(nonTargEventInd), 'UniformOutput', false);
        targRates = vertcat(targRates{:});
        nonTargRates = vertcat(nonTargRates{:});
        
        % signrank to compare epoch to baseline.
        [sigSwitch2, ~, ci] = ttest(targRates, targBaseline);
        if isnan(sigSwitch2)
          baselineMat(allUnitInd, ep_i) = nan;
        elseif ~isnan(sigSwitch2) && sigSwitch2
          baselineMat(allUnitInd, ep_i) = mean(ci);
        end
        
        % Running tests - ttest2 to compare soc v non soc
        if any(nonTargEventInd)
          [sigSwitch, ~, ci] = ttest2(targRates, nonTargRates);
          if isnan(sigSwitch)
            targVNonTargMat(allUnitInd, ep_i) = nan;
          elseif ~isnan(sigSwitch) && sigSwitch
            targVNonTargMat(allUnitInd, ep_i) = mean(ci);
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
  
  % Fixation selectivity
  selTable.fixationSel = fixSel;
  
  % baseline difference
  for ep_i = 1:epochCompareCount
    fieldName =  sprintf('baselineDiff_%s', labels{epochCompare(ep_i)});
    selTable.(fieldName) = targVNonTargMat(:,ep_i);
  end
  
  % soc v non-soc differences
  for ep_i = 1:epochCompareCount
    fieldName =  sprintf('%sSel_%s', target, labels{epochCompare(ep_i)});
    selTable.(fieldName) = baselineMat(:,ep_i);
  end
  
end

end

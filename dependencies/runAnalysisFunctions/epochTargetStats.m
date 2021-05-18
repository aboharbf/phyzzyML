function [selTable] = epochTargetStats(spikesByEvent, selTable, eventIDs, paradigm, epochStatsParams)
% This function performs a few statistical tests to determine whether there
% are meaningful difference between different conditions of the task during
% different predefined time periods.

% Bin spikes
timeBins = epochStatsParams.times;
timeLabels = epochStatsParams.labels;
spikeCountsByImageByEpoch = cell(size(timeBins,1),1);
for epoch_i = 1:size(timeBins,1)
  [spikeCounts, ~, ~] = spikeCounter(spikesByEvent, timeBins(epoch_i, 1), timeBins(epoch_i, 2));
  spikeCountsByImageByEpoch{epoch_i} = spikeCounts;
end

%% Determine Selectivity per Epoch selectivity, comparing target and non-target

pStruct = epochStatsParams.(paradigm);
targNames = pStruct.targNames;
targLabelList = pStruct.targ;
targetEpochsParadigm = pStruct.targetEpochs;
oneVsAllSwitch = pStruct.oneVsAll;

labels = epochStatsParams.labels;
alpha = 0.05;

pullRates = @(x) x.rates;
chanCount = length(spikesByEvent{1});
unitCounts = cellfun('length', spikesByEvent{1});

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
  
  [targVBase_pValMat, targVBase_diffMat, targVnonTarg_pValMat, targVnonTarg_diffMat] = deal(zeros(sum(unitCounts), epochCompareCount));
  [selTablePref, ~] = deal(cell(sum(unitCounts), epochCompareCount));
  
  for chan_i = 1:chanCount
    for unit_i = 1:unitCounts(chan_i)
      
      % Collect baseline for the target of the test.
      targBaseline = cellfun(pullRates, spikeCountsByImageByEpoch{1}{chan_i}{unit_i}(targEventInd), 'UniformOutput', false);
      targBaseline = vertcat(targBaseline{:});
      
      for ep_i = 1:epochCompareCount
        % Identify the epoch data to collect
        epInd = epochCompare(ep_i);
        
        if oneVsAllSwitch(group_i)
          % Perform a test to see if the epoch-target combo differs from
          % baseline, and if its activity against non-target epochs is
          % different.
          
          % Rates for tests
          targRates = cellfun(pullRates, spikeCountsByImageByEpoch{epInd}{chan_i}{unit_i}(targEventInd), 'UniformOutput', false);
          targRates = vertcat(targRates{:});
          
          % Compare target epoch vs Baseline.
          if epochStatsParams.nonParametric
            [pVal1, ~, ~] = signrank(targRates, targBaseline);
          else
            [~, pVal1, ~] = ttest(targRates, targBaseline);
          end
          
          % store difference b/t means and pVal
          targVBase_pValMat(allUnitInd, ep_i) = pVal1;
          targVBase_diffMat(allUnitInd, ep_i) = mean(targRates) - mean(targBaseline);
          
          % If there are non-target events, compare target vs non-target
          % epochs. Otherwise save a NaN.
          if any(nonTargEventInd)
            % Collect the rates for the non-Target
            nonTargRates = cellfun(pullRates, spikeCountsByImageByEpoch{epInd}{chan_i}{unit_i}(nonTargEventInd), 'UniformOutput', false);
            nonTargRates = vertcat(nonTargRates{:});
            
            % Run the statistical test to compare target vs non-target.
            if epochStatsParams.nonParametric
              [pVal, ~, ~] = ranksum(targRates, nonTargRates);
            else
              [~, pVal, ~] = ttest2(targRates, nonTargRates);
            end
            
            % Store the mean difference
            targVnonTarg_pValMat(allUnitInd, ep_i) = pVal;
            targVnonTarg_diffMat(allUnitInd, ep_i) = mean(targRates) - mean(nonTargRates);
            
          else
            
            % Place NaN if there aren't any non-targets to compare with.
            targVnonTarg_pValMat(allUnitInd, ep_i) = NaN;
            targVnonTarg_diffMat(allUnitInd, ep_i) = 0;
            
          end
          
        else
          % In this case, you have more than a binary choice, and you want
          % to know if epochs belonging to different members of a group are
          % different from each other. Perform the test with ANOVA.
          
          % Grab labels, initialize things.
          groupLabels = targLabelList{group_i};
          realMeans = zeros(length(groupLabels),1);
          [groupRates, groupLabelIn] = deal([]);
          
          % Cycle through each label, collecting rates and creating a label
          % vector alongside it.
          for label_i = 1:length(groupLabels)
            tmp = cellfun(pullRates, spikeCountsByImageByEpoch{epInd}{chan_i}{unit_i}(groupLabelsByEvent == label_i), 'UniformOutput', false);
            data2Add = vertcat(tmp{:});
            realMeans(label_i) = mean(data2Add);
            groupRates = [groupRates; data2Add];
            groupLabelIn = [groupLabelIn; repmat(string(groupLabels{label_i}), [length(data2Add),1])];
          end
          
          % run the ANOVA. 
          [targVnonTarg_pValMat(allUnitInd, ep_i), ~, C] = anovan(groupRates, {groupLabelIn}, 'display', 'off');
         
          % Do multiple comparisons
          [multTable, est] = multcompare(C, 'display', 'off');
          multTable = multTable(multTable(:,end) <= alpha, :);
          if ~isempty(multTable)
            groupsWithDiffs = multTable(:,1:2);
            groupsWithDiffs = unique(groupsWithDiffs(:));
            % See what the highest rate is across groups
            estSig = est(groupsWithDiffs,:);
            
            % Find the highest rate participating in a significant
            % comparison.
            
            [~, frInd] = sort(estSig(:,1), 'descend');
            selTablePref{allUnitInd, ep_i} = groupLabels{groupsWithDiffs(frInd(1))};
          else
            selTablePref{allUnitInd, ep_i} = 'None';
          end
                    
        end
        
      end
      
      % Increment index
      allUnitInd = allUnitInd + 1;
      
    end
  end
  
  % Save things to the selTable for the group
  for ep_i = 1:epochCompareCount
    
    % Some pVal mat is always produced - save it
    fieldName =  sprintf('epochSel_%s_%s_pVal', targName, labels{epochCompare(ep_i)});
    selTable.(fieldName) = targVnonTarg_pValMat(:,ep_i);

    if oneVsAllSwitch(group_i)
      % in one vs all, t tests are performed, store pVals and differences
      fieldName =  sprintf('epochSel_%s_%s_diff', targName, labels{epochCompare(ep_i)});
      selTable.(fieldName) = targVnonTarg_diffMat(:,ep_i);
      
      % Comparisons Against baseline
      fieldName =  sprintf('epochSel_%s_%s_vBase_pVal', targName, labels{epochCompare(ep_i)});
      selTable.(fieldName) = targVBase_pValMat(:,ep_i);
      
      fieldName =  sprintf('epochSel_%s_%s_vBase_diff', targName, labels{epochCompare(ep_i)});
      selTable.(fieldName) = targVBase_diffMat(:,ep_i);
      
    else
      % in non-one vs all tests, ANOVA is performed, store the preferred
      % group in the selTable.
      
      fieldName =  sprintf('epochSel_%s_%s_prefStim', targName, labels{epochCompare(ep_i)});
      selTable.(fieldName) = selTablePref(:,ep_i);
      
    end
    
  end
  
end

end

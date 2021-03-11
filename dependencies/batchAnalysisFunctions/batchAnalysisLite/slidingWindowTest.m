function spikeDataBank = slidingWindowTest(spikeDataBank,params, figStruct)
disp('Starting sliding window test...');
% Perform sliding scale ANOVA, calculate Omega at each bin.

if ~exist(params.outputDir, 'dir')
  mkdir(params.outputDir)
end

% Generate variables needed
runList = fields(spikeDataBank);
binSize = params.binSize;
binStep = params.binStep;
target = params.target;
Omega = params.Omega;
stimParamFile = params.stimParamFile;
nestStrCmp = @(x, y) any(strcmp(x, y));

% Step 1 - for each bin (epoch), calculate the rates and counts. save them
% into a larger epochRates structure, which can be stored in the folder for
% this function for subsequent runs. 

if ~exist(fullfile(params.outputDir, 'epochRates.mat'), 'file')
  epochRates = cell(length(runList),1);
  catagoryInd = cell(length(runList),1);
  % Step 1 - generate bin times and spike rates, and proper memberships to groups.
  for run_ind = 1:length(runList)
    runStruct = spikeDataBank.(runList{run_ind});
    starts = (runStruct.start:binStep:(runStruct.end - binSize))';
    ends = (runStruct.start+binSize:binStep:(runStruct.end))';
    spikeDataBank.(runList{run_ind}).epochs = [starts,ends];
    epochRates{run_ind} = cell(length(starts),1);
    for bin_ind = 1:length(starts)
      [epochRates{run_ind}{bin_ind}, ~, ~] = spikeCounter(spikeDataBank.(runList{run_ind}).spikesByEvent, starts(bin_ind), ends(bin_ind));
    end
    for group_ind = 1:length(target)
      catagoryInd{run_ind}(:,group_ind) = cell2mat(cellfun(nestStrCmp, runStruct.eventCategories, repmat(target(group_ind), [length(runStruct.eventCategories),1]), 'UniformOutput',0));
    end
  end
  save(fullfile(params.outputDir, 'epochRates'), 'epochRates', 'catagoryInd');
else
  load(fullfile(params.outputDir, 'epochRates'), 'epochRates', 'catagoryInd');
  assert(length(epochRates) == length(runList), 'epochRates length does not match runList length');
  fprintf('Sliding windows rates already calculated, continuing... \n');
end

% Step 2 - Perform ANOVA, Omega calculation across bins, and store values.
% Save the output in the function folder for later retrieval. 

statsType = {'groupRates', 'nonGroupRates', 'pVec', 'cohensVec', 'nullPVec', 'nullCohensVec', 'sigBins'};

if ~exist(fullfile(params.outputDir, 'allStatsArray.mat'), 'file')
  chanUnitStruct = structfun(@(x) x.psthByImage, spikeDataBank,'UniformOutput',0);
  chanUnitStruct = struct2cell(chanUnitStruct);
  allStatsArray = initNestedCellArray(chanUnitStruct, 'cell', [length(epochRates{1}), length(target), length(statsType)], 3);
  clear chanUnitStruct;
  scrambleCount = params.scrambleCount;
  
  tic
  if license('test','Distrib_Computing_Toolbox')
    
    epochRatesPar = parallel.pool.Constant(epochRates);
    catagoryIndPar =  parallel.pool.Constant(catagoryInd);
%     epochRatesPar.Value = epochRates; for debugging.
%     catagoryIndPar.Value =  catagoryInd;
    parfor run_ind = 1:length(runList)
      % Initialize relevant Structures
      if length(unique(catagoryIndPar.Value{run_ind})) > 1 %Only run tests where you have members in each group.
        binCount = length(epochRatesPar.Value{run_ind});
        for chan_ind = 1:length(epochRatesPar.Value{run_ind}{1})
          for unit_ind = 1:length(epochRatesPar.Value{run_ind}{1}{chan_ind})
            for bin_ind = 1:binCount
              unitResponsePerEvent = epochRatesPar.Value{run_ind}{bin_ind}{chan_ind}{unit_ind};
              catagortyIndSlice = catagoryIndPar.Value{run_ind};
              %unitData{event}.rates = trial*1
              for target_ind = 1:length(target)
                if length(unique(catagortyIndSlice(:,target_ind))) == 1
                  break % If there is only 1 label type, there is no comparison to be made.
                end
                [trialSpikes, trialLabels]  = deal([]);
                % grab the relevant events
                targetInd = catagortyIndSlice(:,target_ind);
                targetSpikes = unitResponsePerEvent(targetInd);
                otherSpikes = unitResponsePerEvent(~targetInd);
                % Initialize relevant vecotrs
                spikeGroups = {targetSpikes otherSpikes};
                spikeGroupLabels ={(target{target_ind}) (['non-' target{target_ind}])};
                % Cluster and reshape the arrays properly
                for group_i = 1:length(spikeGroups)
                  tmp = spikeGroups{group_i};
                  tmp = [tmp{:}];
                  dataVec = vertcat(tmp.rates);
                  labelVec = repmat(spikeGroupLabels(group_i), length(dataVec),1);
                  trialSpikes = vertcat(trialSpikes,dataVec);
                  trialLabels = vertcat(trialLabels, labelVec);
                end
                % Check for social v non-social
                % 3rd ind in 3D mat statsType = {'pVec', 'cohensVec', 'nullPVec', 'nullCohensVec'};
                pop1 = trialSpikes(strcmp(trialLabels,spikeGroupLabels{1}));
                pop2 = trialSpikes(strcmp(trialLabels,spikeGroupLabels{2}));
                allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('groupRates',statsType)} = pop1;
                allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('nonGroupRates',statsType)} = pop2;
                %[allStatsArray{run_ind}{chan_ind}{unit_ind}(bin_ind,target_ind, strcmp('pVec',statsType)), pStatsTable, ~] = anovan(trialSpikes,{trialLabels},'model','interaction','varnames',{'SvNS'}, 'alpha', 0.05,'display','off');
                [allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('sigBins',statsType)}, allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('pVec',statsType)}, ~, pStatsTable] = ttest2(pop1, pop2); %,'model','interaction','varnames',{'SvNS'}, 'alpha', 0.05,'display','off');
                allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('cohensVec',statsType)} = (mean(pop1) - mean(pop2))/pStatsTable.sd;
                [nullPVec, nullCohensVec] = deal(zeros(1,scrambleCount));
                for rand_ind = 1:scrambleCount
                  trialLabelsTmp = trialLabels(randperm(length(trialLabels)));
                  pop1 = trialSpikes(strcmp(trialLabelsTmp,spikeGroupLabels{1}));
                  pop2 = trialSpikes(strcmp(trialLabelsTmp,spikeGroupLabels{2}));
                  [~, nullPVec(rand_ind), ~, pStatsTableTmp] = ttest2(pop1, pop2);
                  nullCohensVec(rand_ind) = (mean(pop1) - mean(pop2))/pStatsTableTmp.sd;
                end
                allStatsArray{run_ind}{chan_ind}{unit_ind}(bin_ind,target_ind, strcmp('nullPVec',statsType)) = {nullPVec};
                allStatsArray{run_ind}{chan_ind}{unit_ind}(bin_ind,target_ind, strcmp('nullCohensVec',statsType)) = {nullCohensVec};
              end
            end
          end
        end
      else
        disp('Skipping')
      end
    end
    
  else
    
    for run_ind = 1:length(runList)
      % Initialize relevant Structures
      if length(unique(catagoryInd{run_ind})) > 1 %Only run tests where you have members in each group.
        chanCount = length(epochRates{run_ind}{1});
        binCount = length(epochRates{run_ind});
        for chan_ind = 1:chanCount
          for unit_ind = 1:length(epochRates{run_ind}{1}{chan_ind})
            for bin_ind = 1:binCount
              unitResponsePerEvent = epochRates{run_ind}{bin_ind}{chan_ind}{unit_ind};
              catagortyIndSlice = catagoryInd{run_ind};
              %unitData{event}.rates = trial*1
              for target_ind = 1:length(target)
                if length(unique(catagortyIndSlice(:,target_ind))) == 1
                  break % If there is only 1 label type, there is no comparison to be made.
                end
                [trialSpikes, trialLabels]  = deal([]);
                % grab the relevant events
                targetInd = catagortyIndSlice(:,target_ind);
                targetSpikes = unitResponsePerEvent(targetInd);
                otherSpikes = unitResponsePerEvent(~targetInd);
                % Initialize relevant vecotrs
                spikeGroups = {targetSpikes otherSpikes};
                spikeGroupLabels ={(target{target_ind}) (['non-' target{target_ind}])};
                % Cluster and reshape the arrays properly
                for group_i = 1:length(spikeGroups)
                  tmp = spikeGroups{group_i};
                  tmp = [tmp{:}];
                  dataVec = vertcat(tmp.rates);
                  labelVec = repmat(spikeGroupLabels(group_i), length(dataVec),1);
                  trialSpikes = vertcat(trialSpikes,dataVec);
                  trialLabels = vertcat(trialLabels, labelVec);
                end
                % Check for social v non-social
                % 3rd ind in 3D mat statsType = {'pVec', 'cohensVec', 'nullPVec', 'nullCohensVec'};
                pop1 = trialSpikes(strcmp(trialLabels,spikeGroupLabels{1}));
                pop2 = trialSpikes(strcmp(trialLabels,spikeGroupLabels{2}));
                allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('groupRates',statsType)} = pop1;
                allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('nonGroupRates',statsType)} = pop2;
                %[allStatsArray{run_ind}{chan_ind}{unit_ind}(bin_ind,target_ind, strcmp('pVec',statsType)), pStatsTable, ~] = anovan(trialSpikes,{trialLabels},'model','interaction','varnames',{'SvNS'}, 'alpha', 0.05,'display','off');
                [allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('sigBins',statsType)}, allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('pVec',statsType)}, ~, pStatsTable] = ttest2(pop1, pop2); %,'model','interaction','varnames',{'SvNS'}, 'alpha', 0.05,'display','off');
                allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('cohensVec',statsType)} = (mean(pop1) - mean(pop2))/pStatsTable.sd;
                [nullPVec, nullCohensVec] = deal(zeros(1,scrambleCount));
                for rand_ind = 1:scrambleCount
                  trialLabelsTmp = trialLabels(randperm(length(trialLabels)));
                  pop1 = trialSpikes(strcmp(trialLabelsTmp,spikeGroupLabels{1}));
                  pop2 = trialSpikes(strcmp(trialLabelsTmp,spikeGroupLabels{2}));
                  [~, nullPVec(rand_ind), ~, pStatsTableTmp] = ttest2(pop1, pop2);
                  nullCohensVec(rand_ind) = (mean(pop1) - mean(pop2))/pStatsTableTmp.sd;
                end
                allStatsArray{run_ind}{chan_ind}{unit_ind}(bin_ind,target_ind, strcmp('nullPVec',statsType)) = {nullPVec};
                allStatsArray{run_ind}{chan_ind}{unit_ind}(bin_ind,target_ind, strcmp('nullCohensVec',statsType)) = {nullCohensVec};
              end
            end
          end
        end
      else
        disp('Skipping')
      end
    end
    
  end
  save(fullfile(params.outputDir, 'allStatsArray'), 'allStatsArray');
  fprintf('Done in %d hours \n', toc/3600)
else
  fprintf('Loading statsArray from saved data... \n')
  load(fullfile(params.outputDir, 'allStatsArray.mat'), 'allStatsArray');
  assert(length(allStatsArray) == length(runList), 'allStatsArray length does not match runList length');
end

runStruct = spikeDataBank.(runList{1});
starts = (runStruct.start:binStep:(runStruct.end - binSize))';
xlabelPlot = starts(1:8:length(starts));
stimOnInd = find(starts == 0);
stimEndInd = find(starts == 2800);

if params.plotTest
  for run_i = 1:length(allStatsArray)
    statsArray = allStatsArray{run_i};
    % Plot the results for each unit seen
    for chan_ind = 1:length(statsArray)
      for unit_ind = 1:length(statsArray{chan_ind})
        if unit_ind == 1
          ANOVAvarName = ['Ch' num2str(chan_ind) ' Unsorted - SocVsNonSoc'];
        elseif unit_ind == length(statsArray{chan_ind})
          ANOVAvarName = ['Ch' num2str(chan_ind) ' MUA - SocVsNonSoc'];
        else
          ANOVAvarName = ['Ch' num2str(chan_ind) ' U' num2str(unit_ind-1) ' - SocVsNonSoc'];
        end
        % {'groupRates', 'nonGroupRates', 'pVec', 'cohensVec', 'nullPVec', 'nullCohensVec'};
        titleString = sprintf("%s - Sliding Window Group mean comparison - %s", runList{run_i} ,ANOVAvarName);
        fileString = fullfile(spikeDataBank.(runList{run_i}).figDir, [titleString{1}]);
        if ~exist([fileString '.fig'], 'file') || (~exist([fileString '.png'], 'file') && params.exportFig)
          fullfile(spikeDataBank.(runList{run_i}).figDir, titleString)
          groupMeans = cellfun(@(x) mean(x), statsArray{chan_ind}{unit_ind}(:,:,strcmp('groupRates',statsType)));
          if ~(nnz(groupMeans) == 0)
            groupSEMs = cellfun(@(x) (std(x)/sqrt(length(x))), statsArray{chan_ind}{unit_ind}(:,:,strcmp('groupRates',statsType)));
            nonGroupMeans = cellfun(@(x) mean(x), statsArray{chan_ind}{unit_ind}(:,:,strcmp('nonGroupRates',statsType)));
            nonGroupSEMs = cellfun(@(x) (std(x)/sqrt(length(x))), statsArray{chan_ind}{unit_ind}(:,:,strcmp('nonGroupRates',statsType)));
            pVec = cell2mat(statsArray{chan_ind}{unit_ind}(:,:,strcmp('pVec',statsType)));
            %cohensVec = cell2mat(statsArray{chan_ind}{unit_ind}(:,:,strcmp('cohensVec',statsType)));
            if size(pVec,2) ~= length(target)
              tmp = statsArray{chan_ind}{unit_ind}(:,:,strcmp('pVec',statsType));
              emptyInd = cellfun(@(x) isempty(x), tmp);
              tmp(emptyInd) = deal({1});
              pVec = cell2mat(tmp);
            end
            % plot
            h = figure('NumberTitle', 'off', 'Name', titleString,'units','normalized','outerposition',[0 0 params.plotSize]);
            sgtitle(titleString)
            for target_i = 1:length(target)
              subplot(length(target),1,target_i)
              title(target{target_i});
              lineData = [groupMeans(:, target_i), nonGroupMeans(:, target_i)];
              errorData = [groupSEMs(:, target_i), nonGroupSEMs(:, target_i)];
              mseb(1:length(lineData), lineData', errorData');
              legend([target(target_i), ['non-' target{target_i}]], 'location', 'northwest','AutoUpdate','off')
              hold on
              xlabel('Bin Start Time from Stimulus Onset (ms)')
              ylabel('Firing Rates (Hz)')
              xlim([0,length(lineData)]);
              xticks(1:8:length(starts))
              xticklabels(xlabelPlot);
              plot([stimOnInd, stimOnInd], ylim(),'k','linewidth',4);
              plot([stimEndInd, stimEndInd], ylim(),'k','linewidth',4);
              
              % Star Significant Bins
              sigBins = find(pVec(:, target_i) < 0.05);
              lineMax = max(lineData,[],2);
              for star_i = 1:length(sigBins)
                text(sigBins(star_i),(lineMax(sigBins(star_i))*1.05),'*','FontSize',20)
              end
            end
            % Save and close
            % (outDir, filename, figData, saveFig, exportFig, saveData, varargin )
            saveFigure(spikeDataBank.(runList{run_i}).figDir, titleString, [], figStruct, []);
            close(h)
          end
        end
      end
    end
  end
end

fprintf('Sliding Window Calculations finished... \n')

% Step 3 - Count Results of ANOVA across all units and count stretches of significant
% bins.
totalUnitCount = 0;
totalChannelCount = 0;
groupingType = {'Unsorted', 'Unit', 'MUA'};
% target = {'socialInteraction','agents','interaction'};, pulled from params.
dataType = {'binTotal','sigBins','nullSigBins'};
statType = {'Counts', 'Sig Run Lengths', 'Mean Stretch Lengths', 'Max run lengths', 'Sig Run Starts'};
cohensStatType = {'Bin Value', 'Maximum values'};

% Find Bins in each Epoch of the stimulus.
starts = (runStruct.start:binStep:(runStruct.end - binSize))';
epochNames = {'Early','Middle','Late'};
epochBinTimes = [[find(starts == 50) find(starts == 1025)];...
                [find(starts == 1050) find(starts == 2025)];...
                [find(starts == 2050) find(starts == 2850)]];
epochBinSig = zeros(length(epochNames), length(groupingType) , length(target));

testMat = zeros(length(epochRates{1}), length(target), length(groupingType), length(dataType)); % (binCount, target, group, dataType)
targetRunLens = cell(length(target), length(groupingType), length(statType));
unitSigEpochData = {'Count', 'Run Ch Unit Index'};
unitSigEpoch = cell(length(epochNames), length(groupingType) , length(target), length(unitSigEpochData));
cohensData = cell(length(target), length(groupingType), length(cohensStatType));

for run_ind = 1:length(allStatsArray)
  for chan_ind = 1:length(allStatsArray{run_ind})
    unitCount = length(allStatsArray{run_ind}{chan_ind});
    totalUnitCount = totalUnitCount + (unitCount - 2);
    totalChannelCount = totalChannelCount + 1;
    for unit_ind = 1:unitCount
      % {'groupRates', 'nonGroupRates', 'pVec', 'cohensVec', 'nullPVec', 'nullCohensVec'};
      uPVec = cell2mat(allStatsArray{run_ind}{chan_ind}{unit_ind}(:,:,strcmp('pVec',statsType)));
      if size(uPVec,2) ~= length(target)
        tmp = allStatsArray{run_ind}{chan_ind}{unit_ind}(:,:,strcmp('pVec',statsType));
        emptyInd = cellfun(@(x) isempty(x), tmp);
        tmp(emptyInd) = deal({0});
        uPVec = cell2mat(tmp);
      end
      
      if ~isempty(uPVec)
        % Check for significant bins
        uPVec(uPVec == 0) = nan;
        tmpSigCount = double(uPVec < 0.05);
        tmpNullSigCount = cellfun(@(x) sum(x < 0.05), allStatsArray{run_ind}{chan_ind}{unit_ind}(:,:,strcmp('nullPVec',statsType)));
        cohensVec = allStatsArray{run_ind}{chan_ind}{unit_ind}(:,:,strcmp('cohensVec',statsType));
%         omegaVec(isnan(omegaVec)) = 0;
        % Add to relevant structures
        if unit_ind == 1
          group_i = 1;
        elseif unit_ind == unitCount
          group_i = 3;
        else
          group_i = 2;
        end
        % Save significant bins and null bins correctly.
        %{'binTotal','sigBins','nullSigBins'}
        testMat(:, :, group_i, strcmp(dataType, 'binTotal')) = testMat(:, :, group_i, strcmp(dataType, 'binTotal')) + ~isnan(uPVec);
        testMat(:, :, group_i, strcmp(dataType, 'sigBins')) = testMat(:, :, group_i, strcmp(dataType, 'sigBins')) + tmpSigCount;
        testMat(:, :, group_i, strcmp(dataType, 'nullSigBins')) = testMat(:, :, group_i, strcmp(dataType, 'nullSigBins')) + tmpNullSigCount;
        % Keep count of consecutive bins
        for targ_i = 1:length(target)
          trace = tmpSigCount(:,targ_i);
          trace(isnan(trace)) = 0;
          starts = find(diff([0; trace]) == 1);
          ends = find(diff([trace; 0]) == -1)+1;
          runLengths = ends - starts;
          % {'Counts', 'Sig Run Lengths', 'Sig Run Starts', 'Mean Stretch Lengths', 'Max run lengths'};
          targetRunLens{targ_i, group_i, strcmp(statType, 'Counts')} = [targetRunLens{targ_i, group_i, strcmp(statType, 'Counts')}; sum(trace)];
          targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Lengths')} = [targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Lengths')}; runLengths];
          targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Starts')} = [targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Starts')}; starts];
          targetRunLens{targ_i, group_i, strcmp(statType, 'Mean Stretch Lengths')} = [targetRunLens{targ_i, group_i, strcmp(statType, 'Mean Stretch Lengths')}; mean(runLengths)];
          targetRunLens{targ_i, group_i, strcmp(statType, 'Max run lengths')} = [targetRunLens{targ_i, group_i, strcmp(statType, 'Max run lengths')}; max(runLengths)];
          
          cohensTrace = cell2mat(cohensVec(:,targ_i)');
          cohensData{targ_i, group_i, 1} = [cohensData{targ_i, group_i, 1}; cohensTrace];
          cohensData{targ_i, group_i, 2} = [cohensData{targ_i, group_i, 2}; max(cohensTrace)];
          
          %See if the unit has a significant stretch of bins during each
          %epoch
          if any(runLengths >= 3)
            
            runsOfInterest = runLengths >= 3;
            startsOfInterest = starts(runsOfInterest);
            for epoch_i = 1:length(epochNames)
              % Check where run belongs
              % unitSigEpoch = (epochNames), (groupingType) ,(target));
              startBin = epochBinTimes(epoch_i,1);
              endBin = epochBinTimes(epoch_i,2);
              newSigRuns = ((startsOfInterest >= startBin) + (startsOfInterest <= endBin) == 2);
              if any(newSigRuns)
                %{'Count', 'Run Ch Unit Index'};
                unitSigEpoch{epoch_i, group_i, targ_i, strcmp(unitSigEpochData, 'Count')} = [unitSigEpoch{epoch_i, group_i, targ_i, strcmp(unitSigEpochData, 'Count')};  sum(newSigRuns)];
                unitSigEpoch{epoch_i, group_i, targ_i, strcmp(unitSigEpochData, 'Run Ch Unit Index')} = [unitSigEpoch{epoch_i, group_i, targ_i, strcmp(unitSigEpochData, 'Run Ch Unit Index')};  [run_ind chan_ind unit_ind]];
              end
            end

          end
        end
      end
      
    end
  end
end

% for targ_i = 1:length(target)
%   for group_i = 1:length(groupingType)
%     binLengths = targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Lengths')};
%     binStarts = targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Starts')};
%     longBinInd = binLengths >= 3;
%     binStarts(~longBinInd) = 0; 
%     binLengths(~longBinInd) = 0; 
%     for epoch_i = 1:length(epochNames)
%       startBin = epochBinTimes(epoch_i,1);
%       endBin = epochBinTimes(epoch_i,2);
%       binsOfInterest = ((binStarts >= startBin) + (binStarts <= endBin) == 2);
%       epochBinSig(epoch_i, group_i, targ_i) = sum(binLengths(binsOfInterest));
%     end
%   end
% end

% Go through the unitSigEpoch structure, generate a count of how many
% units, MUA had a significant run of at least 3 bins in the early, middle,
% and late phases of stimulus presentation.
% unitSigEpoch = ('Early,'Mid','Late'), ('Unsort','Unit','MUA') ,(target), ('Count', 'Run Ch Unit Index'));
unitSigEpochSoc = squeeze(unitSigEpoch(:, 2, 1, :));
MUASigEpochSoc = squeeze(unitSigEpoch(:, 3, 1, :));

% totalFractionUnits - how many units have at least 1 bout?
totalUnitswitSigRun = length(unique(vertcat(unitSigEpochSoc{:,2}), 'row'));
totalMUAswitSigRun = length(unique(vertcat(MUASigEpochSoc{:,2}), 'row'));

% Plot Results - using values found in ANOVAmat and targetRunLens
% Fig 1 - 'Social Interactions, Run Length Stats' 4 * 3 grid, 1 subplot per
% group and stat.
statTypePlotLabel = {'Number of Significant Bins', 'Length of Significant Stretches'};
statTypePlotLabel2 = {'Average Significant Stretch Length','Max Significant Stretch Length'};
statTypePlots = {statTypePlotLabel statTypePlotLabel2};
xLabelPlotArray = {'Counts', 'Run Lengths', 'Run Lengths', 'Run Lengths'};

for targ_i = 1:length(target)
  for stat_i = 1:length(statTypePlots)
    statTypePlotLabelTmp = statTypePlots{stat_i};
    figTitle = sprintf('Summary Statistic for Sliding T tests %d - %s', stat_i,  target{targ_i});
    h = figure('NumberTitle', 'off', 'Name', figTitle,'units','normalized','outerposition',[0 0 params.plotSize]);
    sgtitle(figTitle);
    for group_i = 2:length(groupingType)
      for plot_i = 1:length(statTypePlotLabelTmp)
        plotInd = (((group_i - 2) * length(statTypePlotLabelTmp))) + plot_i;
        subplot(length(groupingType)-1, length(statTypePlotLabelTmp), plotInd);
        dataInd = (((stat_i - 1) * length(statTypePlotLabelTmp))) + plot_i;
        plotData = targetRunLens{targ_i, group_i, dataInd};
        histogram(plotData);
        if plot_i == 1
          ylabel(groupingType{group_i})
        end
        xlabel(xLabelPlotArray{dataInd})
        plotTitle = sprintf('%s (90th p = %s)', statTypePlotLabelTmp{plot_i}, num2str(round(prctile(plotData, 90), 2)));
        title(plotTitle);
      end
    end
    saveFigure(params.outputDir, ['1.' figTitle], [], figStruct);
    close(h)
  end
end

% Figure 1.5 - The Population Landscape - find runs of 3+ Bins, plot
% frequen

runStruct = spikeDataBank.(runList{1});
starts = (runStruct.start:binStep:(runStruct.end - binSize))';
xlabelPlot = starts(1:8:length(starts));
stimOnInd = find(starts == 0);
stimEndInd = find(starts == 2800);

for targ_i = 1:length(target)
    figTitle = sprintf('Population Summary for Sliding T tests - %s',  target{targ_i});
    h = figure('NumberTitle', 'off', 'Name', figTitle, 'units', 'normalized', 'outerposition', [0 0 params.plotSize]);
    sgtitle(figTitle);
    for group_i = 2:length(groupingType)
      s(group_i) = subplot(2,1,group_i - 1);
      stretchStarts = targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Starts')};
      stretchLen = targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Lengths')};
      
      keepInd = stretchLen>=3;
      stretchStarts = stretchStarts(keepInd);
      stretchLen = stretchLen(keepInd);
      
      histoBins = zeros(1, 163);
      %histoBins = [];
      for stretch_i = 1:length(stretchStarts)
        histoBins(stretchStarts(stretch_i):stretchStarts(stretch_i) + stretchLen(stretch_i)) = histoBins(stretchStarts(stretch_i):stretchStarts(stretch_i) + stretchLen(stretch_i))  + 1;
        %histoBins = [histoBins, stretchStarts(stretch_i):stretchStarts(stretch_i) + stretchLen(stretch_i) - 1];
      end
      plot(histoBins)
      title(groupingType{group_i})
      xlim([0 161])
      hold on
      xlabel('Bin Start Time from Stimulus Onset (ms)')
      ylabel('Frequency (Hz)')
      xlim([0,length(histoBins)]);
      xticks(1:8:length(starts))
      xticklabels(xlabelPlot);
      if group_i == 3
        linkaxes([s(2), s(3)],'xy');
      end
      plot([stimOnInd, stimOnInd], ylim(),'k','linewidth',4);
      plot([stimEndInd, stimEndInd], ylim(),'k','linewidth',4);
    end
    saveFigure(params.outputDir, ['2.' figTitle], [], figStruct, []);
    close(h)
end

cohensStatTypePlotLabel = ["Cohen's D Values", "Maximum Cohen's D Values"];
% Figure 2 - the mean Omega Curve.
for targ_i = 1:length(target)
  figTitle = sprintf("Summary Statistic for Cohen's D - %s", target{targ_i});
  h = figure('NumberTitle', 'off', 'Name', figTitle,'units','normalized','outerposition',[0 0 params.plotSize]);
  sgtitle(figTitle);
  for group_i = 2:length(groupingType)
    for plot_i = 1:length(cohensStatType)
      plotInd = ((group_i - 2) * length(cohensStatType)) + plot_i;
      subplot(length(groupingType)-1, length(cohensStatType), plotInd);
      plotData = cohensData{targ_i, group_i, plot_i};
      plotData = [plotData(:)];
      histogram(plotData);
      if plot_i == 1
        ylabel(groupingType{group_i})
      end
      xlabel(cohensStatType{plot_i});
      plotTitle = sprintf('%s (90th p = %s)', cohensStatTypePlotLabel{plot_i}, num2str(round(prctile(plotData, 90), 2)));
      title(plotTitle);
    end
  end
  saveFigure(params.outputDir, strjoin(["2.", figTitle]), [], figStruct, []);
  close(h)
end

end

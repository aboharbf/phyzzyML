function meanPSTH(spikePathBank, params, figStruct)
% Function which combines stimulus presentations across all runs in the spikeDataBank.
% Inputs include spikeDataBank and list of parameters.
disp('Starting mean PSTH Analysis...');

plotIndParams = params.NDTParams.spikeToRasterParams.plotIndParams; % Forgive me, Father.

meanPSTHStruct = struct();

% Step 1 - Unpack variables, generate those needed for plotting.
allStimuliVec = unique(vertcat(spikePathBank.stimuli{:}));
allCategoryVec = unique(horzcat(spikePathBank.categories{:}));
groupingType = {'Unsorted', 'Unit', 'MUA'};

% Make stimPresCount for stimuli and categories

% Generate grid for indexing into individual runs and extracting relevant
% PSTHes.
runList = spikePathBank.Properties.RowNames;
small2BigInd = zeros(length(allStimuliVec), length(runList));
small2BigIndCat = zeros(length(allCategoryVec), length(runList));
for run_ind = 1:length(runList)
  [~, small2BigInd(:,run_ind)] = ismember(allStimuliVec, spikePathBank.stimuli{run_ind});
  [~, small2BigIndCat(:,run_ind)] = ismember(allCategoryVec, spikePathBank.categories{run_ind});
end

stimCatPresCounts = sum(logical(small2BigInd),2);
stimCatCounts = sum(logical(small2BigIndCat),2);

stimCatPresCounts = [stimCatPresCounts; stimCatCounts];
stimCatVec = [allStimuliVec; allCategoryVec'];
params.meanPSTHParams.topStimPresThreshold = 1;

% Step 2 - modify PSTHes - rewardEpoch removal, firstXRuns, animations,
% topStimOnly
if any([params.meanPSTHParams.plotTopStim, params.meanPSTHParams.removeFollowing, params.meanPSTHParams.stimInclude, params.meanPSTHParams.removeRewardEpoch, params.meanPSTHParams.firstXRuns])

  %Initialize a keep index
  keepInd = true(size(stimCatPresCounts));
  
  if params.meanPSTHParams.plotTopStim
    keepInd = keepInd & (stimCatPresCounts >= params.meanPSTHParams.topStimPresThreshold);
  end

  if params.meanPSTHParams.removeFollowing
    keepInd = keepInd & (~contains(stimCatVec, 'Follow'));
  end

  % Animation related processing - allows for onlyAnim, no Anim, or all.
  if params.meanPSTHParams.stimInclude ~= 0

    % If removing or focusing on anims, identify them in the list.
    animParam.stimParamsFilename = params.meanPSTHParams.stimParamsFilename;
    animParam.plotLabels = {'animSocialInteraction', 'animControl'};
    animParam.outLogic = 1;
    animParam.removeEmpty = 1;
    [tmp, ~, ~] = plotIndex(allStimuliVec, animParam);
    animInd = logical(sum(tmp,2));

    % If excluding or only including animations, update the meanPSTHStruct
    % accordingly to allow for unity between this function and those which rely
    % on its outputs.
    % Relies on stimParam file having all the things plotIndex looks at,
    % categories aren't there.
%     switch params.meanPSTHParams.stimInclude
%       case 1
%          keepInd = keepInd & animInd;
%       case 2
%          keepInd = keepInd & ~animInd;
%     end
  end
  
  % Once keep ind has been generated across all switches, edit appropriate
  % structures
  stimCatPresCounts = stimCatPresCounts(keepInd);
  stimCatVec = stimCatVec(keepInd,:);
%   stimPSTH = stimPSTH(keepInd, :, :);
%   eventMat = eventMat(keepInd, :, :);
  
  meanPSTHStruct.stimCatPresCounts = stimCatPresCounts;
  meanPSTHStruct.stimCatVec{1} = stimCatVec;
%   meanPSTHStruct.stimPSTH = stimPSTH;


  % If desired, remove the reward epoch from the base matrix, and preserve
  % the full matrix for where needed.
  if params.meanPSTHParams.removeRewardEpoch
    
    if params.meanPSTHParams.normalize
      varCheck = 'psthByImageProcNoReward';
    else
      varCheck = 'psthByImageNoReward';
    end
    
    % check to see if the reward-less PSTH exists
    if ~any(strcmp(params.spikePathLoadParams.fileVars, varCheck))
      % Find out which variables to process
      variables2Extract = {'psthByImage', 'psthErrByImage', 'psthByCategory', 'psthErrByCategory', 'psthParams'};
      if params.meanPSTHParams.normalize
        variables2Extract(1:4) = strcat(variables2Extract(1:4), 'Proc');
      else
        
      end
      variables2Save(1:4) = strcat(variables2Extract(1:4), 'NoReward');
      
      % Collect the vectors of activity from each analysis directory.
      [psthByImage, psthErrByImage, psthByCategory, psthErrByCategory, psthParamsPerRun] = spikePathLoad(spikePathBank, variables2Extract, params.spikePathLoadParams);
      analyzedDataBatchPaths = fullfile(spikePathBank.analyzedDir, params.spikePathLoadParams.batchAnalysisOutputName);
      
      for run_i = 1:length(runList)
        % Identify the period of the activity to keep.
        nonRewardPeriod = 1:1+(psthParamsPerRun{run_i}.psthPre+psthParamsPerRun{run_i}.psthImDur);
        
        % Initialize new vectors for shorter vectors
        [A, B] = deal(initNestedCellArray(psthByImage{run_i}, 'NaN'));
        [C, D] = deal(initNestedCellArray(psthByCategory{run_i}, 'NaN'));
        
        % Cycle through and shorten the vectors.
        for chan_i = 1:length(psthByImage{run_i})
          for unit_i = 1:length(psthByImage{run_i}{chan_i})
            A{chan_i}{unit_i} = psthByImage{run_i}{chan_i}{unit_i}(:, nonRewardPeriod); 
            B{chan_i}{unit_i} = psthErrByImage{run_i}{chan_i}{unit_i}(:, nonRewardPeriod);             
            C{chan_i}{unit_i} = psthByCategory{run_i}{chan_i}{unit_i}(:, nonRewardPeriod); 
            D{chan_i}{unit_i} = psthErrByCategory{run_i}{chan_i}{unit_i}(:, nonRewardPeriod); 
          end
        end
        
        tmpVar = {'A', 'B', 'C', 'D'};
        
        % Save the new variables to the output file in each directory
        for var_i = 1:length(variables2Save)
          eval(sprintf('%s = %s;', variables2Save{var_i}, tmpVar{var_i}))
          if ~exist(analyzedDataBatchPaths{run_i})
            save(analyzedDataBatchPaths{run_i}, variables2Save{var_i})
          else
            save(analyzedDataBatchPaths{run_i}, variables2Save{var_i}, '-append')
          end
        end
      end
      
      % Update the param files, save to the original.
      if ~any(strcmp(params.spikePathLoadParams.fileVars, varCheck))
        params.spikePathLoadParams.files = [params.spikePathLoadParams.files, params.spikePathLoadParams.batchAnalysisOutputName];
      end
      params.spikePathLoadParams.fileVars = [params.spikePathLoadParams.fileVars, variables2Save];
      params.spikePathLoadParams.fileInds = [params.spikePathLoadParams.fileInds, repmat(4, [1,length(variables2Save)])];
      batchAnalysisParams = params;
      spikePathFile = batchAnalysisParams.spikePathLoadParams.batchAnalysisOutput;
      save(spikePathFile, 'batchAnalysisParams', '-append')

    end
  end

  % If you want to only count the first X runs of a stimulus, 
%   if params.meanPSTHParams.firstXRuns
%     keepMat = cellfun(@(x) x <= params.meanPSTHParams.firstXRuns, stimPSTH(:, :, strcmp(meanPSTHStruct.IndStructs{3}, 'presCount')), 'UniformOutput', 0);
% 
%     for ii = 1:size(stimPSTH, 1)
%       for jj = 1:size(stimPSTH, 2)
%         for kk = 1:size(stimPSTH, 3)
%           stimPSTH{ii, jj, kk} = stimPSTH{ii, jj, kk}(keepMat{ii, jj},:);
% 
%           if jj == 3 && kk == 3 %MUA and presCount
%             % Update stimPresCounts according to the new numbers.
%             stimCatPresCounts(ii) = max(stimPSTH{ii, jj, kk});
%           end
% 
%         end
%       end
%     end
%     
%   end

end

% Generate the figure directory
dirTags = {' NoReward', 'fixAligned', 'Normalized', 'Pres Threshold', 'AnimOnly', 'NoAnim', sprintf('first%d', params.meanPSTHParams.firstXRuns)};
dirTagSwitch = [params.meanPSTHParams.removeRewardEpoch, params.meanPSTHParams.fixAlign, params.meanPSTHParams.normalize, params.meanPSTHParams.plotTopStim, params.meanPSTHParams.stimInclude == 1, params.meanPSTHParams.stimInclude == 2, params.meanPSTHParams.firstXRuns];
params.meanPSTHParams.outputDir = [params.meanPSTHParams.outputDir, strjoin(dirTags(logical(dirTagSwitch)),' - ')];

if ~exist(params.meanPSTHParams.outputDir, 'dir')
  mkdir(params.meanPSTHParams.outputDir)
end

% Plotting related variables
eventColors = 'kbrg';
groupIterInd = 2:3; % An index for which groups to look at. 1 = Unsorted, 2 = Units, 3 = MUA. Combine as fit.

% [plotMat, briefStimList, params.meanPSTHParams] = plotIndex(allStimuliVec, params.meanPSTHParams);
% 
if params.meanPSTHParams.normalize
  normTag = ' - normalized';
else
  normTag = '';
end
% 
% allStimuliNames = cellfun(@(x) extractBetween(x, 1, length(x)-4), allStimuliVec);
% allStimuliNames = strrep(allStimuliNames, '_', ' ');
% 
% stimPresMat = cellfun(@(x) size(x,1),stimPSTH(:,:,end));
% meanPSTHStruct.stimPresMat = stimPresMat;
% broadLabelInd = plotMat;
% dataInd2Plot = strcmp(meanPSTHStruct.IndStructs{3}, 'PSTH') | strcmp(meanPSTHStruct.IndStructs{3}, 'PSTH Err');

% Variable preperation - extract variables depending on switches
% variables2Extract = {'psthByImage', 'psthErrByImage', 'psthByCategory', 'psthErrByCategory', 'psthParams'};
variables2Extract = {'psthByImage', 'psthByCategory', 'psthParams'};
if params.meanPSTHParams.normalize
  variables2Extract(1:2) = strcat(variables2Extract(1:2), 'Proc');
end

if params.meanPSTHParams.removeRewardEpoch
  variables2Extract(1:2) = strcat(variables2Extract(1:2), 'NoReward');
end

variables2Extract = [variables2Extract, 'selTable'];

% [psthByImage, psthErrByImage, psthByCategory, psthErrByCategory, psthParamsPerRun] = spikePathLoad(spikePathBank, variables2Extract, params.spikePathLoadParams);
[psthByImage, psthByCategory, psthParamsPerRun, selTable] = spikePathLoad(spikePathBank, variables2Extract, params.spikePathLoadParams);
paradigmList = unique(spikePathBank.paradigmInd);
runList = extractAfter(spikePathBank.Properties.RowNames, 'S');

% Plot 1 - All Stimuli means in the same plot & Plot 2 - Catagory Plot - 'All Chasing Stimuli, mean PSTH'
if params.meanPSTHParams.allStimPSTH || params.meanPSTHParams.catPSTH || params.meanPSTHParams.analysisGroupPSTH
  
  % Create a few strings for later, for the sake of calling eval instead of
  % writing the same code many times.
  psth2Collect = {'psthByImage', 'psthByCategory', 'psthByCategory'};
  spikePathBankField = {'stimuli', 'categories', 'categories'};
  spikePathBankPlot = {'Stimuli', 'Categories', 'AnalysisGroups'};
  
  selectionInd = logical([params.meanPSTHParams.allStimPSTH params.meanPSTHParams.catPSTH params.meanPSTHParams.analysisGroupPSTH]);
  psth2Collect = psth2Collect(selectionInd);
  spikePathBankField = spikePathBankField(selectionInd);
  spikePathBankPlot = spikePathBankPlot(selectionInd);
  plotNum = find(selectionInd);
  
  % For every paradigm
  for ii = 1:length(paradigmList)
    
    % Cycle through the run data, collecting the traces of interest
    paradigmInd = spikePathBank.paradigmInd == paradigmList(ii);
    paradigmName = spikePathBank.paradigmName{find(paradigmInd,1)};
    psthParamTimes = psthParamsPerRun{find(paradigmInd,1)};
    if params.meanPSTHParams.removeRewardEpoch
      psthParamTimes.psthPost = 0;
    end
    
    spikePathBankParadigm = spikePathBank(paradigmInd, :);
    
    % check to see if analysisGroups need ot be plotted as well
    if params.meanPSTHParams.analysisGroupPSTH
      analysisGroups2Plot = [paradigmList(ii); fields(params.meanPSTHParams.analysisGroups.(paradigmName))];
    else
      analysisGroups2Plot = {paradigmList(ii)};
    end
    
    for plot_i = 1:length(analysisGroups2Plot)
      for data_i = 1:length(psth2Collect)
        
        traces2Extract = eval(sprintf('%s(paradigmInd)', psth2Collect{data_i})); % way to extract Image or Category data.
        if isnumeric(analysisGroups2Plot{plot_i})
          [meanTraces, errTraces, counts, paradigmStimPlot] = extractMeanTraces(traces2Extract, spikePathBankParadigm.(spikePathBankField{data_i}), [], 0);
        else
          [meanTraces, errTraces, counts, paradigmStimPlot] = extractMeanTraces(traces2Extract, spikePathBankParadigm.(spikePathBankField{data_i}), params.meanPSTHParams.analysisGroups.(paradigmName).(analysisGroups2Plot{plot_i}), 0);
        end
        
        % check for empty cells, as this messes w/ later plotting
        if any(cellfun('isempty', meanTraces))
          keepInd = ~cellfun('isempty', meanTraces(:,1));
          meanTraces = meanTraces(keepInd, :);
          errTraces = errTraces(keepInd, :);
          counts = counts(keepInd, :);
          paradigmStimPlot = paradigmStimPlot(keepInd, :);
        end
        
        if isempty(meanTraces)
          continue;
        end
        
        for group_ind = groupIterInd
          
          % pull and label a few variables.
          allStimuliNames = paradigmStimPlot;
          plotData = vertcat(meanTraces{:,group_ind});
          
          % If there is no data to plot, skip
          if all(all(isnan(plotData)))
            continue
          end
          
          plotErr = vertcat(errTraces{:,group_ind});
          stimCounts = counts(:,group_ind);
          
          % If you want to add the count to each title, do so here.
          if params.meanPSTHParams.traceCountLabel
            allStimuliLabel = cell(length(allStimuliNames),1);
            for stim_ind = 1:length(allStimuliNames)
              allStimuliLabel{stim_ind} = [allStimuliNames{stim_ind} ', n = ' num2str(stimCounts(stim_ind))];
            end
          else
            allStimuliLabel = allStimuliNames;
          end
          
          % If you want to sort based on presentation number, do so here.
          if params.meanPSTHParams.sortPresCountSort
            [~, newOrder] = sort(stimCounts);
            allStimuliLabel = allStimuliLabel(newOrder);
            plotData = plotData(newOrder, :);
            plotErr = plotErr(newOrder, :);
          end
          
          % Set the figure title
          if isnumeric(analysisGroups2Plot{plot_i})
            catPSTHTitle = sprintf('Paradigm %s - All Catagories, Mean PSTH %s, %s', paradigmName, normTag, groupingType{group_ind});
          else
            catPSTHTitle = sprintf('Paradigm %s - %s, Mean PSTH, paradigm %s %s, %s', paradigmName, analysisGroups2Plot{plot_i}, normTag, groupingType{group_ind});
          end
          
          h = figure('NumberTitle', 'off', 'Name', catPSTHTitle,'units','normalized','outerposition',[0 0 params.meanPSTHParams.plotSizeAllStimPSTH]);
          
          % Things w/ only a few traces can be plotted as lines.
          %           if size(plotData,1) > 8
          [~, cbh] = plotPSTH(plotData, plotErr, axes(), psthParamTimes, 'color', catPSTHTitle, allStimuliLabel);
          %           else
          %             [~, cbh] = plotPSTH(plotData, plotErr, axes(), psthParamTimes, 'line', catPSTHTitle, allStimuliLabel);
          %           end
          
          % Change colorbar axis label
          if params.meanPSTHParams.normalize == 1
            cbh.Label.String = 'Signal Change relative to Baseline (%)';
          elseif params.meanPSTHParams.normalize == 2
            cbh.Label.String = 'Z scored relative to fixation';
            ylabel('Normalized Activity (Baseline Z scored)');
          else
            ylabel('Activity (Firing Rate)');
          end
          
          % Add the title, some more labels, and save
          title(catPSTHTitle);
          xlabel('Time from Stimulus Onset (ms)');
          saveFigure(params.meanPSTHParams.outputDir, sprintf('%d. %s', plotNum(data_i), catPSTHTitle), [], figStruct, [])
          clear allStimuliLabel
        end
      end
    end
  end
end

% Plot 2 + 3 - analysisGroup Plot - 'All head Turning vs non-Head turning'
if params.meanPSTHParams.analysisGroupPSTH
  % Extract the groups to plot per paradigm.
  
  for ii = 1:length(paradigmList)
    % Cycle through the run data, collecting the traces of interest
    
    paradigmInd = spikePathBank.paradigmInd == paradigmList(ii);
    paradigmName = spikePathBank.paradigmName{find(paradigmInd,1)};
    psthParamTimes = psthParamsPerRun{find(paradigmInd,1)};
    spikePathBankParadigm = spikePathBank(paradigmInd, :);
    
    for data_i = 1:length(psth2Collect)
      
      traces2Extract = eval(sprintf('%s(paradigmInd)', psth2Collect{data_i})); % way to extract Image or Category data.
      [meanTraces, ~, counts, paradigmStimPlot] = extractMeanTraces(traces2Extract, spikePathBankParadigm.(spikePathBankField{data_i}), [], 1);
      
      for group_ind = groupIterInd
        
        allStimuliNames = paradigmStimPlot;
        plotData = vertcat(meanTraces{:,group_ind});
        stimCounts = counts(:,group_ind);
        
        if params.meanPSTHParams.traceCountLabel
          allStimuliLabel = cell(length(allStimuliNames),1);
          for stim_ind = 1:length(allStimuliNames)
            allStimuliLabel{stim_ind} = [allStimuliNames{stim_ind} ', n = ' num2str(stimCounts(stim_ind))];
          end
        else
          allStimuliLabel = allStimuliNames;
        end
        
        if params.meanPSTHParams.sortPresCountSort
          [~, newOrder] = sort(stimCounts);
          allStimuliLabel = allStimuliLabel(newOrder);
          plotData = plotData(newOrder, :);
        end
        
        catPSTHTitle = sprintf('2. All %s, Mean PSTH for %s %s, %s', spikePathBankPlot{data_i}, paradigmName, normTag, groupingType{group_ind});
        h = figure('NumberTitle', 'off', 'Name', catPSTHTitle,'units','normalized','outerposition',[0 0 params.meanPSTHParams.plotSizeAllStimPSTH]);
        [~, cbh] = plotPSTH(plotData, [], axes(), psthParamTimes, 'color', catPSTHTitle, allStimuliLabel);
        if params.meanPSTHParams.normalize == 1
          cbh.Label.String = 'Signal Change relative to Baseline (%)';
        elseif params.meanPSTHParams.normalize == 2
          cbh.Label.String = 'Z scored relative to fixation';
          ylabel('Normalized Activity (Baseline Z scored)');
        else
          ylabel('Activity (Firing Rate)');
        end
        title(catPSTHTitle);
        xlabel('Time from Stimulus Onset (ms)');
        saveFigure(params.meanPSTHParams.outputDir, sprintf('%d. %s', data_i, catPSTHTitle), [], figStruct, [])
        clear allStimuliLabel
      end
    end
    
  end
    
end

params.comboEvents = {'subSel_headTurn_all', 'subSel_allTurn', 'socIntSel_any', 'headTurnSel_any'};
params.comboSubEvents = {{'subSel_headTurn_left', 'subSel_headTurn_right'}, {'subSel_headTurn_left', 'subSel_headTurn_right', 'subSel_bodyTurn'}, ...
  {'socIntSel_stimOnset', 'socIntSel_stimPres', 'socIntSel_reward'}, {'headTurn_baseDiff_stimPres', 'headTurn_baseDiff_reward'}};

% Plot 4 - Create a plot with means for selective units between conditions, and non-selective units.
if params.meanPSTHParams.selPSTH
  UnitTypes = params.selParam.UnitTypes;
  UnitTypePlot = params.selParam.UnitTypePlot;
  stimParamData = load(params.stimParamsFilename);
  stimParamStimuli = cellfun(@(x) x(1), stimParamData.paramArray);
  
  for ii = 1:length(paradigmList)
    
    % Cycle through the run data, collecting the traces of interest
    paradigmInd = spikePathBank.paradigmInd == paradigmList(ii);
    spikePathBankParadigm = spikePathBank(paradigmInd, :);
    paradigmName = spikePathBankParadigm.paradigmName{1};
    psthParamTimes = psthParamsPerRun{find(paradigmInd,1)};
    paradigmSelTable = selTable(paradigmInd);
    psthByImageParadigm = psthByImage(paradigmInd);
    runListParadigm = runList(paradigmInd);
    plotIndParams.plotLabels = {{'agents', 'socialInteraction'}};
    
    
    % Combine across paradigmSelTable
    paradigmSelTable = vertcat(paradigmSelTable{:});
    
    % Replace channel numbers for this annoying run
    paradigmSelTable = replaceChanNum_1123(paradigmSelTable);

    paradigmSelTable = expandSelTableComboEvents(paradigmSelTable, params.selParam);
    variableList = paradigmSelTable.Properties.VariableNames(contains(paradigmSelTable.Properties.VariableNames', 'socIntSel'));
    
    % Create an array of selectivity type, unitType, Sign, and dataType
    dataType = {'agents', 'socialInteraction'};
    normalizedActivity = cell(length(variableList), length(UnitTypes), 2, length(dataType));
    [normalizedActivity{:}] = deal(nan(1e3 ,size(psthByImageParadigm{1}{1}{1},2)));
    
    % For every event, find the unit traces selective for it.
    totCount = 1;
    for event_i = 1:length(variableList)
      for sign_i = 1:2
        
        % Seperate Increases and decresed in firing
        if sign_i == 1
          eventSelective = paradigmSelTable.(variableList{event_i}) > 0;
        else
          eventSelective = paradigmSelTable.(variableList{event_i}) < 0;
        end
        paradigmSelTableEvent = paradigmSelTable(eventSelective,:);
        
        % Cycle for normal units and MUA.
        for unitType_i = 1:length(UnitTypes)
          paradigmSelTableEventUnit = paradigmSelTableEvent(contains(paradigmSelTableEvent.unitType, UnitTypes{unitType_i}), :);
          runListUnit = strcat(paradigmSelTableEventUnit.dateSubj, paradigmSelTableEventUnit.runNum);
          
          chanUnit = double(extractAfter(paradigmSelTableEventUnit.channel, 'Ch'));
          unitUnit = paradigmSelTableEventUnit.unitType;
          
          % for each unit, Grab the Activity traces, and store
          for unit_i = 1:size(paradigmSelTableEventUnit, 1)
            if strcmp(unitUnit(unit_i), 'MUA')
              psthRun = psthByImageParadigm{strcmp(runListParadigm, runListUnit{unit_i})}{chanUnit(unit_i)}{end};
            else
              unit_Index = double(extractAfter(paradigmSelTableEventUnit.unitType(unit_i), 'U')) + 1;
              psthRun = psthByImageParadigm{strcmp(runListParadigm, runListUnit{unit_i})}{chanUnit(unit_i)}{unit_Index};
            end
            
            % Grab the traces, and the traces they belong to.
            runStim = spikePathBankParadigm.stimuli{strcmp(runListParadigm, runListUnit{unit_i})};
            traceIndexPerStimuli = plotIndex(runStim, plotIndParams);
            
            % store traces
            for trace_i = 1:2
              nanInd = find(isnan(normalizedActivity{event_i, unitType_i, sign_i, trace_i}(:,1)), 1);
              traces2Store = traceIndexPerStimuli == trace_i;
              
              normalizedActivity{event_i, unitType_i, sign_i, trace_i}(nanInd:nanInd+sum(traces2Store)-1, :) = psthRun(traces2Store,:);
              fprintf('KEEP ON GOIN: %d /n\n', totCount);
              totCount = totCount + 1;
            end
            
          end
        end
      end
    end
    
    % After retrieving all the traces, remove the NaNs from the stored
    % values
    for pp = 1:size(normalizedActivity, 1)
      for jj = 1:size(normalizedActivity, 2)
        for kk = 1:size(normalizedActivity, 3)
          for ll = 1:size(normalizedActivity, 4)
            NaNInd = find(isnan(normalizedActivity{pp, jj, kk, ll}(:,1)));
            
            normalizedActivity{pp, jj, kk, ll} = normalizedActivity{pp, jj, kk, ll}(1:NaNInd-1, :);
          end
        end
      end
    end
    
    % Plot Time - % normalizedActivity - selectivity type, unitType, Sign, and dataType
    signPlot = {'Increased Activity', 'Decreased Activity'};
    variableListPlot = strrep(variableList, '_', ' ');
    legendPlots = {'Non-SI Agents', 'SI Agents'};
    for event_i = 1:length(variableList)
      for unitType_i = 1:length(UnitTypePlot)
        for sign_i = 1:2
          
          % Extract 2 traces
          traces2Plot = squeeze(normalizedActivity(event_i, unitType_i, sign_i, :));
          
          [meanTraces, stdTraces] = deal(zeros(2, size(traces2Plot{1}, 2)));
          % Take means across each
          for jj = 1:length(traces2Plot)
            meanTraces(jj,:) = mean(traces2Plot{jj}, 1);
            stdTraces(jj,:) = std(traces2Plot{jj}, 0, 1)./sqrt(size(traces2Plot{jj},1));
          end         
          
          % Plot the Figure
          psthTitle = sprintf('Normalized Population Activity of %s selective %s (%s)', variableListPlot{event_i}, UnitTypePlot{unitType_i}, signPlot{sign_i});
          psthParamTimes.psthPost = 0;
          h = figure('Name', psthTitle, 'units', 'normalized', 'outerposition', [0 0 params.meanPSTHParams.plotSizeLineCatPlot]);
          plotPSTH(meanTraces, stdTraces, [], psthParamTimes, 'line', psthTitle, legendPlots)
          ylabel('Change in Activity relative to Baseline, Normalized')
          h.Children(2).FontSize = 16;
          saveFigure(params.meanPSTHParams.outputDir, ['4. ' psthTitle], [], figStruct, [])
          
        end
      end
    end
    
  end
end

end

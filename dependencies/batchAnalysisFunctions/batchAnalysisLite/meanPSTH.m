function meanPSTH(spikePathBank, params, figStruct)
% Function which combines stimulus presentations across all runs in the spikeDataBank.
% Inputs include spikeDataBank and list of parameters.
disp('Starting mean PSTH Analysis...');

meanPSTHStruct = struct();

% Step 1 - Unpack variables, generate those needed for plotting.
allStimuliVec = unique(vertcat(spikePathBank.stimuli{:}));
allCategoryVec = unique(vertcat(spikePathBank.categories{:}));
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
stimCatVec = [allStimuliVec; allCategoryVec];
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

% [psthByImage, psthErrByImage, psthByCategory, psthErrByCategory, psthParamsPerRun] = spikePathLoad(spikePathBank, variables2Extract, params.spikePathLoadParams);
[psthByImage, psthByCategory, psthParamsPerRun] = spikePathLoad(spikePathBank, variables2Extract, params.spikePathLoadParams);
paradigmList = unique(spikePathBank.paradigmInd);

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

% Plot 3 - analysisGroup Plot - 'All head Turning vs non-Head turning'
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
        
        catPSTHTitle = sprintf('All %s, Mean PSTH for %s %s, %s', spikePathBankPlot{data_i}, paradigmName, normTag, groupingType{group_ind});
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

% % Plot 3 - Stimuli Plot - 'All chasing 1 PSTHs, sorted by...'
% if params.meanPSTHParams.allRunStimPSTH 
%   % Make a PSTH of each stimulus across all its repetitions.
%   sortType = {'Run of Day', 'Grid Hole', 'Recording Depth', 'Run Ind'};
%   % {'Run of Day', 'Grid Hole', 'Recording Depth', 'Run Ind'};
%   % Must be the same order as the final parts of stimPSTH
%   [sortMat, sortLabel] = deal(cell(size(stimPSTH,1), size(stimPSTH,2), length(sortType)));
%   
%   % Generate the PSTH sorting indicies
%   for stim_i = 1:size(stimPSTH,1)
%     for group_i = 1:size(stimPSTH,2)
%       for sort_i = 1:length(sortType) % the sortings that need processing
%         
%         switch sort_i
%           case {1,2}
%           % Store indicies in 1st, labels in 2nd %Grid Holes --> Indicies
%           [tmpLabels , sortMat{stim_i, group_i, sort_i}] = sort(stimPSTH{stim_i, group_i,strcmp(dataType, sortType{sort_i})});
%           
%           % Generate Labeling index
%           uniqueLabels = unique(tmpLabels);
%           labeledIndex = zeros(length(uniqueLabels),1);
%           for uni_i = 1:length(uniqueLabels)
%             if sort_i == 1
%               labeledIndex(uni_i) = find(tmpLabels == uniqueLabels(uni_i), 1);
%             elseif sort_i == 2
%               labeledIndex(uni_i) = find(strcmp(tmpLabels, uniqueLabels(uni_i)),1);
%             end
%           end
%           case 3 %Recording Depth --> Indicies
%           tmpDepths = [stimPSTH{stim_i, group_i,strcmp(dataType, sortType{sort_i})}];
%           [tmpLabels , sortMat{stim_i, group_i, sort_i}] = sort(tmpDepths');
%           labeledIndex = 1:5:length(tmpLabels);
%           case 4
%           sortMat{stim_i, group_i, sort_i} = 1:length(stimPSTH{stim_i, group_i, strcmp(dataType, sortType{sort_i})});
%           tmpLabels = runList(stimPSTH{stim_i, group_i,strcmp(dataType, sortType{sort_i})});
%           labeledIndex = 1:5:length(tmpLabels);
%         end
%         
%         % Use the labeledIndex to generate the proper label array with 0's
%         % everywhere else.
%         sortLabelTmp = cell(length(tmpLabels),1);
%         for tmp_i = 1:length(labeledIndex)
%           if sort_i == 2 || sort_i == 4
%             sortLabelTmp{labeledIndex(tmp_i)} = tmpLabels{labeledIndex(tmp_i)};
%           else
%             sortLabelTmp{labeledIndex(tmp_i)} = tmpLabels(labeledIndex(tmp_i));
%           end
%         end
%         sortLabel{stim_i, group_i, sort_i} = tmpLabels;
%         
%       end
%     end
%   end
%   
%   % Iterate through PSTH, generating plots
%   if params.allRunStimPSTH
%     for stim_i = 1:length(stimPSTH)
%       % Get the relevant frameMotion/eventData
%       eventDataStim = eventData(allStimuliVec{stim_i},:);
%       stimTimePerFrame = frameMotionData(strcmp(allStimuliVec{stim_i}, frameMotionDataNames)).timePerFrame;
%       for group_i = groupIterInd
%         stimData = stimPSTH{stim_i,group_i,1};
%         
%         for sort_i = 4%1:length(sortType)
%           
%           figTitle = sprintf('%s - %s PSTHs, Sorted by %s, %s', allStimuliNames{stim_i}, groupingType{group_i} ,sortType{sort_i}, normTag);
%           h = figure('NumberTitle', 'off', 'Name', figTitle,'units','normalized','outerposition',[0 0 params.plotSizeAllRunStimPSTH]);
%           sgtitle(figTitle)
%           
%           sortIndex = sortMat{stim_i, group_i, sort_i};
%           
%           % Break up PSTHes with many lines into subplots.
%           TracesPerPlot = ceil(size(sortIndex,2)/3);
%           if TracesPerPlot < 20
%             subplot2Plot = 1;
%           elseif TracesPerPlot < 45
%             subplot2Plot = 2;
%           else
%             subplot2Plot = 3;
%           end
%           TracesPerPlot = ceil(size(sortIndex,2)/subplot2Plot);
%           
%           plotStarts = 1:TracesPerPlot:size(sortIndex,2);
%           plotEnds = [plotStarts(2:end)-1, size(sortIndex,2)];
%           
%           subplotAxes = gobjects(subplot2Plot,1);
%           cbHandle = gobjects(subplot2Plot,1);
%           for plot_i = 1:subplot2Plot
%             plotLabels = sortLabel{stim_i, group_i, sort_i}(plotStarts(plot_i):plotEnds(plot_i));
%             plotData = stimData(plotStarts(plot_i):plotEnds(plot_i), :);
%             %sortIndex = sortMat{stim_i, group_i, sort_i}(plotStarts(plot_i):plotEnds(plot_i));
%             psthAxes = subplot(1,subplot2Plot,plot_i);
%             [subplotAxes(plot_i), cbHandle(plot_i)] = plotPSTH(plotData, [], psthAxes, params, 'color', [], plotLabels); %(sortIndex,:)
%             cbHandle(plot_i).Label.FontSize = 10;
%             
%             % Add eventData line if applicable
%             [eventsPres, legendObjs] = deal([]);
%             
%             for event_i = 1:length(eventList)
%               if ~isempty(eventDataStim.(eventList{event_i}){1})
%                 eventsPres = [eventsPres; eventList(event_i)];
%                 hold on
%                 singleEventDataStim = eventDataStim.(eventList{event_i}){1};
%                 for ev_i = 1:size(singleEventDataStim, 1)
%                   startX = singleEventDataStim.startFrame(ev_i) * stimTimePerFrame;
%                   endX = singleEventDataStim.endFrame(ev_i) * stimTimePerFrame;
%                   lineLeg = plot([startX startX], ylim(), 'Color', eventColors(event_i), 'LineWidth', 2);
%                   plot([endX endX], ylim(), 'Color', eventColors(event_i), 'LineWidth', 2);
%                   if ev_i == 1
%                     legendObjs = [legendObjs; lineLeg];
%                   end
%                 end
%               end
%             end
%             
%             if plot_i == subplot2Plot
%               % Add Legend for event lines, if present
%               if ~isempty(legendObjs)
%                 legend(legendObjs, eventsPres, 'location', 'northeastoutside', 'Fontsize', 8);
%               end
%               % Label the Y
%               if params.normalize == 1
%                 cbHandle(plot_i).Label.String = 'Signal Change relative to Baseline (%)';
%               elseif params.normalize == 2
%                 cbHandle(plot_i).Label.String = 'Z score relative to fixation';
%               end
%             else
%               delete(cbHandle(plot_i))
%             end
%             
%             set(gca,'FontSize',10,'TickLength',[.01 .01],'LineWidth',.25);
%           end
%           linkprop(subplotAxes, 'CLim');
%           
%           saveFigure(params.outputDir, ['3. ' figTitle], [], figStruct, []);          
%         end
%       end
%     end
%   end
% end
% 
% % Plot 4 - Line plot with Line per Catagory
% if params.lineCatPlot
%   % Find the end of the 'stimuli catagory' count.
%   catCount = find(strcmp(params.plotLabels, 'scene'));
%   if isempty(catCount)
%     catCount = find(strcmp(params.plotLabels, 'grooming'));
%   end
%   singleCatPlotMat = plotMat(:,1:catCount);
%   singleCatplotLabels = params.plotLabels(1:catCount);
%   singleCatplotLabelsSocialInd = params.plotLabelSocialInd(1:catCount);
%   %tmp = distinguishable_colors(catCount);
%   plotColors = cell(catCount,1);
%   socCount = sum(singleCatplotLabelsSocialInd);
%   nonSocCount = sum(~singleCatplotLabelsSocialInd);
%   
%   socialColorsHSV = repmat(rgb2hsv(params.socialColor), [socCount, 1]);
%   nonSocialColorsHSV = repmat(rgb2hsv(params.nonSocialColor), [nonSocCount, 1]);
%   % HSV values cap at 240
%   hsvGradient = 25/240;
%   hsvGradientSoc = hsvGradient * (1:socCount);
%   hsvGradientNonSoc = hsvGradient * (1:nonSocCount);
%   
%   socialColorsHSV(:,3) = socialColorsHSV(:,3) - hsvGradientSoc';
%   nonSocialColorsHSV(:,3) = nonSocialColorsHSV(:,3) - hsvGradientNonSoc';
%   plotColorTmp = [socialColorsHSV; nonSocialColorsHSV];
%   plotColorTmp = hsv2rgb(plotColorTmp);
%   
%   for col_i = 1:catCount
%     plotColors{col_i} = plotColorTmp(col_i,:);
%   end
%   
%   for group_ind = groupIterInd
%     catPSTHTitle = sprintf('%s Mean PSTH per Catagory Label %s', groupingType{group_ind}, normTag);
%     h = figure('NumberTitle', 'off', 'Name', catPSTHTitle,'units','normalized','outerposition',[0 0 params.plotSizeLineCatPlot]);
%     hold on
%     h.Children.FontSize = 15;
%     % Generate line plots w/ error bars.
%     lineProps.width = 2;
%     lineProps.col = plotColors;
%     lineProps.patch.FaceAlpha = '0.5';
%     
%     [groupData, groupErr] = deal(cell(size(singleCatPlotMat,2),1));
%     for cat_ind = 1:size(singleCatPlotMat,2)
%       tmp = vertcat(stimPSTH{logical(singleCatPlotMat(:,cat_ind)), group_ind});
%       groupData{cat_ind} = mean(tmp,1);
%       groupErr{cat_ind} = std(tmp)/sqrt(size(tmp,1));
%     end
%     
%     groupData = vertcat(groupData{:});
%     
%     if params.fixAlign
%       groupMean = mean(mean(groupData(:,500:800)));
%       for g_ind = 1:size(groupData,1)
%         groupData(g_ind,:) = groupData(g_ind,:) - mean(groupData(g_ind,500:800)) + groupMean;
%       end
%     end
%     
%     mseb(-800:(size(groupData,2)-801),groupData, vertcat(groupErr{:}), lineProps);
%     xlim([-800 (size(groupData,2)-801)])
%     legend(singleCatplotLabels, 'AutoUpdate','off','location', 'northeastoutside');
%     ylim manual
%     line([0 0], ylim(), 'Linewidth', 3, 'color', 'k');
%     line([2800 2800], ylim(), 'Linewidth',3,'color','k');
%     xlabel('Time from Stimulus Onset (ms)');
%     ylabel('Normalized Activity (Baseline Z scored)');
%     title(catPSTHTitle);
%     saveFigure(params.outputDir, ['4. ' catPSTHTitle], [], figStruct, []);
%   end
%   
% end
% 
% % Plot 5 - Means across broad catagorizations (like Social vs non Social)
% if params.lineBroadCatPlot  
%   % Generate line plots across different labeling schemes.
%   agentInd = plotMat(:,strcmp(params.plotLabels,'agents'));
%   socialInd = plotMat(:,strcmp(params.plotLabels,'socialInteraction'));
%   headTurnInd = plotMat(:,strcmp(params.plotLabels,'headTurn'));
%   allTurnInd = logical(plotMat(:,strcmp(params.plotLabels,'allTurn')));
%   %headTurnOldInd = logical(plotMat(:,strcmp(params.plotLabels,'headTurnClassic')));
% 
%   % Agent videos without head turning.
%   agentNHInd = agentInd & ~headTurnInd;
%   agentNTInd = agentInd & ~allTurnInd;
%   
%   % Social videos 
%   socialHTInd = socialInd & headTurnInd;        %   Head turning vs Not
%   socialNonHTInd = socialInd & ~headTurnInd;    %   
%   socialTInd = socialInd & allTurnInd;
%   socialNonTInd = socialInd & ~allTurnInd;
%       
%   % Social Agent vs Non-Social Agent
%   socialAgentInd = socialInd;
%   nonSocialAgentInd = agentInd & ~socialInd;
%   
%   figIncInd = {agentInd, socialInd, headTurnInd, allTurnInd, socialHTInd, socialAgentInd, socialTInd};
%   figExcInd = {~agentInd, ~socialInd, agentNHInd, agentNTInd, socialNonHTInd, nonSocialAgentInd, socialNonTInd};
%   figTitInd = {'Agent containing Stimuli Contrast',...
%     'Social Interactions Contrast',...
%     'Agents engaging in Head Turning Contrast',...
%     'Agents engaging in All Turning Contrast',...
%     'Social Interactions with Head turning Contrast',...
%     'Agents engaging in Social Interactions Contrast'...
%     'Social Interactions with all Turning Contrast'};
%   figLegends = {{'Agent containing Stimuli','Non-Agent containing Stimuli'},...
%     {'Social Interaction Stimuli','non-Social Interaction Stimuli'},...
%     {'Agents with Head Turning','Agents without Head Turning'},...
%     {'Agents with Any Turning','Agents without Any Turning'},...
%     {'Socially Interacting Agents with Head Turning','Socially Interacting Agents without Head Turning'},...
%     {'Agents engaging in Social Interactions ','Agents not engaging in Social Interactions'}...
%     {'Socially Interacting Agents with Any Turning','Socially Interacting Agents without Any Turning'}};
%   lineProps.col = {params.socialColor, params.nonSocialColor};
%   assert(length(figIncInd) == length(figLegends), 'Update figTitInd and figLegends to match figIncInd')
%   
%   meanPSTHStruct.lineBroadCatPlot.figIncInd = figIncInd;
%   meanPSTHStruct.lineBroadCatPlot.figExcInd = figExcInd;
%   meanPSTHStruct.lineBroadCatPlot.figTitInd = figTitInd;
%   meanPSTHStruct.lineBroadCatPlot.figLegends = figLegends;
%   meanPSTHStruct.lineBroadCatPlot.lineprops = lineProps;
%   
%   % Append counts to the stimuliNames
%   allStimuliLabels = arrayfun(@(x) sprintf('%s (n = %d)', allStimuliNames{x}, stimCatPresCounts(x)), 1:length(stimCatPresCounts), 'UniformOutput',0)';
%   
%   for fig_ind = 1:length(figIncInd)
%     % Retrieve the labels for each plot
%     line1Array = stimPSTH(figIncInd{fig_ind},:,1);
%     line2Array = stimPSTH(figExcInd{fig_ind},:,1);
%     line1Labels = allStimuliLabels(figIncInd{fig_ind});
%     line2Labels = allStimuliLabels(figExcInd{fig_ind});
%     line1EventMat = eventMat(figIncInd{fig_ind},:,:);
%     line2EventMat = eventMat(figExcInd{fig_ind},:,:);
%     
%     % If there are two lines to plot, cycle through them.
%     if ~isempty(line1Array) && ~isempty(line2Array)
%       for group_ind = groupIterInd
%         % Generate Data
%         line1Data = vertcat(line1Array{:,group_ind});
%         line2Data = vertcat(line2Array{:,group_ind});
%         lineMean = [mean(line1Data, 1); mean(line2Data, 1)];
%         lineErr = [std(line1Data)/sqrt(size(line1Data,1)); std(line2Data)/sqrt(size(line2Data,1))];
%         
%         if params.fixAlign
%           fixMean = mean(mean(lineMean(:,500:800),1));
%           for line_i = 1:size(lineMean,1)
%             lineMean(line_i,:) = lineMean(line_i,:) - mean(lineMean(line_i,500:800)) + fixMean;
%           end
%         end
%         
%         % Prepare figure
%         plotTitle = sprintf('%s - %s, %s', ['mean PSTH of ' figTitInd{fig_ind}], groupingType{group_ind}, normTag);
%         h = figure('NumberTitle', 'off', 'Name', plotTitle,'units','normalized','outerposition',[0 0 params.plotSizeLineBroadCatPlot]);
%         params.lineProps = lineProps;
%         
%         if params.addSubEventBars
%           figLegItr = [figLegends{fig_ind}, strrep(eventList, '_', ' ')];
%         else
%           figLegItr = figLegends{fig_ind};
%         end       
%         
%         [ ~,  ~, ~, legendH] = plotPSTH(lineMean, lineErr, [], params, 'line', plotTitle, figLegItr);
%         axesH = findobj(gcf, 'Type', 'Axes');
%         axesH.FontSize = 15;
%         hold on
%         legendH.Location = 'northeast';
%         ylabel('Normalized Activity');
%         
%         % Add event plots underneath the traces using eventMat and
%         % add_bar_to_plots.m
%         if params.addSubEventBars
%           % Prepare data matrices to plot
%           line1DataMat = squeeze(sum(line1EventMat, 1))';
%           line2DataMat = squeeze(sum(line2EventMat, 1))';
%           eventCount = size(line1DataMat, 1);
%           comboMat = [line1DataMat; line2DataMat];
%           
%           % Shift to fit and pad pre and post.
%           comboMat = comboMat(:, 1:params.psthImDur);
%           comboMat = [zeros(size(comboMat,1), params.psthPre), comboMat, zeros(size(comboMat,1), params.psthPost)];
%           cMatSize = size(comboMat);
%           comboMatShow = ~comboMat == 0;
%           
%           % Package for the add_bars_to_plots
%           comboMatArray = mat2cell(comboMat, repmat(eventCount, [2,1]), cMatSize(2));
%           comboMatShowArray = mat2cell(comboMatShow, repmat(eventCount, [2,1]), cMatSize(2));
%           % Prepare labels/colors
%           colorMat = lineProps.col;
%           
%           % Add the plots
%           [newBarImg, barDummyHands] = add_bars_to_plots([], [], comboMatArray, colorMat, comboMatShowArray, []);
%           plotTitle = horzcat(plotTitle, '_eventBar');
%         end
%         
%         % Save the figure
%         saveFigure(params.outputDir, ['5. ' plotTitle], [], figStruct, []);
%         
%         % If desired, Generate plots of constituient traces.
%         if params.splitContrib
%           lineData = {line1Array(:, group_ind), line2Array(:, group_ind)};
%           lineLabels = {line1Labels, line2Labels};
%           for line_i = 1:length(lineData)
%             j = figure();
%             figTit = sprintf('%s - %s', figLegends{fig_ind}{line_i}, groupingType{group_ind});
%             meanLines = cellfun(@(x) mean(x), lineData{line_i}, 'UniformOutput', 0);
%             stdLines = cellfun(@(x) std(x)/sqrt(size(x,1)), lineData{line_i}, 'UniformOutput', 0);
%             meanLines = vertcat(meanLines{:});
%             stdLines = vertcat(stdLines{:});
%             
%             % Sort things based on magnitutde of the peak
%             peaks = max(meanLines,[],2);
%             [~, sortOrder] = sort(peaks, 'descend');
%             
%             plotPSTH(meanLines(sortOrder,:), stdLines(sortOrder,:), [], params, 'line', figTit, lineLabels{line_i}(sortOrder,:));
%             % Save the figure
%             saveFigure(params.outputDir, ['5.1. ' figTit], [], figStruct, [])            
%           end
%         end
%         
%       end
%     end
%   end
%   
% end

end

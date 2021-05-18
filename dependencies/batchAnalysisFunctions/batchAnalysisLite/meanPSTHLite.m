function meanPSTHLite(spikePathBank, params, figStruct)
% Function which combines stimulus presentations across all runs in the spikeDataBank.
% Inputs include spikeDataBank and list of parameters.
disp('Starting mean PSTH Analysis...');

plotTopStim = params.meanPSTHParams.plotTopStim;
removeFollowing = params.meanPSTHParams.removeFollowing;
fixAlign = params.meanPSTHParams.fixAlign;
normalizeTraces = params.meanPSTHParams.normalize;
stimInclude = params.meanPSTHParams.stimInclude;
removeRewardEpoch = params.meanPSTHParams.removeRewardEpoch;
firstXRuns = params.meanPSTHParams.firstXRuns;

plotIndParams = params.NDTParams.spikeToRasterParams.plotIndParams; % Forgive me, Father.
paradigmList = unique(spikePathBank.paradigmName);

% For later reference - create a list of unique channels per run.
spikePathBank = channelPerRunArray(spikePathBank);

for par_i = 1:length(paradigmList)
  
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  spikePathBankParadigm = spikePathBank(pInd, :);
  
  % Step 1 - Unpack variables, generate those needed for plotting.
  allStimuliVec = unique(vertcat(spikePathBankParadigm.stimuli{:}));
  groupingType = {'Unsorted', 'Unit', 'MUA'};
    
  % Generate grid for indexing into individual runs and extracting relevant
  % PSTHes.
  runList = spikePathBankParadigm.Properties.RowNames;
  small2BigInd = zeros(length(allStimuliVec), length(runList));
  for run_ind = 1:length(runList)
    [~, small2BigInd(:,run_ind)] = ismember(allStimuliVec, spikePathBankParadigm.stimuli{run_ind});
  end
  
  stimCatPresCounts = sum(logical(small2BigInd),2);
    
  % Step 2 - modify PSTHes - rewardEpoch removal, firstXRuns, animations,
  % topStimOnly
  if plotTopStim
    %Initialize a keep index
    keepInd = true(size(allStimuliVec));
    
    if params.meanPSTHParams.plotTopStim
      keepInd = keepInd & (stimCatPresCounts >= params.meanPSTHParams.topStimPresThreshold);
    end
    
    % Once keep ind has been generated across all switches, edit appropriate
    % structures
    stimCatPresCounts = stimCatPresCounts(keepInd);
    allStimuliVec = allStimuliVec(keepInd);
        
  end
  
  % Generate the figure directory
  dirTags = {'NoReward', 'fixAligned', 'Normalized', 'Pres Threshold', 'AnimOnly', 'NoAnim', sprintf('first%d', params.meanPSTHParams.firstXRuns)};
  dirTagSwitch = [removeRewardEpoch fixAlign normalizeTraces plotTopStim stimInclude == 1 stimInclude == 2 firstXRuns];
  outputDir = [params.meanPSTHParams.outputDir, strjoin(dirTags(logical(dirTagSwitch)),' - ')];
  
  if ~exist(outputDir, 'dir')
    mkdir(outputDir)
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
  
%   if params.meanPSTHParams.removeRewardEpoch
%     variables2Extract(1:2) = strcat(variables2Extract(1:2), 'NoReward');
%   end
    
  % [psthByImage, psthErrByImage, psthByCategory, psthErrByCategory, psthParamsPerRun] = spikePathLoad(spikePathBank, variables2Extract, params.spikePathLoadParams);
  [psthByImage, psthByCategory, psthParamsPerRun] = spikePathLoad(spikePathBankParadigm, variables2Extract, params.spikePathLoadParams);
  
  selTable = spikePathBankParadigm.selTable;
  paradigmInd = unique(spikePathBankParadigm.paradigmInd);
  runList = extractAfter(spikePathBankParadigm.Properties.RowNames, 'S');
  
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
    
    % Cycle through the run data, collecting the traces of interest
    paradigmName = spikePathBankParadigm.paradigmName{1};
    psthParamTimes = psthParamsPerRun{1};
    if params.meanPSTHParams.removeRewardEpoch
      psthParamTimes.psthPost = 0;
    end
    
    % check to see if analysisGroups need ot be plotted as well
    if params.meanPSTHParams.analysisGroupPSTH
      analysisGroups2Plot = [paradigmInd; fields(params.meanPSTHParams.analysisGroups.(paradigmName))];
    else
      analysisGroups2Plot = {paradigmInd};
    end
    
    for plot_i = 1:length(analysisGroups2Plot)
      for data_i = 1:length(psth2Collect)
        
        traces2Extract = eval(psth2Collect{data_i}); % way to extract Image or Category data.
        if isnumeric(analysisGroups2Plot{plot_i})
          [meanTraces, errTraces, counts, paradigmStimPlot] = extractMeanTraces(traces2Extract, spikePathBankParadigm.(spikePathBankField{data_i}), [], 0, 1);
        else
          [meanTraces, errTraces, counts, paradigmStimPlot] = extractMeanTraces(traces2Extract, spikePathBankParadigm.(spikePathBankField{data_i}), params.meanPSTHParams.analysisGroups.(paradigmName).(analysisGroups2Plot{plot_i}), 0, 1);
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
          [~, cbh] = plotPSTH(plotData, plotErr, axes(), psthParamTimes, 'color', catPSTHTitle, allStimuliLabel);
          title(catPSTHTitle);
          xlabel('Time from Stimulus Onset (ms)');
          
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

          saveFigure(outputDir, sprintf('%d. %s', plotNum(data_i), catPSTHTitle), [], figStruct, [])
          clear allStimuliLabel
        end
      end
    end
  end
  
  % Plot 2 + 3 - analysisGroup Plot - 'All head Turning vs non-Head turning'
  if params.meanPSTHParams.analysisGroupPSTH
    % Extract the groups to plot per paradigm.
    
    for ii = 1:length(paradigmInd)
      % Cycle through the run data, collecting the traces of interest
      
      paradigmInd = spikePathBankParadigm.paradigmInd == paradigmInd(ii);
      paradigmName = spikePathBankParadigm.paradigmName{find(paradigmInd,1)};
      psthParamTimes = psthParamsPerRun{find(paradigmInd,1)};
      spikePathBankParadigm = spikePathBankParadigm(paradigmInd, :);
      
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
          saveFigure(outputDir, sprintf('%d. %s', data_i, catPSTHTitle), [], figStruct, [])
          clear allStimuliLabel
        end
      end
      
    end
    
  end
  
  % Plot 4 - Create a plot with means for selective units between conditions, and non-selective units.
  if params.meanPSTHParams.selPSTH
    UnitTypes = params.selParam.UnitTypes;
    UnitTypePlot = params.selParam.UnitTypePlot;
    stimParamData = load(params.stimParamsFilename);
          
      % Cycle through the run data, collecting the traces of interest
      psthParamTimes = psthParamsPerRun{1};

      psthByImageParadigm = psthByImage;
      runListParadigm = runList;
      plotIndParams.plotLabels = {{'agents', 'socialInteraction'}};
      
      % Combine
      paradigmSelTable = spikePathBankParadigm.selTable;
      paradigmSelTable = vertcat(paradigmSelTable{:});
      tableVars = paradigmSelTable.Properties.VariableNames';
      
      % Identify variables with a Cohen's D - these also have selectivity
      % indices.
      variableList = tableVars(contains(tableVars, 'subSel') & contains(tableVars, 'cohensD') & ~contains(tableVars, 'reward'));
      eventVarList = strcat(extractBefore(variableList, '_cohensD'), '_selInd');
      
      % Create an array of selectivity type, unitType, Sign, and dataType
      dataType = {'agents', 'socialInteraction'};
      normalizedActivity = cell(length(variableList), length(UnitTypes), 2, length(dataType));
      [normalizedActivity{:}] = deal(nan(1e5 ,size(psthByImageParadigm{1}{1}{1},2)));
      
      % For every event, find the unit traces selective for it.
      for event_i = 1:length(variableList)
        for sign_i = 1:2
          
          % Seperate Increases and decresed in firing
          if sign_i == 1
            eventSelective = paradigmSelTable.(variableList{event_i}) > 0 & paradigmSelTable.(eventVarList{event_i});
          else
            eventSelective = paradigmSelTable.(variableList{event_i}) < 0 & paradigmSelTable.(eventVarList{event_i});
          end
          paradigmSelTableEvent = paradigmSelTable(eventSelective,:);

          % Cycle for normal units and MUA.
          for unitType_i = 1:length(UnitTypes)
            unitTypeIndex = contains(paradigmSelTableEvent.unitType, UnitTypes{unitType_i});
            paradigmSelTableEventUnit = paradigmSelTableEvent(unitTypeIndex, :);
            runListUnit = strcat(paradigmSelTableEventUnit.dateSubj, paradigmSelTableEventUnit.runNum);
            runListChannels = spikePathBank{strcat('S', runListUnit),"channelNumbers"};
            
            chanUnit = double(extractAfter(paradigmSelTableEventUnit.channel, 'Ch'));
            unitUnit = paradigmSelTableEventUnit.unitType;
            
            % for each unit, Grab the Activity traces, and store
            for unit_i = 1:length(chanUnit)
              % Identify the trace of activity to collect
              runInd = strcmp(runListParadigm, runListUnit{unit_i});
              unitDataInd = runListChannels{unit_i} == chanUnit(unit_i); % runs without all the channels need to have the channel number converted to an index.
              
              % Extract the trace
              if strcmp(unitUnit(unit_i), 'MUA')
                psthRun = psthByImageParadigm{runInd}{unitDataInd}{end};
              else
                unit_Index = double(extractAfter(paradigmSelTableEventUnit.unitType(unit_i), 'U')) + 1;
                psthRun = psthByImageParadigm{runInd}{unitDataInd}{unit_Index};
              end
              
              % Grab the traces, and the traces they belong to.
              runStim = spikePathBankParadigm.stimuli{runInd};
              traceIndexPerStimuli = plotIndex(runStim, plotIndParams);
              
              % store traces
              for trace_i = 1:2
                nanInd = find(isnan(normalizedActivity{event_i, unitType_i, sign_i, trace_i}(:,1)), 1);
                traces2Store = traceIndexPerStimuli == trace_i;
                normalizedActivity{event_i, unitType_i, sign_i, trace_i}(nanInd:nanInd+sum(traces2Store)-1, :) = psthRun(traces2Store,:);
              end
            end
            
          end
        end
      end
      
      % After retrieving all the traces, remove the NaNs from the stored
      % values
      normalizedActivity = cellfun(@(x) x(1:find(isnan(x(:,1)))-1, :), normalizedActivity, 'UniformOutput', false);
              
      % Plot Time - % normalizedActivity - selectivity type, unitType, Sign, and dataType
      signPlot = {'Increased Activity', 'Decreased Activity'};
      variableListPlot = extractBetween(variableList, 'subSel_', '_c');
      variableListPlot = strrep(variableListPlot, '_', ' ');
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
            psthTitle = sprintf('Norm Pop %s - %s selective (%s)', UnitTypePlot{unitType_i}, variableListPlot{event_i}, signPlot{sign_i});
            h = figure('Name', psthTitle, 'units', 'normalized', 'outerposition', [0 0 params.meanPSTHParams.plotSizeLineCatPlot]);
            plotPSTH(meanTraces, stdTraces, [], psthParamTimes, 'line', psthTitle, legendPlots)
            ylabel('Change in Activity relative to Baseline, Normalized')
            h.Children(2).FontSize = 16;
            saveFigure(outputDir, ['4. ' psthTitle], [], figStruct, [])
            
          end
        end
      end
      
  end
  
end

end



function holder()

eventMat = generateEventImage(eventData, params.psthImDur);

line1EventMat = eventMat(figIncInd{fig_ind},:,:);
line2EventMat = eventMat(figExcInd{fig_ind},:,:);


if params.addSubEventBars
  % Prepare data matrices to plot
  line1DataMat = squeeze(sum(line1EventMat, 1))';
  line2DataMat = squeeze(sum(line2EventMat, 1))';
  eventCount = size(line1DataMat, 1);
  comboMat = [line1DataMat; line2DataMat];
  
  % Shift to fit and pad pre and post.
  comboMat = comboMat(:, 1:params.psthImDur);
  comboMat = [zeros(size(comboMat,1), params.psthPre), comboMat, zeros(size(comboMat,1), params.psthPost)];
  cMatSize = size(comboMat);
  comboMatShow = ~comboMat == 0;
  
  % Package for the add_bars_to_plots
  comboMatArray = mat2cell(comboMat, repmat(eventCount, [2,1]), cMatSize(2));
  comboMatShowArray = mat2cell(comboMatShow, repmat(eventCount, [2,1]), cMatSize(2));
  % Prepare labels/colors
  colorMat = lineProps.col;
  
  % Add the plots
  [newBarImg, barDummyHands] = add_bars_to_plots([], [], comboMatArray, colorMat, comboMatShowArray, []);
  plotTitle = horzcat(plotTitle, '_eventBar');
end
end
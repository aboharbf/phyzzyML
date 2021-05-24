function [subEventPSTHStruct] = subEventPSTH(spikePathBank, params, figStruct)
% Function compiles all the PSTHes of subEvents and plots them.
% One version uses spikeDataBank and pregenerated 'subEventSig' to pool and
% plot. 2nd segment uses stimPSTH to take slices out of stimuli PSTHes to
% generat plots.
disp('running subEventPSTH()...');
baselineSubtract = 1;
runList = spikePathBank.Properties.RowNames;
removeExcessNaNs = @(cellWNans) cellWNans(~isnan(cellWNans(:,1)), :);
outputDir = params.subEventPSTHParams.outputDir;

% Load and extract key things from subEvents. Plots 1 and 2 only need the
% eventList, 3 and after use more features.
load(params.subEventPSTHParams.eventData);

% Extract data from meanPSTHStruct
paradigmList = unique(spikePathBank.paradigmName);

for par_i = 1:length(paradigmList)
  
  paradigmIndex = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  spikePathBankParadigm = spikePathBank(paradigmIndex,:);
  
  [psthByImagePerRun, eventDataArray, dataTypeArray, groupTypeArray, eventList, subEventSig, psthParams2Use] = generateEventDataArray(spikePathBankParadigm, params);
    
  % groupType = {'Unsorted', 'Unit', 'MUA'};
  
  % Temp to account for lack of parameters in data
  psthParam.psthPre = abs(subEventSig.psthWindow(1));
  psthParam.psthImDur = subEventSig.psthWindow(2);
  psthParam.psthPost = size(eventDataArray{1,1,1},2)-subEventSig.psthWindow(2)-abs(subEventSig.psthWindow(1))-1;
  
  % Plotting variables/code
  if params.subEventPSTHParams.normalize
    normTag = ' - normalized';
  else
    normTag = '';
  end
  
  % Plot 1 - subEvent PSTHs, sorted by X
  % Make a PSTH of each event across all its repetitions.
  % {'Run of Day', 'Grid Hole', 'Recording Depth', 'Run Ind'};
  % Must be the same order as the final parts of stimPSTH
  
  sortType = {'Run Index'};
  [sortMat, sortLabel] = deal(cell(size(eventDataArray,1), size(eventDataArray,2), length(sortType)));
  
  % Groups to creating sorting indicies for.
  groupIterInd = 3; 
  
  % Generate the PSTH sorting indicies
  for event_i = 1:size(sortMat,1)
    for group_i = groupIterInd
      for sort_i = 1:length(sortType) % the sortings that need processing
        
        switch sort_i
          case 1
            sortMat{event_i, group_i, sort_i} = [1:length(eventDataArray{event_i, group_i, strcmp(dataTypeArray, sortType{sort_i})})]';
            tmpLabels = eventDataArray{event_i, group_i,strcmp(dataTypeArray, sortType{sort_i})};
            labeledIndex = 1:5:length(tmpLabels);
        end
        
        % Use the labeledIndex to generate the proper label array with 0's
        % everywhere else.
        if params.subEventPSTHParams.sparseLabels
          finalLabels = cell(length(tmpLabels),1);
          for tmp_i = 1:length(labeledIndex)
            finalLabels{labeledIndex(tmp_i)} = tmpLabels{labeledIndex(tmp_i)};
          end
        else
          finalLabels = tmpLabels;
        end
        
        % Store vector into larger array for later rference
        sortLabel{event_i, group_i, sort_i} = finalLabels;
        
      end
    end
  end
  
  % Plots 1 and 2 treat all instances of a particular subEvent the same.
  % Plots 3 and 4 Do not.
  
  % Plot 1 - Individual PSTHes for every instance of a subEvent across Runs, Sorted
    eventListPlot = strrep(eventList, '_', ' ');
  if params.subEventPSTHParams.allRunStimPSTH
    
    % Iterate through PSTH, generating plots
    for event_i = 1:size(eventDataArray,1)
      for sort_i = 1
        for group_i = groupIterInd
          eventNamePlot = strrep(eventList{event_i}, '_', ' ');
          figTitle = sprintf('%s - %s PSTHs, Sorted by %s, %s', eventNamePlot, groupTypeArray{group_i} , normTag);
          h = figure('NumberTitle', 'off', 'Name', figTitle, 'units', 'normalized', 'outerposition', [0 0 params.subEventPSTHParams.plotSizeAllRunStimPSTH]);
          sgtitle(figTitle)
          
          sortIndex = sortMat{event_i, group_i, sort_i};
          
          % Break up PSTHes with many lines into subplots.
          TracesPerPlot = ceil(length(sortIndex)/3);
          if TracesPerPlot < 20
            subplot2Plot = 1;
          elseif TracesPerPlot < 45
            subplot2Plot = 2;
          else
            subplot2Plot = 3;
          end
          TracesPerPlot = ceil(length(sortIndex)/subplot2Plot);
          
          plotStarts = 1:TracesPerPlot:length(sortIndex);
          plotEnds = [plotStarts(2:end)-1, length(sortIndex)];
          
          subplotAxes = gobjects(subplot2Plot,1);
          cbHandle = gobjects(subplot2Plot,1);
          for plot_i = 1:subplot2Plot
            plotLabels = sortLabel{event_i, group_i, sort_i}(plotStarts(plot_i):plotEnds(plot_i));
            plotData = eventDataArray{event_i,group_i,1}(plotStarts(plot_i):plotEnds(plot_i), :);
            %sortIndex = sortMat{event_i, group_i, sort_i}(plotStarts(plot_i):plotEnds(plot_i));
            psthAxes = subplot(1,subplot2Plot,plot_i);
            [subplotAxes(plot_i), cbHandle(plot_i)] = plotPSTH(plotData, [], psthAxes, psthParam, 'color', [], plotLabels); %(sortIndex,:)
            cbHandle(plot_i).Label.FontSize = 12;
            if plot_i == subplot2Plot
              if params.subEventPSTHParams.normalize == 1
                cbHandle(plot_i).Label.String = 'Signal Change relative to Baseline (%)';
              elseif params.subEventPSTHParams.normalize == 2
                cbHandle(plot_i).Label.String = 'Z score relative to fixation';
              end
            else
              delete(cbHandle(plot_i))
            end
            set(gca,'FontSize',10,'TickLength',[.01 .01],'LineWidth',.25);
          end
          linkprop(subplotAxes, 'CLim');
          
          saveFigure(outputDir, ['1. ' figTitle], [], figStruct, []);
        end
      end
    end
    
  end
  
  % Plot 2 - mean subEvent PSTH vs Null, Line plot
  if params.subEventPSTHParams.meanSubEventPSTH
    for event_i = 1:length(eventListPlot)
      
      % Check if there are any actual traces
      emptyInd = cellfun('isempty', squeeze(eventDataArray(event_i,:,:)));
      if all(emptyInd(:))
        continue
      end
      
      % Prepare Plot Titles
      figTitle = sprintf('mean subEvent PSTH for %s %s', eventListPlot{event_i}, normTag);
      
      if baselineSubtract
        figTitle = strcat(figTitle, ' baselineSub');
      end
      
      h = figure('NumberTitle', 'off', 'Name', figTitle, 'units','normalized','outerposition',[0 0 0.5 1]);
      sgtitle(figTitle)
      
      for group_i = 1:size(eventDataArray,2)
        subplot(length(groupTypeArray), 1, group_i);
        [plotData, plotErr] = raw2meanSEM(eventDataArray(event_i, group_i, 1:2), 0);
        
        if baselineSubtract
          baselinePeriod = (psthParam.psthPre/2):psthParam.psthPre;
          plotData = plotData - mean(plotData(:,baselinePeriod), 2);
        end
        
        plotPSTH(plotData, plotErr, [], psthParam, 'line', groupTypeArray{group_i}, {'Event', 'Null'});
        
      end
      
      saveFigure(outputDir, ['2. ' figTitle], [], figStruct, []);
      
    end
  end
  
  % Variables for Plot 3 and 4 - These plots rely on activity coming from the
  % the entire stimulus PSTH, instead of individual run data loaded from sigStruct
  % pre-cut. The code below is used to generate all the needed variables to
  % extract and plot that data.
  
  % Extract Event related data
  allStimuliVec = unique(vertcat(spikePathBankParadigm.stimuli{:}));
  allStimuliNames = strrep(allStimuliVec, '_', ' ');
  allStimuliNames = extractBefore(allStimuliNames, '.avi');
  eventList = eventData.Properties.VariableNames;                     % Create structure w/o Blinks/Saccades;
  eventListPlot = strrep(eventList, '_', ' ');
    
  % Remove non-represented variables
  eventDataPar = eventData(allStimuliVec, :);
  eventLists = eventDataPar.Properties.VariableNames;
  stimEventCountarray = cellfun(@(x) size(x, 1), table2cell(eventDataPar));
  stimEventLogicArray = stimEventCountarray ~= 0;
  stimEvent2Plot = any(stimEventLogicArray,2);        % Only worry about plotting stimuli below which have events.
  stimEvent2PlotCount = sum(stimEventCountarray,2);   % Only worry about plotting stimuli below which have events.

  % Decide which variables to extract from the runs (Any amount of
  % Proc/Normalizing should be considered).
  %   variables2Extract = {'psthByImage', 'psthByCategory', 'psthParams'};
  
  psthPreData = psthParams2Use.psthPre;   % Used to find the correct spot in larger vector.
  psthPre = params.subEventPSTHParams.psthPre;
  psthImDur = params.subEventPSTHParams.psthImDur;
  psthPost = params.subEventPSTHParams.psthPost;
  traceLength = size(psthByImagePerRun{1}{1}{1}, 2);
  
  meanPSTHStruct.addLegend = false;
  meanPSTHStruct.psthPre = psthParams2Use.psthPre;
  meanPSTHStruct.psthImDur = psthParams2Use.psthImDur;
  meanPSTHStruct.psthPost = psthParams2Use.psthPost;
  
  % Compile a stimPSTH Structure.
  [stimPSTH, dataType] = compilePsthPerStimArray(spikePathBankParadigm.stimuli, psthByImagePerRun, allStimuliVec);
   
  % Process that further into an array of subEvent data.
  [subEventData, dataType2] = compileSubEventPSTHArray(stimPSTH, dataType, eventDataPar, eventLists, groupIterInd, psthPreData, params);
  
  % Plot 3 - all subEvent PSTH from stimuli
  if params.subEventPSTHParams.stimPlusEvents_extracted && strcmp(paradigmList{par_i}, 'NaturalSocial')
    for stim_i = 1:size(stimPSTH,1)
      if stimEvent2Plot(stim_i) % If a stimuli has events...
        % Get the relevant eventDataPar
        eventDataParStim = eventDataPar(allStimuliVec{stim_i},:);
        
        % Use the number of event in a stimulus to generate plot indices.
        numPlotsRight = stimEvent2PlotCount(stim_i);
        subplotInds = 2:2:(numPlotsRight*2);
        
        for group_i = groupIterInd
          subplotCount = 0; % Reset this per figure
          
          fullStimTraces = stimPSTH{stim_i, group_i, strcmp(dataType, 'PSTH')};
          
          % If traces exceed some number, simply retrieve a slice.
          traces2plot = 30;
          if size(fullStimTraces, 1) > traces2plot
            fullStimTracesPlot = fullStimTraces(randi(size(fullStimTraces, 1), [traces2plot, 1]), :);
          else
            fullStimTracesPlot = fullStimTraces;
          end
          
          % If Desired, remove reward activity.
          removeReward = 0;
          rewRemThres = 4;
          if removeReward
            rewardPeaks = any(fullStimTraces(:, meanPSTHStruct.psthPre + meanPSTHStruct.psthImDur:end) > rewRemThres, 2);
            [tmpMeans, tmpSEM] = raw2meanSEM([{fullStimTraces}; {fullStimTraces(~rewardPeaks, :)}], 0);
          else
            rewardPeaks = zeros(size(fullStimTraces,1),1);
          end
          
          % Plot, Left side
          sgFigTitle = sprintf('Events in %s - %s', allStimuliNames{stim_i}, groupTypeArray{group_i});
          figure('NumberTitle', 'off', 'Name', sgFigTitle, 'units','normalized','outerposition',[0 0 params.subEventPSTHParams.plotSizeAllRunStimPSTH_extracted])
          sgtitle(sgFigTitle)
          
          % if there are no traces removed due to a reward filter, first
          % column will be 1 plot.
          rowCountLeft = 1 + any(rewardPeaks);
          
          % Normal activity
          subplot(rowCountLeft,2,1);
          plotPSTH(fullStimTracesPlot, [], [], meanPSTHStruct, 'line', 'All Activity Traces', [])
          
          % Reward Only Activity
          if removeReward && any(rewardPeaks)
            subplot(rowCountLeft,2,3)
            plotPSTH(tmpMeans, tmpSEM, [], meanPSTHStruct, 'line', sprintf('Activity Traces w.o reward activity > %d', rewRemThres), [])
            legend({'All Traces - Reward Peak Traces', 'All Traces'}, 'location', 'northwest')
          end
          
          for event_i = 1:length(eventList)
            % Grab event data for each event, adjusting it by the psthPre.
            evTable = eventDataParStim.(eventList{event_i}){1};
            eventStarts = round(evTable.startTime);
            eventEnds = round(evTable.endTime);
            eventLengths = eventEnds - eventStarts;
            
            for inst_i = 1:length(eventStarts)
              % Expand the desired stop and start by the amount outside of
              % the event you'd like to see.
              startT = eventStarts(inst_i) + psthPreData - psthPre;
              endT = min(eventEnds(inst_i) + psthPreData + psthPost, traceLength);            % If endT exceeds the length of the trace, then adjust it.
              endDist = min(traceLength - eventEnds(inst_i), psthPost);
              
              eventPSTHData = fullStimTraces(~rewardPeaks, startT:endT);
              [traceMean, traceSEM] = raw2meanSEM({eventPSTHData},0);
              
              % Plot
              vertLinePos = 200;
              figTitle = sprintf('%s (Red Line @ %d ms)', eventListPlot{event_i}, vertLinePos);
              subplotCount = subplotCount + 1;
              psthAxes = subplot(numPlotsRight, 2, subplotInds(subplotCount));
              
              subplotParams.psthPre = psthPre;
              subplotParams.psthImDur = eventLengths(inst_i);
              subplotParams.psthPost = endDist;
              subplotParams.addLegend = false;
              
              plotPSTH(traceMean, traceSEM, psthAxes, subplotParams, 'line', [], []);
              
              % Modify X axis
              psthAxes.XTick = sort([0 eventLengths(inst_i) vertLinePos]);
              psthAxes.XTickLabel = sort([eventStarts(inst_i) eventEnds(inst_i) eventStarts(inst_i) + vertLinePos]);
              
              % Add Red line @ 200 after stim onset to give sense of length
              line([vertLinePos vertLinePos], ylim(),'LineWidth', 5, 'color', 'r');
              
              if ~(subplotCount == numPlotsRight)
                psthAxes.XLabel.String = '';
              end
              
              title(figTitle);
            end
          end
          
          if params.subEventPSTHParams.stimPlusEvents_extracted
            saveFigure(outputDir, ['3. ' sgFigTitle], [], figStruct, []);
          end
          
        end
        
      end
    end
  end
  
  % Plot 3.1 - Event PSTH Color Plots, based on 'eventStorage' from plot 3.
  if params.subEventPSTHParams.eventPsthColorPlots
    for event_i = 1:size(subEventData,1)
      for inst_i = 1:length(subEventData{event_i, 1})
        % Extract traces for an instance of an event.
        instTraces = subEventData{event_i, strcmp(dataType2, 'subEventPSTH')}{inst_i};
        instRunLab = runList(subEventData{event_i, strcmp(dataType2, 'RunInd')}{inst_i});
        
        % Traces are padded with NaNs to make them the same length. Remove
        % this padding.
        NaNCutOffInd = find(isnan(instTraces(1,:)), 1);
        if ~isempty(NaNCutOffInd)
          psthArray = instTraces(:,1:NaNCutOffInd);
        else
          psthArray = instTraces;
        end
        
        figTitle = sprintf('PSTH for Instances of %s in %s', eventListPlot{event_i}, subEventData{event_i, 2}{inst_i});
        figure('NumberTitle', 'off', 'Name', figTitle, 'units','normalized','outerposition',[0 0 params.subEventPSTHParams.plotSizeAllRunStimPSTH_extracted])
        
        % Add the grand mean
        psthArray = [psthArray; nanmean(psthArray)];
        psthLabels = [instRunLab; 'Mean'];
        
        % Find and subtract the mean of each trace
        for trace_i = 1:size(psthArray, 1)
          psthArray(trace_i, :) = psthArray(trace_i, :) - mean(psthArray(trace_i, 1:psthPre));
        end
        
        % Plot using plotPSTH, Line version
        tmpParams.psthPre = psthPre;
        tmpParams.psthImDur = size(psthArray,2) - psthPre - 1;
        tmpParams.psthPost = 0;
        tmpParams.addLegend = 1;
        
        psthHand = plotPSTH(psthArray, [], [], tmpParams, 'color', figTitle, psthLabels);
        psthHand.XLabel.String = 'Time from Event Onset (ms)';
        
        saveFigure(outputDir, sprintf('3.1. %s - %', paradigmList{par_i}, figTitle), [], figStruct, []);
        
      end
    end
  end
  
  % Plot 3.2 - Event PSTH Mean Line Plots
  if params.subEventPSTHParams.eventPsthMeanLinePlots
    for event_i = 1:size(eventListPlot,1)
      for group_i = 3
        % Extract relevant info for this event.
        allTraces = subEventData{event_i, 1}';
        allTraceLabels = [subEventData{event_i, 2}, {''}]';
        
        % Use function to calculate means and SEM. include grandMean.
        grandMeanSwitch = 1;
        [meanStack, ~, nStack] = raw2meanSEM(allTraces, grandMeanSwitch);
        
        % do quick loop to find the appropriate offsets.
        traceCountPlusMean = length(allTraces)+1;
        traceBaseline = zeros(traceCountPlusMean, 1);
        for trace_i = 1:traceCountPlusMean
          if trace_i ~= traceCountPlusMean
            traceBaseline(trace_i) = mean(allTraces{trace_i}(:, 1:psthPre), 'all');
          else
            grandStack = vertcat(allTraces{:});
            traceBaseline(trace_i) = mean(grandStack(:,1:psthPre), 'all');
          end
        end
        
        % Apply the offset to the means
        for mean_i = 1:size(meanStack, 1)
          meanStack(mean_i) = meanStack(mean_i) - traceBaseline(mean_i);
        end
        
        %Complete Trace labels.
        allTraceLabels{end} = sprintf('All %s mean', eventListPlot{event_i});
        allTraceLabels = strcat(allTraceLabels, ' (', strtrim(string(num2str(nStack))), ')');
        
        % Plot using plotPSTH, Line version
        tmpParams.psthPre = psthPre;
        tmpParams.psthImDur = size(meanStack,2) - psthPre - 1;
        tmpParams.psthPost = 0;
        tmpParams.addLegend = 1;
        
        % Generate figure and titles.
        figTitle = sprintf('Event means for all instances of %s in paradigm %s', eventListPlot{event_i}, paradigmList{par_i});
        figure('NumberTitle', 'off', 'Name', figTitle, 'units','normalized','outerposition',[0 0 params.subEventPSTHParams.plotSizeAllRunStimPSTH_extracted])
        hold on
        
        plotPSTH(meanStack, [], [], tmpParams, 'line', figTitle, allTraceLabels)
        
      end
    end
  end
  
  % Plot 4 - mean subEvent PSTH, based on extracted slices
  % Generates Plots for each event by overlaying particular event types and
  % generating null traces by grabbing the same slice of time in other
  % events.
  
  groupIterInd = {1, 2, 3, [1, 2]};
  
  if params.subEventPSTHParams.meanSubEventPSTH_extracted
    % This code generate an event aligned mean activity
    
    % Generate a list of events to plot - both simple and comboEvents
    simpleEvents = eventListPlot;
    comboEvents = {'headTurn All', 'all Turns'};
    
    % For the sake of data gathering, denote which are combo events.
    comboInd = [false(1, length(eventListPlot)), true(1, length(comboEvents))];
    eventListPlotCombo = [eventListPlot,  comboEvents];
    
    % Generate indicies for later data gathering. Simple events are 1
    % number, Combo events are more.
    comboGroups = mat2cell(1:length(simpleEvents), 1, ones(length(simpleEvents),1)');
    comboGroups{strcmp(eventListPlotCombo, 'headTurn All')} = [find(strcmp(eventListPlotCombo, 'headTurn right')),  find(strcmp(eventListPlotCombo, 'headTurn left'))];
    comboGroups{strcmp(eventListPlotCombo, 'all Turns')} = [comboGroups{strcmp(eventListPlotCombo, 'headTurn All')}, find(strcmp(eventListPlotCombo, 'bodyTurn'))];
    
    % For combination type events, like 'all Turns', I need to compile things before.
    % 3rd D - 1 = real data, 2 = NullData.
    dataType3 = {'PSTH', 'PSTHnull'};
    eventPSTHData = cell(length(eventListPlotCombo), 3, length(dataType3));
    
    % Initialize structures for holding both traces and their null
    % equivalents.
    for event_i = 1:length(eventListPlotCombo)
      % If the event exists in this paradigm and is not a combo event...
      if ~comboInd(event_i)
        
        % Initialize the entries for this event using the appropriate
        % template.
        [eventPSTHData{event_i, :, 1}] = deal(nan(1e5, psthPre + psthImDur + psthPost + 1));
        [eventPSTHData{event_i, :, 2}] = deal(nan(4e5, psthPre + psthImDur + psthPost + 1));
        
        for group_i = 1:3
          % Simple events, defined by eventDataPar.
          stimLogicArrayEvent = stimEventLogicArray(:,event_i);
          stim2Check = find(stimLogicArrayEvent);
          
          for stim_i = stim2Check'
            % Go through stimuli with events
            eventSpecTable = eventDataPar{stim_i, event_i}{1};
            eventStarts = round(eventSpecTable.startTime);
            eventEnds = round(eventSpecTable.endTime);
            
            for inst_i = 1:length(eventStarts)
              % For every instance of an event, identify start and stop times
              startT = (eventStarts(inst_i) + psthPreData) - psthPre;               % Event start - desired preEvent time.
              endT = (eventStarts(inst_i) + psthPreData) + psthImDur + psthPost;    % Event duration + desired postEvent time.
              
              % MAX LENGTH TEMPLATE NEEDS EDITING - USE PSTH Params rather
              % than event lengths in this plot. Might conflict across
              % Instances.
              
              % Use those times to extract the desired slice, preparing it
              % for later combining and for plotting.
              activitySlice = stimPSTH{stim_i, group_i, strcmp(dataType, 'PSTH')}(:, startT:endT);
              actSize = size(activitySlice);
              newActInd = find(isnan(eventPSTHData{event_i, group_i, strcmp(dataType3, 'PSTH')}(:,1)), 1);
              eventPSTHData{event_i, group_i, strcmp(dataType3, 'PSTH')}(newActInd:newActInd + actSize(1)-1, 1:actSize(2)) = activitySlice;
              
              % Gather some number of traces for the null distribution
              eventPSTHNullDataTmp = vertcat(stimPSTH{~stimLogicArrayEvent, group_i, strcmp(dataType, 'PSTH')});             % Generate the null trace, store null values for later.
              trace2Gather = min(size(activitySlice,1)*4, size(eventPSTHNullDataTmp,1));
              eventPSTHNullDataTmp = eventPSTHNullDataTmp(1:trace2Gather, startT:endT);
              actSize = size(eventPSTHNullDataTmp);
              newActInd = find(isnan(eventPSTHData{event_i, group_i, strcmp(dataType3, 'PSTHnull')}(:,1)), 1);
              eventPSTHData{event_i, group_i, strcmp(dataType3, 'PSTHnull')}(newActInd:newActInd + actSize(1)-1, 1:actSize(2)) = eventPSTHNullDataTmp;
            end
          end
        end
      end
    end
    
    % Remove Excess NaNs, horizontally and vertically.
    eventPSTHData = cellfun(removeExcessNaNs, eventPSTHData(1:end-length(comboEvents), :, :), 'UniformOutput', false);
    eventPSTHData(cellfun('isempty', eventPSTHData)) = deal({zeros(0,0)});
    
    % Identify which events have data to plot
    comboInd1 = strcmp(eventListPlotCombo, 'headTurn right') | strcmp(eventListPlotCombo, 'headTurn left');
    comboInd2 = comboInd1 | strcmp(eventListPlotCombo, 'bodyTurn');
    eventPlotInd = all(cellfun('isempty', eventPSTHData), 2:3);
    events2Plot = ~[eventPlotInd; any(eventPlotInd(comboInd1)); any(eventPlotInd(comboInd2))]';
    
    for event_i = find(events2Plot)
      for group_i = 1:length(groupIterInd)
        
        if length(groupIterInd{group_i}) ~= 1
          unitTag = strjoin(groupTypeArray(groupIterInd{group_i}), ' & ');
        else
          unitTag = groupTypeArray{groupIterInd{group_i}};
        end
        
        figTitle = sprintf('Event aligned PSTH - %s %s', unitTag, eventListPlotCombo{event_i});
        h = figure('NumberTitle', 'off', 'Name', figTitle,'units','normalized','outerposition', [0.15 0.22 0.44 0.54]);
        
        %comboInd means the set is a combination of previous entries.
        plotTraceData = vertcat(eventPSTHData{comboGroups{event_i}, groupIterInd{group_i}, 1});
        plotTraceNull = vertcat(eventPSTHData{comboGroups{event_i}, groupIterInd{group_i}, 2});
        
        % Prepare plot data
        plotData = mean(plotTraceData, 1);
        plotErr = std(plotTraceData)/sqrt(size(plotTraceData,1));
        
        plotNull = mean(plotTraceNull, 1);
        plotErrNull = std(plotTraceNull)/sqrt(size(plotTraceNull,1));
        
        % Subtract baseline if desired
        baselineSubtract = 1;
        if baselineSubtract
          baselinePeriod = (psthParam.psthPre/2):psthParam.psthPre;
          plotData = plotData - mean(plotData(baselinePeriod));
          plotNull = plotNull - mean(plotNull(baselinePeriod));
        end
        
        msebData = [plotData; plotNull];
        msebError = [plotErr; plotErrNull];
%         plotLabels = [eventListTitle(event_i); sprintf('non-%s', eventListTitle{event_i})];
        plotLabels = [eventListPlotCombo(event_i); 'Scrambled Event Times'];
        
%         subplot(length(groupIterInd), 1, group_i)
        [psthAxes, ~, ~, legendH] = plotPSTH(msebData, msebError, [], params.subEventPSTHParams, 'line', [], plotLabels);
        
        title(figTitle);
        psthAxes.FontSize = 16;
        legendH.Location = 'northwest';
        xlabel('Time from Event Onset (ms)');
        ylabel('Normalized Activity');
        
        saveFigure(outputDir, sprintf('4. %s - %s', paradigmList{par_i}, figTitle), [], figStruct, []);
        
      end
      
      
    end
    
  end
  
end

% Put desired outputs in the struct.
subEventPSTHStruct = struct();

end

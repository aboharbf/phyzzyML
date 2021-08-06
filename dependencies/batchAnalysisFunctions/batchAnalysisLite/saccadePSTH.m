function saccadePSTH(spikePathBank, params)
% A function which retrieves all the saccade

smoothTraces = 1;
smoothingWidth = 50;

collapseBins = 1;    % Triggers the 'binStuff' function on the saved data, which simply rebins to a larger window
binningSize = 100;       % Defines the size of the larger window.
binStep = 100;       % Defines the step of the window.
load(params.subEventPSTHParams.eventData, 'eventData'); %loads 'eventData' var.

outputDir = fullfile(params.outputDir, 'saccadeDynamics');

% Collect unit selectivity
paradigmList = unique(spikePathBank.paradigmName);
monkeyList = {'Sam', 'Mo', 'Combo'};

% Generate smoothing filter
if smoothTraces
  filterPoints = -3*smoothingWidth:3*smoothingWidth;
  smoothingFilter = exp(-1*filterPoints.^2/(2*smoothingWidth^2));
  smoothingFilter = smoothingFilter/sum(smoothingFilter);
  smoothingWidthNom = (smoothingWidth * 3); % * (binningSize * collapseBins); % Nominal smoothing is larger if bins are collapsed.
  smoothTag = sprintf('(Smoothed, %d ms)', smoothingWidthNom);
else
  smoothTag = '';  
end

if collapseBins
  colTag = sprintf('(Binned, %d ms)', binningSize);
else
  colTag = '';
end

for monkey_i = 1:length(monkeyList)
  for par_i = 1:length(paradigmList)
    
    pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
    if ~strcmp(monkeyList{monkey_i}, 'Combo')
      mInd = contains(spikePathBank.analyzedDir, monkeyList{monkey_i});
      spikePathBankParadigm = spikePathBank(pInd & mInd, :);
    else
      spikePathBankParadigm = spikePathBank(pInd, :);
    end
    
%     spikePathBankParadigm = spikePathBank(pInd, :);
    
    % Grab necessary data from each run.
    [saccadeStackParamsPerRun, saccadeMatArrayPerRun, saccadeMatLabelPerRun] = spikePathLoad(spikePathBankParadigm, {'saccadeStackParams', 'saccadeMatArray', 'saccadeMatLabel'}, params.spikePathLoadParams);
    
    % Find Unique entries
    [uniqueSaccadeEvents, firstInd] = unique(vertcat(saccadeMatLabelPerRun{:}));
    tmp = vertcat(saccadeMatArrayPerRun{:});
    firstDataEntry = tmp(firstInd);
    eventLengths = cellfun('length', firstDataEntry);
    uniqueLengths = unique(eventLengths);
    
    % Initialize Storage vectors, and some quick parameters for the
    % plotPSTH function.
    [plotDataStacks, plotLabelStacks, plotParamStack, plotTitleStack] = deal(cell(length(uniqueLengths),1));
    shorterParams.psthPre = saccadeStackParamsPerRun{1}.preEventTime;
    shorterParams.psthImDur = saccadeStackParamsPerRun{1}.postEventTime;
    shorterParams.psthPost = 0;
    plotParamStack{1} = shorterParams;
    plotTitleStack{1} = sprintf('Task Event Saccade Mean');
    
    longerParams.psthPre = saccadeStackParamsPerRun{1}.preStimTime;
    longerParams.psthImDur = 2800;
    longerParams.psthPost = saccadeStackParamsPerRun{1}.postStimTime;
    plotParamStack{2} = longerParams;
    plotTitleStack{2} = sprintf('%s Stimulus Saccade Mean %s %s', monkeyList{monkey_i}, smoothTag, colTag);
    
    for len_i = 1:length(uniqueLengths)
      events2Plot = eventLengths == uniqueLengths(len_i);
      eventLabels2Plot = uniqueSaccadeEvents(events2Plot);
      crossEventTraceStack = cell(size(eventLabels2Plot));
      
      for event_i = 1:length(eventLabels2Plot)
        % For each event, collect all the traces from each run corresponding
        % to it.
        eventTraceStack = [];
        for run_i = 1:length(saccadeMatArrayPerRun)
          traceInd = strcmp(saccadeMatLabelPerRun{run_i}, eventLabels2Plot{event_i});
          eventTraceStack = [eventTraceStack; saccadeMatArrayPerRun{run_i}{traceInd}];
        end
        
        % Take the mean across all instances of the event
        crossEventTraceStack{event_i} = eventTraceStack;
        
      end
      
      plotDataStacks{len_i} = crossEventTraceStack;
      plotLabelStacks{len_i} = eventLabels2Plot;
      
    end
    
    % Collect another form of data which isn't in the initial structures -
    % subEvent data
    stimuliNames = plotLabelStacks{uniqueLengths > 2800};
    stimuliData = plotDataStacks{uniqueLengths > 2800};
    eventDataStim = eventData(stimuliNames, :);
    eventsInd = any(cellfun(@(x) size(x,1) ~= 0, table2cell(eventDataStim)));
    eventNames = eventDataStim.Properties.VariableNames(eventsInd);
    eventDataStim = eventDataStim(:, eventsInd);
    
    % Initialize data
    subEventData = cell(size(eventNames'));
    plotParamStack = [plotParamStack; plotParamStack(1)];
    window2Pull = [-plotParamStack{1}.psthPre:plotParamStack{1}.psthImDur] + plotParamStack{2}.psthPre;
    for event_i = 1:length(eventNames)
      eventTimes = eventDataStim{:, eventNames{event_i}};
      stimWEvent = find(cellfun(@(x) size(x, 1) ~= 0, eventTimes));
      for stim_i = stimWEvent'
        eventTimeStim = round(eventTimes{stim_i}.startTime);
        for ev_i = 1:length(eventTimeStim)
          subEventData{event_i} = [subEventData{event_i}; stimuliData{stim_i}(:, window2Pull + eventTimeStim(ev_i))];
        end
      end
      
    end
    
    % Stack, and add necessary cells to arrays
    plotDataStacks =  [plotDataStacks; {subEventData}];
    plotLabelStacks{3} = strrep(eventNames, '_', ' ')';
    plotTitleStack{3} = sprintf('subEvent Saccade Mean');
    
    % Turn the arrays into single mean traces, with counts
    plotDataCounts = cellfun(@(x1) cellfun(@(x2) size(x2,1), x1), plotDataStacks, 'UniformOutput', false);
    plotDataStacks2 = cellfun(@(x1) cellfun(@(x2) mean(x2,1), x1, 'UniformOutput', false), plotDataStacks, 'UniformOutput', false);
    plotDataStacks = cellfun(@(x) vertcat(x{:}), plotDataStacks2, 'UniformOutput', false);
    
    % Post processing of the data - Smooth first.
    if smoothTraces
      for ii = 1:length(plotDataStacks)
        for jj = 1:size(plotDataStacks{ii},1)
          plotDataStacks{ii}(jj,:) = conv(mean(plotDataStacks{ii}(jj,:), 1), smoothingFilter, 'same');
        end
      end
    end
    
    % Then bin
    if collapseBins
      plotDataStacks = cellfun(@(x) binStuff(x, binningSize, binStep), plotDataStacks, 'UniformOutput', false);
    end
    
    % Plot the data
    figTitle = sprintf('Saccade Dynamics - %s - %s %s %s', monkeyList{monkey_i}, paradigmList{par_i}, colTag, smoothTag);
    h = figure('NumberTitle', 'off', 'Name', figTitle,'units','normalized','outerposition', [0.0042 0.0074 0.6250 0.9139]);
    subplotSize = [4, {1:3}, 5];
    for plot_i = 1:length(plotDataStacks)
      % Select plots
      subplot(5 , 1, subplotSize{plot_i})
      
      % Process labels to make them look better
      [plotLabels, resortInd] = resortAndTrimLabels(plotLabelStacks{plot_i});
      
      % Plot
      psthAx = plotPSTH(plotDataStacks{plot_i}(resortInd, :), [], [], plotParamStack{plot_i}, 'color', plotTitleStack{plot_i}, plotLabels);
      
      % This label collides w/ the title below.
      if plot_i ~= 3
        psthAx.XLabel.String = '';
      end
      
    end
    
    % Save the figure
    saveFigure(outputDir, figTitle, [], params.figStruct, [])
    
  end
end

end
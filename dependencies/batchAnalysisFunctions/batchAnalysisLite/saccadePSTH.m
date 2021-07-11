function saccadePSTH(spikePathBank, params)
% A function which retrieves all the saccade

smoothingWidth = 1;
smoothTraces = true;
collapseBins = true;
binningSize = 30;

% Collect unit selectivity
paradigmList = unique(spikePathBank.paradigmName);

% Generate smoothing filter
if smoothTraces
  filterPoints = -3*smoothingWidth:3*smoothingWidth;
  smoothingFilter = exp(-1*filterPoints.^2/(2*smoothingWidth^2));
  smoothingFilter = smoothingFilter/sum(smoothingFilter);
  smoothTag = sprintf('(Smoothing %d ms)', smoothingWidth);
else
  smoothTag = '';  
end

for par_i = 1:length(paradigmList)
  
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  spikePathBankParadigm = spikePathBank(pInd, :);
  
  [saccadeStackParamsPerRun, saccadeMatArrayPerRun, saccadeMatLabelPerRun] = spikePathLoad(spikePathBankParadigm, {'saccadeStackParams', 'saccadeMatArray', 'saccadeMatLabel'}, params.spikePathLoadParams);
  
  % Find Unique entries
  [uniqueSaccadeEvents, firstInd] = unique(vertcat(saccadeMatLabelPerRun{:}));
  tmp = vertcat(saccadeMatArrayPerRun{:});
  firstDataEntry = tmp(firstInd);
  eventLengths = cellfun('length', firstDataEntry);
  uniqueLengths = unique(eventLengths);
  
  [plotDataStacks, plotLabelStacks, plotParamStack, plotTitleStack] = deal(cell(length(uniqueLengths),1));
  
  shorterParams.psthPre = saccadeStackParamsPerRun{1}.preEventTime;
  shorterParams.psthImDur = saccadeStackParamsPerRun{1}.postEventTime;
  shorterParams.psthPost = 0;
  plotParamStack{1} = shorterParams;
  plotTitleStack{1} = sprintf('Task Event Saccade Mean %s', smoothTag);
  
  longerParams.psthPre = saccadeStackParamsPerRun{1}.preStimTime;
  longerParams.psthImDur = 2800;
  longerParams.psthPost = saccadeStackParamsPerRun{1}.postStimTime;
  plotParamStack{2} = longerParams;
  plotTitleStack{2} = sprintf('Stimulus Saccade Mean %s', smoothTag);

  for len_i = 1:length(uniqueLengths)
    events2Plot = eventLengths == uniqueLengths(len_i);
    eventLabels2Plot = uniqueSaccadeEvents(events2Plot);
    
    if collapseBins
      crossEventTraceStack = nan(sum(events2Plot), ceil(uniqueLengths(len_i)/binningSize));
    else
      crossEventTraceStack = nan(sum(events2Plot), uniqueLengths(len_i));
    end
    
    for event_i = 1:length(eventLabels2Plot)
      % For each event, collect all the traces from each run corresponding
      % to it.
      eventTraceStack = [];
      for run_i = 1:length(saccadeMatArrayPerRun)
        traceInd = strcmp(saccadeMatLabelPerRun{run_i}, eventLabels2Plot{event_i});
        eventTraceStack = [eventTraceStack; saccadeMatArrayPerRun{run_i}{traceInd}];
      end
      
      if collapseBins
        eventTraceStack = binStuff(eventTraceStack, binningSize);
      end
      
      if smoothTraces
        crossEventTraceStack(event_i,:) = conv(mean(eventTraceStack, 1), smoothingFilter, 'same');
      else
        crossEventTraceStack(event_i,:) = mean(eventTraceStack, 1);
      end
    end
    
    plotDataStacks{len_i} = crossEventTraceStack;
    plotLabelStacks{len_i} = eventLabels2Plot;
    
  end
  
  figure()
  subplotSize = [4, {1:3}];
  for plot_i = 1:length(plotDataStacks)
    subplot(4,1, subplotSize{plot_i})
    plotPSTH(plotDataStacks{plot_i}, [], [], plotParamStack{plot_i}, 'color', plotTitleStack{plot_i}, plotLabelStacks{plot_i});
  end
  
end

disp('y')

end
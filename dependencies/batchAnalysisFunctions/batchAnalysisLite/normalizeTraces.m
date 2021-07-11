function subEventSigStructPerRunNormed = normalizeTraces(subEventSigStructPerRun, psthByImagePerRunUnProc, stimTimingByRun, spikeAlignParamsByRun)
% Cycle through each subEventSigStruct and modify the traces to be
% normalized.

% Use the taskData and spikeAlignParams to calculate the window of
% activity to sample inpsthByImageFixAlign, then clear
ITIwindow = cell(size(spikeAlignParamsByRun));
for run_i = 1:length(spikeAlignParamsByRun)
  ITIwindow{run_i} = spikeAlignParamsByRun{run_i}.preAlign - stimTimingByRun{run_i}.ISI : spikeAlignParamsByRun{run_i}.preAlign;
  if stimTimingByRun{run_i}.ISI == 0
    % in circumstances where ISI is set to 0 (which the computer does not actually achieve) take a slice of 200 ms instead.
    ITIwindow{run_i} = spikeAlignParamsByRun{run_i}.preAlign - 200 : spikeAlignParamsByRun{run_i}.preAlign;
  end
end

subEventSigStructPerRunNormed = cell(size(subEventSigStructPerRun));
for run_i = 1:length(subEventSigStructPerRun)
  
  % Pull Data
  psthData = psthByImagePerRunUnProc{run_i};
  subEventData = subEventSigStructPerRun{run_i};
  
  % Find the first 500 ms of each available trace
  for chan_i = 1:length(psthData)
    for unit_i = 1:length(psthData{chan_i})
      
      % find the baseline across traces
      baselineData = psthByImagePerRunUnProc{run_i}{chan_i}{unit_i}(:, ITIwindow{run_i});
      meanData = mean(baselineData(:));
      stdData = std(baselineData(:));
      
      % Ignore if mean during background is zero, unit is likely very low
      % firing, consider thresholding to remove such instances.
      if stdData == 0
        continue
      end
      
      % Modify the traces
      subEventData.subEventPSTH{chan_i,1}{unit_i} = (subEventData.subEventPSTH{chan_i,1}{unit_i} - meanData)/stdData;
      subEventData.subEventNullPSTH{chan_i,1}{unit_i} = (subEventData.subEventNullPSTH{chan_i,1}{unit_i} - meanData)/stdData;
  
    end
  end
  
  subEventSigStructPerRunNormed{run_i} = subEventData;

end

end
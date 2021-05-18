function subEventSigStructPerRunNormed = normalizeTraces(subEventSigStructPerRun, psthByImagePerRunUnProc)
% Cycle through each subEventSigStruct and modify the traces to be
% normalized.

subEventSigStructPerRunNormed = cell(size(subEventSigStructPerRun));
for run_i = 1:length(subEventSigStructPerRun)
  
  % Pull Data
  psthData = psthByImagePerRunUnProc{run_i};
  subEventData = subEventSigStructPerRun{run_i};
  
  % Find the first 500 ms of each available trace
  for chan_i = 1:length(psthData)
    for unit_i = 1:length(psthData{chan_i})
      
      % find the baseline across traces
      baselineData = psthData{chan_i}{unit_i}(:, 1:400);
      meanData = mean(baselineData(:));
      stdData = std(baselineData(:));
      
      % Modify the traces
      subEventData.subEventPSTH{chan_i,1}{unit_i} = (subEventData.subEventPSTH{chan_i,1}{unit_i} - meanData)/stdData;
      subEventData.subEventNullPSTH{chan_i,1}{unit_i} = (subEventData.subEventNullPSTH{chan_i,1}{unit_i} - meanData)/stdData;
  
    end
  end
  
  subEventSigStructPerRunNormed{run_i} = subEventData;

end

end
function [psthByImagePerRun, eventDataArray, dataTypeArray, groupTypeArray, eventList, subEventSig, psthParams2Use] = generateEventDataArray(spikePathBankParadigm, params)
% Extracts subEvent data across all the runs.

% Decide which variables to extract from the runs (Any amount of
% Proc/Normalizing should be considered).
variables2Extract = {'psthByImage', 'psthParams', 'spikeAlignParams', 'stimTiming'};
if params.meanPSTHParams.normalize
  variables2Extract(1) = strcat(variables2Extract(1), 'Proc');
end

paradigmRunList = extractAfter(spikePathBankParadigm.Row, 'S');
variables2Extract = [variables2Extract, 'subEventSigStruct'];

[psthByImagePerRun, psthParamsPerRun, spikeAlignParamsByRun, stimTimingByRun, subEventSigStructPerRun] = spikePathLoad(spikePathBankParadigm, variables2Extract, params.spikePathLoadParams);

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

allSubEventSigStructs = vertcat(subEventSigStructPerRun{:});
eventList = unique(vertcat(allSubEventSigStructs(:).events))';
psthParams2Use = psthParamsPerRun{1};

% If you want to normalize these subEvents, you'll need unprocessed forms
% of the psth from every run
if params.subEventPSTHParams.normalize
  psthByImagePerRunFixAlignedUnProc = spikePathLoad(spikePathBankParadigm, {'psthByImageFixAlign'}, params.spikePathLoadParams);
end

% Gather data across all the spikeDataBank Structure, generating a cell
% array with the indicies eventDataArray{event_i}{group_i}{data_i}
% event_i = eventList
groupTypeArray = {'Unsorted', 'Unit', 'MUA'};
dataTypeArray = {'PSTH', 'Null PSTH', 'pVal', 'Cohens D', 'Run Index'};
eventDataArray = cell(length(eventList), length(groupTypeArray), length(dataTypeArray));

tmpFile = fullfile(params.subEventPSTHParams.outputDir, sprintf('%s_tempFile.mat', spikePathBankParadigm.paradigmName{1}));

if ~exist(tmpFile, 'file')
  
  if ~exist(params.subEventPSTHParams.outputDir, 'dir')
    mkdir(params.subEventPSTHParams.outputDir);
  end
  
  tic
  for run_i = 1:length(paradigmRunList)
    subEventSig = subEventSigStructPerRun{run_i};
    subEventData = subEventSig.subEventPSTH(:,1);
    subEventNull = subEventSig.subEventNullPSTH(:,1);
    for event_i = 1:length(subEventSig.events)
      eventInd = strcmp(subEventSig.events{event_i}, eventList);
      for chan_i = 1:length(subEventData)
        for unit_i = 1:length(subEventData{chan_i})
          % Grouping Index
          if unit_i == 1
            groupInd = 1;
          elseif unit_i == length(subEventData{chan_i})
            groupInd = 3;
          else
            groupInd = 2;
          end
          
          % Pull data from struct.
          PSTHData = subEventData{chan_i}{unit_i}(event_i,:);
          if sum(PSTHData) < 1
            continue % Units which are empty
          end
          PSTHNull = subEventNull{chan_i}{unit_i}(event_i,:);
          pVal = subEventSig.testResults{chan_i}{unit_i}{event_i};
          cohensD = subEventSig.cohensD{chan_i}{unit_i}{event_i};
          runLabel = paradigmRunList(run_i);
          
          % Normalize with respect to baseline from all the stimuli.
          if params.subEventPSTHParams.normalize
            preFixActivity = psthByImagePerRunFixAlignedUnProc{run_i}{chan_i}{unit_i}(:,ITIwindow{run_i});
            preFixActivity = reshape(preFixActivity, [], 1);
            
            preFixMean = mean(preFixActivity);
            preFixSD = std(preFixActivity);
            
            if preFixMean ~= 0 && preFixSD ~= 0
              PSTHData = (PSTHData - preFixMean)/preFixSD;
              PSTHNull = (PSTHNull - preFixMean)/preFixSD;
            end
          end
          
          % Sort these into an array correctly.
          newEntry = cell(size(dataTypeArray));
          newEntry{strcmp(dataTypeArray, 'PSTH')} = PSTHData;
          newEntry{strcmp(dataTypeArray, 'Null PSTH')} = PSTHNull;
          newEntry{strcmp(dataTypeArray, 'pVal')} = pVal;
          newEntry{strcmp(dataTypeArray, 'Cohens D')} = cohensD;
          newEntry{strcmp(dataTypeArray, 'Run Index')} = runLabel;
          
          % Compile relevant data from run
          for data_i = 1:length(dataTypeArray)
            eventDataArray{eventInd, groupInd, data_i} = [eventDataArray{eventInd, groupInd, data_i}; newEntry{data_i}];
          end
        end
      end
    end
  end
  
  % save a temp file
  save(tmpFile, 'eventDataArray', 'subEventSig')
  fprintf('Done compiling subEvent data: %s s \n', num2str(toc))
  
else
  load(tmpFile, 'eventDataArray', 'subEventSig')  
end

end
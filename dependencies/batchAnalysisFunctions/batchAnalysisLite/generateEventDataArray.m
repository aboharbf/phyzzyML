function [psthByImagePerRun, eventDataArray, dataTypeArray, groupTypeArray, eventList, subEventSig, psthParams2Use] = generateEventDataArray(spikePathBankParadigm, params)
% Extracts subEvent data across all the runs.

% Decide which variables to extract from the runs (Any amount of
% Proc/Normalizing should be considered).
variables2Extract = {'psthByImage', 'psthParams'};
if params.meanPSTHParams.normalize
  variables2Extract(1) = strcat(variables2Extract(1), 'Proc');
end

paradigmRunList = extractAfter(spikePathBankParadigm.Row, 'S');
variables2Extract = [variables2Extract, 'subEventSigStruct'];

[psthByImagePerRun, psthParamsPerRun, subEventSigStructPerRun] = spikePathLoad(spikePathBankParadigm, variables2Extract, params.spikePathLoadParams);

allSubEventSigStructs = vertcat(subEventSigStructPerRun{:});
eventList = unique(vertcat(allSubEventSigStructs(:).events))';
psthParams2Use = psthParamsPerRun{1};

% If you want to normalize these subEvents, you'll need unprocessed forms
% of the psth from every run
if params.subEventPSTHParams.normalize
  psthByImagePerRunUnProc = spikePathLoad(spikePathBankParadigm, {'psthByImage'}, params.spikePathLoadParams);
end

% Gather data across all the spikeDataBank Structure, generating a cell
% array with the indicies eventDataArray{event_i}{group_i}{data_i}
% event_i = eventList
groupTypeArray = {'Unsorted', 'Unit', 'MUA'};
dataTypeArray = {'PSTH', 'Null PSTH', 'pVal', 'Cohens D', 'Run Index'};
eventDataArray = cell(length(eventList), length(groupTypeArray), length(dataTypeArray));

tmpFile = fullfile(params.subEventPSTHParams.outputDir, sprintf('%s_tempFile.mat', spikePathBankParadigm.paradigmName{1}));

if ~exist(tmpFile, 'file')
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
            fixStart = 1; %params.fixBuffer;
            fixEnd = psthParamsPerRun{run_i}.ITI;
            fixActivity = psthByImagePerRunUnProc{run_i}{chan_i}{unit_i}(:,fixStart:fixEnd);
            fixActivity = reshape(fixActivity, [], 1);
            
            fixMean = mean(fixActivity);
            fixSD = std(fixActivity);
            
            if fixMean ~= 0 && fixSD ~= 0
              PSTHData = (PSTHData - fixMean)/fixSD;
              PSTHNull = (PSTHNull - fixMean)/fixSD;
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
  fprintf('Done compiling subEvent data: %d s \n', num2str(toc))
  
else
  load(tmpFile, 'eventDataArray', 'subEventSig')  
end

end
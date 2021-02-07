function  spikePathBank_to_rasterData(spikePathBank, rasterDataPath, params)
% Turns the spikeDataBank struct of processBatchAnalysis into the format
% needed to process by Neural Decoding Toolbox.

% Creates file with the following 3 fields
% raster_data - 1 * unit cell array, each cell containing trial * bins.
% raster_labels - struct with each field contain 1*trial labels. Labels
% like 'stimulus' or 'position' go here. These are currently generated
% based on the stimParamFile used by phyzzy combined with plotIndex. The
% labels and their fieldNames are defined in params.plotIndParams.
% raster_site_info - file*1 labels/numbers - session_ID, channel, unit,
% site.

%runData.spikesByEventBinned{stim}{chan}{unit} --> needs to be turned to
%rasterDataUnit

% Load the eventData for attaching significance activity to units
if exist(params.subEventBatchStructPath, 'file')
  tmp = load(params.subEventBatchStructPath);
  sigUnitGrid = tmp.subEventBatchStruct.sigUnitGrid;
  eventList = tmp.subEventBatchStruct.eventList;
end

paradigmList = unique(spikePathBank.paradigmName);
paradigmList = paradigmList(2:end);

for par_i = 1:length(paradigmList)
  rasterDataParadigmPath = fullfile(rasterDataPath, paradigmList{par_i});
  
  if ~exist(rasterDataParadigmPath, 'dir')
    mkdir(rasterDataParadigmPath)
  end

  paradigmSPInd = strcmp(spikePathBank.paradigmName, paradigmList(par_i)); 
  spikePathBankParadigm = spikePathBank(paradigmSPInd,:);
  
  % Load the spike data to be turned into raster data.
  [spikesByEventBinned, psthParams] = spikePathLoad(spikePathBankParadigm, {'spikesByEventBinned', 'psthParams'}, params.spikePathLoadParams);
  runList = spikePathBankParadigm.Properties.RowNames;
  
  for run_i = 1:length(runList)
    
    % print a message the the data is being binned (and add a dot for each file that has been binned
    curr_bin_string = ['Processing Runs: ' num2str(run_i) ' of ' num2str(length(runList))];
    if run_i == 1
      disp(curr_bin_string);
    else
      fprintf([repmat(8,1,bin_str_len) curr_bin_string]);
    end
    bin_str_len = length(curr_bin_string);
    
    % Collect data for this particular run
    runData = struct();
    runData.start = -psthParams{run_i}.psthPre;
    runData.end = psthParams{run_i}.psthImDur + psthParams{run_i}.psthPost;
    
    runData.spikesByEventBinned = spikesByEventBinned{run_i};
    runData.eventIDs = spikePathBankParadigm.stimuli{run_i};
    
    %   runData.stimPresCount = 1;
    %   runData.gridHoles = 1;
    %   runData.recDepth = 1;
    
    % Indicies for activity, accounting for padding.
    dataLength = abs(runData.start) + runData.end;
    vecLength = size(runData.spikesByEventBinned{1}{1}{1},2);
    if dataLength ~= vecLength
      startTime = floor((vecLength - dataLength)/2);
      endTime = dataLength+startTime-1;
    else
      startTime = 1;
      endTime = dataLength;
    end
    
    % find out trials per stimuli to generate useful indicies
    trialsPerStim = cellfun(@(x) size(x{1}{1},1), runData.spikesByEventBinned);
    trialsPerStimStart = cumsum([1; trialsPerStim(1:end-1)]);
    trialsPerStimEnd = cumsum(trialsPerStim);
    trialsPerUnit = sum(trialsPerStim);
    stimCount = 1:length(trialsPerStim);
    
    % Generate the raster file fields which are the same for every unit.
    % raster_labels - labels that exist per trial. generate them per stimuli,
    % then repmat them.
    
    % Generate a vector for stimuli for each trial.
    stimuliVecTmp = cellfun(@(x) convertCharsToStrings(x), runData.eventIDs);
    stimuliVecTmp = arrayfun(@(x) repmat(stimuliVecTmp(x), [trialsPerStim(x), 1]), stimCount, 'UniformOutput', 0)';
    stimuliVec = vertcat(stimuliVecTmp{:});
    
    % Generate a vector for raster labels, with 0/1 for binary choices, or an
    % index for a larger set.
    rasterParams = params.(paradigmList{par_i});
    params.plotIndParams.plotLabels = rasterParams.plotIndParams.plotLabels;
    rasterLabelTmp = plotIndex(runData.eventIDs, params.plotIndParams);
    rasterLabelTmp = arrayfun(@(x) repmat(rasterLabelTmp(x, :), [trialsPerStim(x), 1]), stimCount, 'UniformOutput', 0)';
    rasterLabel = vertcat(rasterLabelTmp{:});
    
    % Generate a vector of stimulus presentation counts for a particular run.
    %   stimPresCountTmp = arrayfun(@(x) repmat(runData.stimPresCount(x), [trialsPerStim(x), 1]), stimCount, 'UniformOutput', 0)';
    %   stimPresCountPerTrial = vertcat(stimPresCountTmp{:});
    
    % Assign to raster labels.
    raster_labels = struct();
    raster_labels.stimuli = stimuliVec';
    %   raster_labels.stimPresCount = stimPresCountPerTrial';
    
    % iterate through rasterLabels.
    for label_i = 1:length(rasterParams.rasterLabels)
      % Convert the indices into strings
      label = rasterParams.plotIndParams.plotLabels{label_i};
      entry = cell(size(rasterLabel,1),1);
      
      if iscell(label)
        % Use entries as indices to the provided labels
        entryInd = rasterLabel(:,label_i);
        for lab_i = 1:length(label)
          entry(entryInd == lab_i) = label(lab_i);
        end
      else
        % use single label, add 'Non' to make 2nd.
        entryInd = logical(rasterLabel(:,label_i));
        entry(entryInd) = {label};
        entry(~entryInd) = {['Non-' label]};
      end
      
      % Store them
      raster_labels.(rasterParams.rasterLabels{label_i}) = entry;
    end
    
    % raster_site_info
    raster_site_info.session_ID = runList{run_i};
    raster_site_info.alignment_event_time = abs(runData.start);
%     raster_site_info.runNum = runData.runNum;
    
    % Construct raster_data and save file per unit.
    for chan_i = 1:length(runData.spikesByEventBinned{1})
      % Information which is consistent per channel
      %     raster_site_info.gridHole = num2str(runData.gridHoles{chan_i});
      %     raster_site_info.recordingDepth = runData.recDepth{chan_i};
      
      for unit_i = 1:length(runData.spikesByEventBinned{1}{chan_i})
        
        % Cycle through the stim, store into raster_data.
        raster_data = zeros(trialsPerUnit, dataLength);
        for stim_i = 1:length(runData.spikesByEventBinned)
          raster_data(trialsPerStimStart(stim_i):trialsPerStimEnd(stim_i),:) = runData.spikesByEventBinned{stim_i}{chan_i}{unit_i}(:, startTime:endTime);
        end
        
        % Generate the raster_site_info
        raster_site_info.UnitType = convertUnitToName(unit_i, length(runData.spikesByEventBinned{1}{chan_i}), 2);
        
        ULabel = convertUnitToName(unit_i, length(runData.spikesByEventBinned{1}{chan_i}), 1);
        chLabel = ['Ch' num2str(chan_i)];
        rasterFileName = [runList{run_i},chLabel, ULabel];
        
        % use generated label to look up eventData, subEvent significance.
        if exist('sigUnitGrid', 'var')
          gridRow = strcmp(sigUnitGrid.sigUnitLabels, rasterFileName(2:end));
          subStructSigData = sigUnitGrid.sigUnitMat(gridRow, :);
          for event_i = 1:length(eventList)
            raster_site_info.(eventList{event_i}) = subStructSigData(event_i);
          end
          
          % Combo events
          comboName = {'headTurn_all', 'allTurn'};
          comboInts = {{'headTurn_right','headTurn_left'};{'headTurn_right', 'headTurn_left', 'bodyTurn'}};
          
          for combo_i = 1:length(comboName)
            anyArray = 0;
            for event_i = 1:length(comboInts{combo_i})
              anyArray = anyArray | subStructSigData(strcmp(eventList', comboInts{combo_i}{event_i}));
            end
            raster_site_info.(comboName{combo_i}) = anyArray;
          end
        end
        
        
        % Save file in output folder
        save(fullfile(rasterDataParadigmPath, rasterFileName), 'raster_data', 'raster_labels', 'raster_site_info')
        
      end
    end
    
  end
end



end
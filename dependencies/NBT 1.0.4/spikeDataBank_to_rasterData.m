function  spikeDataBank_to_rasterData(spikeDataBank, rasterDataPath, params)
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

runList = fields(spikeDataBank);

if ~exist(rasterDataPath, 'dir')
  mkdir(rasterDataPath)
end

%runData.spikesByEventBinned{stim}{chan}{unit} --> needs to be turned to
%rasterDataUnit

for run_i = 1:length(runList)
  
  % print a message the the data is being binned (and add a dot for each file that has been binned
  curr_bin_string = ['Processing Runs: ' num2str(run_i) ' of ' num2str(length(runList))];
  if run_i == 1
    disp(curr_bin_string);
  else
    fprintf([repmat(8,1,bin_str_len) curr_bin_string]);
  end
  bin_str_len = length(curr_bin_string);
  
  runData = spikeDataBank.(runList{run_i});
  
  % Initialize some variables which make the rest easier to gather.
  % Indicies for activity
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
  
  % Generate the raster file fields which are the same for every unit.
  % raster_labels - labels that exist per trial. generate them per stimuli,
  % then repmat them.

  stimuliVecTmp = cellfun(@(x) convertCharsToStrings(x), runData.eventIDs);
  rasterLabelTmp = plotIndex(runData.eventIDs, params.plotIndParams);  
  
  stimuliVecTmp = arrayfun(@(x) repmat(stimuliVecTmp(x), [trialsPerStim(x), 1]), 1:length(trialsPerStim), 'UniformOutput', 0)';
  rasterLabelTmp = arrayfun(@(x) repmat(rasterLabelTmp(x, :), [trialsPerStim(x), 1]), 1:length(trialsPerStim), 'UniformOutput', 0)';
  
  stimPresCountTmp = arrayfun(@(x) repmat(runData.stimPresCount(x), [trialsPerStim(x), 1]), 1:length(trialsPerStim), 'UniformOutput', 0)';
  stimPresCountPerTrial = vertcat(stimPresCountTmp{:});
  
  stimuliVec = vertcat(stimuliVecTmp{:});
  rasterLabel = vertcat(rasterLabelTmp{:});
  
  raster_labels.stimuli = stimuliVec';
  raster_labels.stimPresCount = stimPresCountPerTrial';
  
  for label_i = 1:length(params.rasterLabels)
    % Convert the indices into strings
    label = params.plotIndParams.plotLabels{label_i};
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
    raster_labels.(params.rasterLabels{label_i}) = entry;
  end
  
  % raster_site_info
  raster_site_info.session_ID = runList{run_i};
  raster_site_info.alignment_event_time = abs(runData.start);
  raster_site_info.runNum = runData.runNum;
  
  % Construct raster_data and save file per unit.
  for chan_i = 1:length(runData.spikesByEventBinned{1})
    % Information which is consistent per channel
    raster_site_info.gridHole = num2str(runData.gridHoles{chan_i});
    raster_site_info.recordingDepth = runData.recDepth{chan_i};
    
    for unit_i = 1:length(runData.spikesByEventBinned{1}{chan_i})
      
      % Cycle through the stim, store into raster_data.
      raster_data = zeros(trialsPerUnit, dataLength);
      for stim_i = 1:length(runData.spikesByEventBinned)
        raster_data(trialsPerStimStart(stim_i):trialsPerStimEnd(stim_i),:) = runData.spikesByEventBinned{stim_i}{chan_i}{unit_i}(:, startTime:endTime);
      end
      
      % Generate the raster_site_info
      raster_site_info.UnitType = convertUnitToName(unit_i, length(runData.spikesByEventBinned{1}{chan_i}), 2);
      ULabel = convertUnitToName(unit_i, length(runData.spikesByEventBinned{1}{chan_i}), 1);

      % Save file in output folder
      save(fullfile(rasterDataPath, [runList{run_i}, ULabel]), 'raster_data', 'raster_labels', 'raster_site_info')
      
    end
  end
end



end
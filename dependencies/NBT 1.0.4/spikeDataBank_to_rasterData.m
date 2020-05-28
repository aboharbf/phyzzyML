function  spikeDataBank_to_rasterData(spikeDataBank, rasterDataPath, params)
% Turns the spikeDataBank struct of processBatchAnalysis into the format
% needed to process by Neural Decoding Toolbox. The output of this file can
% only be used as input into a phyzzy version of
% 'create_binned_data_from_raster_data'. After this, the pipelines should
% be the same.

% Creates file with the following 3 fields
% raster_data - 1 * unit cell array, each cell containing trial * bins.
% raster_labels - struct with each field contain 1*trial labels. Labels
% like 'stimulus' or 'position' go here.
% raster_site_info - file*1 labels/numbers - session_ID, channel, unit,
% site.

runList = fields(spikeDataBank);

if ~exist(rasterDataPath, 'dir')
  mkdir(rasterDataPath)
end

%runData.spikesByEventBinned{stim}{chan}{unit} --> needs to be turned to
%rasterDataUnit

for run_i = 1:length(runList)
  runData = spikeDataBank.(runList{run_i});
  
  % Initialize some variables which make the rest easier to gather.
  dataLength = abs(runData.start) + runData.end;
  vecLength = size(runData.spikesByEventBinned{1}{1}{1},2);
  if dataLength ~= vecLength
    startTime = floor((vecLength - dataLength)/2);
    endTime = dataLength+startTime-1;
  else
    startTime = 1;
    endTime = dataLength;
  end
  
  trialsPerStim = cellfun(@(x) size(x{1}{1},1), runData.spikesByEventBinned);
  trialsPerStimStart = cumsum([1; trialsPerStim(1:end-1)]);
  trialsPerStimEnd = cumsum(trialsPerStim);
  trialsPerUnit = sum(trialsPerStim);
  
  stimuliVec = cellfun(@(x) convertCharsToStrings(x), runData.eventIDs);
  tmp = arrayfun(@(x) repmat(stimuliVec(x), [trialsPerStim(x), 1]), 1:length(trialsPerStim), 'UniformOutput', 0)';
  stimuliVec = vertcat(tmp{:});
  

  tmp = arrayfun(@(x) repmat(runData.stimPresCount(x), [trialsPerStim(x), 1]), 1:length(trialsPerStim), 'UniformOutput', 0)';
  stimPresCountPerTrial = vertcat(tmp{:});
  
  % Construct raster_data and save file per unit.
  for chan_i = 1:length(runData.spikesByEventBinned{1})   
      for unit_i = 1:length(runData.spikesByEventBinned{1}{chan_i})
        
        % Cycle through the stim, store into raster_data.
        raster_data = zeros(trialsPerUnit, dataLength);
        for stim_i = 1:length(runData.spikesByEventBinned)
          raster_data(trialsPerStimStart(stim_i):trialsPerStimEnd(stim_i),:) = runData.spikesByEventBinned{stim_i}{chan_i}{unit_i}(:, startTime:endTime);
        end
        
        % Generate the raster_labels - labels that exist per trial.
        raster_labels.stimuli = stimuliVec';
        raster_labels.stimPresCount = stimPresCountPerTrial';
        % to include social vs non-social, need 
        % [plotMat, briefStimList, params]= plotIndex(stimuliList, plotIndParams)
        tmp = plotIndex(raster_labels.stimuli, params.plotIndParams);
        raster_labels.social = tmp;

        
        % Generate the raster_site_info
        raster_site_info.session_ID = runList{run_i};
        raster_site_info.alignment_event_time = abs(runData.start);
        raster_site_info.UnitType = convertUnitToName(unit_i, length(runData.spikesByEventBinned{1}{chan_i}), 2);
        raster_site_info.gridHole = num2str(runData.gridHoles{chan_i});
        raster_site_info.recordingDepth = runData.recDepth{chan_i};
        raster_site_info.runNum = runData.runNum;
        
        ULabel = convertUnitToName(unit_i, length(runData.spikesByEventBinned{1}{chan_i}), 1);

        % Save file in output folder
        save(fullfile(rasterDataPath, [runList{run_i}, ULabel]), 'raster_data', 'raster_labels', 'raster_site_info')

      end
  end
end



end
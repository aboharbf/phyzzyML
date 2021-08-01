function dataHighPrep(spikePathBank, batchAnalysisParams)
% Prepares a folder for processing data to prepare for data high.
% 1. Format raw spike trains into a Matlab struct, D. D(i) refers to the
% ith trial. The D(i).data field contains a matrix of size number of
% neurons x number of timepoints. (neurons * bins?)

saveDir = fullfile(batchAnalysisParams.outputDir, 'dataHigh');
dataStoreDir = fullfile(saveDir, 'preProc');


if ~exist(dataStoreDir, 'dir')
  mkdir(dataStoreDir);
end
    
paradigmList = unique(spikePathBank.paradigmName);
dateSubjClipped = extractBefore(spikePathBank.Row, '00');
nameVec = cellfun(@(x) x(10:end), dateSubjClipped, 'UniformOutput', false);
nameArray = unique(nameVec);

for m_i = 1:length(nameArray)
  for par_i = 1:length(paradigmList)
    
    % Pull the fraction of the table relevant to the paradigm and monkey
    pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
    mInd = contains(spikePathBank.Row, nameArray{m_i});
    spikePathBankParadigm = spikePathBank(pInd & mInd,:);
    
    % Make a directory to save the files.
    parMonkeyDir = fullfile(saveDir, nameArray{m_i}, paradigmList{par_i});
    
    if ~exist(parMonkeyDir, 'dir')
      mkdir(parMonkeyDir);
    end
    
    % grab the spikesByEvents and params.
    [spikesByEventBinnedPerRun, psthParamsPerRun] = spikePathLoad(spikePathBankParadigm, {'spikesByEventBinned', 'psthParams'}, batchAnalysisParams.spikePathLoadParams);
    
    % Identify unique events, and trials per each event, find the min.
    eventListPerRun = spikePathBankParadigm.stimuli;
    uniqueStim = unique([eventListPerRun{:}]);
    stimTrialPerRun = nan(length(uniqueStim), length(eventListPerRun));
    for run_i = 1:length(spikesByEventBinnedPerRun)
      for event_i = 1:length(spikesByEventBinnedPerRun{run_i})
        u_event_i = strcmp(uniqueStim, eventListPerRun{run_i}{event_i});
        stimTrialPerRun(u_event_i, run_i) = size(spikesByEventBinnedPerRun{run_i}{event_i}{1}{1}, 1);
      end
    end
    
    % Stimulus sets vary within Paradigms. Identify which runs have which
    % sets
    stimAccountedFor = false(1, length(eventListPerRun));
    stimSetInd = zeros(size(stimAccountedFor));
    run_i = 1;
    while ~all(stimAccountedFor)
      % see which stim exist in this set;
      stimAccountedFor = stimAccountedFor | ~isnan(stimTrialPerRun(run_i,:));
      % Add to set
      stimSetInd(~isnan(stimTrialPerRun(run_i,:))) = run_i;
      % Increment
      run_i = run_i + 1;
    end
    
    uniqueStimSets = unique(stimSetInd);
    trialsPerStim = min(stimTrialPerRun(:));
    %   eventListSubStim = eventListPerRunSubStim{1};
    spikeDataPerEvent = cell(size(eventListPerRun{1}, 1), length(uniqueStimSets));
    binCount = size(spikesByEventBinnedPerRun{1}{1}{1}{1}, 2);
    [spikeDataPerEvent{:}] = deal(nan(trialsPerStim, binCount, 1000));
    stimSets = [];
    
    for sub_i = 1:length(uniqueStimSets)
      % Find runs with this subset of stimuli
      subStimRunInd = stimSetInd == uniqueStimSets(sub_i);
      eventListPerRunSubStim = eventListPerRun(subStimRunInd);
      stimSets = [stimSets, eventListPerRunSubStim{1}];
      
      spikesByEventBinnedPerRunSubStim = spikesByEventBinnedPerRun(subStimRunInd);
      
      % create the dhStruct(trial_i).data(neurons*trials)
      for run_i = 1:length(eventListPerRunSubStim)
        for event_i = 1:length(eventListPerRunSubStim{run_i})
          % CONCERN: event tables might not line up across different runs.
          % They likely do though, from quick inspection. Won't work with
          % older data where stimuli aren't neatly broken into subsets.
          
          % create a 3D Matrix of trials*bins*MUA - later chop
          % along trials and break into right format for DataHigh.
          spikeData = spikesByEventBinnedPerRunSubStim{run_i}{event_i};
          
          for chan_i = 1:length(spikeData)
            % find the NaN
            TDind = find(isnan(spikeDataPerEvent{event_i, sub_i}(1, 1,:)), 1);
            
            % For now, just grab MUA.
            spikeDataPerEvent{event_i, sub_i}(:, :, TDind) = spikeData{chan_i}{end}(1:trialsPerStim,:);
          end
          
          
        end
      end
    end
    
    % Cut off the excess NaNs
    for ii = 1:size(spikeDataPerEvent, 1)
      for jj = 1:size(spikeDataPerEvent, 2)
        % Find the first NaN
        NaNInd = find(isnan(spikeDataPerEvent{ii, jj}(1, 1,:)), 1);
        
        % Save everything before it.
        spikeDataPerEvent{ii, jj} = spikeDataPerEvent{ii, jj}(:, :, 1:NaNInd-1);
        
      end
    end
    
    % create the dhStruct(trial_i).data(neurons*trials)
    % Example: D(1):
    % type: 'traj' or 'state' (1 x number of epochs)
    % epochStarts: [1 36 78]
    % epochColors: [0.7 0.7 0.7; 1 0 0; 0 1 0]  (number of epochs x RGB)
    % condition: 'left_reach' ('conditionName')
    % data: [7x145 array] (number of neurons x number of data points)
    params2Use = psthParamsPerRun{1}; % All variants of this paradigm should have the same params.
    
    % Data to include: Remove the ITI and Reward period
    dataStartInd = 0;
    dataEndInd = params2Use.psthPre + params2Use.psthImDur;
    
    % Shift the values accordingly
    padding = params2Use.movingWin(1)/2;
    fixEpoch = 1;
    stimOnsetEpoch = params2Use.psthPre;
    stimPresEpoch = stimOnsetEpoch + 500;
    rewardEpoch = params2Use.psthPre + params2Use.psthImDur;
    
    epochNames = {'Fixation', 'stimOnset', 'stimPres', 'reward'};
    %   epochStarts = [1, params2Use.ITI, params2Use.psthPre, params2Use.psthPre + 500, params2Use.psthPre + params2Use.psthImDur];
    epochStarts = [fixEpoch stimOnsetEpoch stimPresEpoch rewardEpoch] + padding;
    condType = 'traj';
    colors = [[0.0 0.0 1.0];[0.0 1.0 0.0];[0.0 0.5 0.0];[1.0 0.0 0.0]];
    
    for subSet_i = 1:size(spikeDataPerEvent, 2)
      
      % Initialize structure, reset count.
      dhStruct = struct();
      totalTrialCount = 1;
      
      for event_i = 1:size(spikeDataPerEvent,1)
        for trial_i = 1:trialsPerStim
          
          dhStruct(totalTrialCount).data = permute(squeeze(spikeDataPerEvent{event_i, subSet_i}(trial_i, :, :)), [2, 1]);
          dhStruct(totalTrialCount).type = condType;
          dhStruct(totalTrialCount).condition = stimSets{event_i, subSet_i};
          dhStruct(totalTrialCount).epochStarts = epochStarts;
          dhStruct(totalTrialCount).epochColors = colors;
          
          totalTrialCount = totalTrialCount + 1;
          
        end
      end
      
      % Generate file name
      dbStructFileName = sprintf('%s_%s_SS%d', nameArray{m_i}, paradigmList{par_i}, subSet_i);
      
      % Save this file into the folder
      save(fullfile(dataStoreDir, dbStructFileName), 'dhStruct')
      
      clear dhStruct
      
    end
    
  end
end
% 2. Call DataHigh(D,‘DimReduce’); in the Matlab command line.

end
function dataHighPrep(spikePathBank, batchAnalysisParams)
% Prepares a folder for processing data to prepare for data high.
% 1. Format raw spike trains into a Matlab struct, D. D(i) refers to the
% ith trial. The D(i).data field contains a matrix of size number of
% neurons x number of timepoints. (neurons * bins?)


saveDir = 'D:\DataAnalysis\batchAnalysis\dataHigh';

if ~exist(saveDir, 'dir')
  mkdir(saveDir);
end

paradigmList = unique(spikePathBank.paradigmName);

for par_i = 2:length(paradigmList)
  % Pull the fraction of the table relevant to the paradigm
  parRunInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  spikePathBankParadigm = spikePathBank(parRunInd,:);
  
  % Make a directory to save the files.
  parDir = fullfile(saveDir, paradigmList{par_i});
  
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
  trialsPerStim = min(min(stimTrialPerRun, [], 2));
%   eventListSubStim = eventListPerRunSubStim{1};
  spikeDataPerEvent = cell(size(eventListPerRun{1}, 1), length(uniqueStimSets));
  binCount = size(spikesByEventBinnedPerRun{1}{1}{1}{1}, 2);
  [spikeDataPerEvent{:}] = deal(nan(trialsPerStim, binCount, 1000));
  stimSets = [];
  
  for sub_i = 1:length(uniqueStimSets)
    % Find runs with this subset of stimuli
    subStimRunInd = stimSetInd == uniqueStimSets(sub_i);
    eventListPerRunSubStim = eventListPerRun(subStimRunInd);
    stimSets = [stimSets, eventListPerRunSubStim{sub_i}];
    
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
  dataStartInd = params2Use.ITI;
  dataEndInd = params2Use.psthPre + params2Use.psthImDur;
  
  % Shift the values accordingly
  fixEpoch = 1;
  stimOnsetEpoch = params2Use.psthPre - params2Use.ITI;
  stimPresEpoch = stimOnsetEpoch + 500;
  
  
  epochNames = {'Fixation', 'stimOnset', 'stimPres'};
%   epochStarts = [1, params2Use.ITI, params2Use.psthPre, params2Use.psthPre + 500, params2Use.psthPre + params2Use.psthImDur];
  epochStarts = [fixEpoch stimOnsetEpoch stimPresEpoch];
  condType = 'traj';
  colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0]];
  
  for subSet_i = 1:size(spikeDataPerEvent, 2)
    
    % Initialize structure, reset count.
    dhStruct = struct();
    totalTrialCount = 1;
    
    for event_i = 1:size(spikeDataPerEvent, 1)
      for trial_i = 1:trialsPerStim
        
        dhStruct(totalTrialCount).data = squeeze(spikeDataPerEvent{event_i, subSet_i}(trial_i, 150 + params2Use.ITI:end-650, 1:32))';
        dhStruct(totalTrialCount).type = condType;
        dhStruct(totalTrialCount).condition = stimSets{event_i, subSet_i};
        dhStruct(totalTrialCount).epochStarts = epochStarts;
        dhStruct(totalTrialCount).epochColors = colors;
        
        
        totalTrialCount = totalTrialCount + 1;
        
      end
    end
    
    % Generate file name
    dbStructFileName = sprintf('%s_%d', paradigmList{par_i}, subSet_i);
    
    % Save this file into the folder
    save(fullfile(saveDir, dbStructFileName), 'dhStruct')
    
  end
  

 
  
  
end

% 2. Call DataHigh(D,‘DimReduce’); in the Matlab command line.

end
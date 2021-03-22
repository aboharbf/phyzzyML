function spikeEyeOverlay(spikePathBank, params, figStruct)
% This function generates a video of each of the stimuli w/ activity
% overlaid in a heatmap.
% Inputs:
% - spikePathBank: with the fields for events, analyzedDir.
% - batchAnalysisParams: with field for spikeLoadParams.
% - figStruct

params.spikePathLoadParams.fileVars = [params.spikePathLoadParams.fileVars, 'spikeEyeData'];
params.spikePathLoadParams.fileInds = [params.spikePathLoadParams.fileInds, 4];

[paradigmList, ~, allParadigmRuns] = unique(spikePathBank.paradigmName);

% Define variables to load for every paradigm here
if params.meanPSTHParams.normalize
  variables2Extract = {'spikeEyeData', 'psthParams', 'eyeDataStruct', 'psthByImageProc'};
else
  variables2Extract = {'spikeEyeData', 'psthParams', 'eyeDataStruct', 'psthByImage'};
end

for par_i = 1:length(paradigmList)
  
  % Generate list of unique stimuli present in this paradigm.
  spikePathBankParadigm = spikePathBank(allParadigmRuns == par_i, :);
  uniqueStim = unique(vertcat(spikePathBankParadigm.stimuli{:}));
    
  % Load the necessary variables.
  [spikeEyeDataPerRun, psthParamsPerRun, eyeDataStruct, psthByImagePerRun] = spikePathLoad(spikePathBankParadigm, variables2Extract, params.spikePathLoadParams);
  
  % Load the first example of this paradigm, and extract parameters needed.
  psthPre = psthParamsPerRun{1}.psthPre;
  psthImDur = psthParamsPerRun{1}.psthImDur;
  psthPost = psthParamsPerRun{1}.psthPost;
  padSize = psthParamsPerRun{1}.movingWin(1)/2;
  traceLength = psthPre + psthImDur + psthPost + (padSize * 2) + 1;
  
  % Initialize things for each stim
  [frameStarts, frameEnds] = deal(cell(size(uniqueStim)));
  runStimMat = false(length(uniqueStim), length(eyeDataStruct));
  for stim_i = 1:length(uniqueStim)
    
    % ID the runs w/ this stim
    runStimMat(stim_i, :) = cellfun(@(x) contains(uniqueStim{stim_i}, x), spikePathBankParadigm.stimuli);
    
    % Grab the first run with the stimulus, retrieve frame data.
    for run_i = find(runStimMat(stim_i, :), 1)
      
      % Use previously calculated 'bins per frame' variables
      eventInd = strcmp(eyeDataStruct{run_i}.attendedObjData.eventIDs, uniqueStim{stim_i});
      if any(eventInd)
        frameStarts{stim_i} = eyeDataStruct{run_i}.attendedObjData.frameStartInd{eventInd};
        frameEnds{stim_i} = eyeDataStruct{run_i}.attendedObjData.frameEndInd{eventInd};
      end
      
    end
    
    
  end
  
  % Gather Spike Information, based on runStimMat entries
  stimPSTH = cell(length(uniqueStim), 3);
  [stimPSTH{:}] = deal(nan(1e5, traceLength));
  
  % Cycle through stim, gathering the traces
  for stim_i = 1:length(uniqueStim)
    runs2Check = find(runStimMat(stim_i,:));
    for run_i = runs2Check
      
      % Find the data to extract
      eventInd = strcmp(eyeDataStruct{run_i}.attendedObjData.eventIDs, uniqueStim{stim_i});
      
      % Extract and sort
      runAct = spikesByEventBinnedPerRun{run_i}{eventInd};
      for ch_i = 1:length(runAct)
        uCount = length(runAct{ch_i});
        for u_i = 1:uCount
          
          % Extract activity trace
          activityTrace = runAct{ch_i}{u_i};
          
          % Don't add if it doesn't meet a threshold (if set)
          if params.spikeEyeOverlay.threshold && sum(activityTrace)/(length(activityTrace)/1000) < params.spikeEyeOverlay.thresholdHz
              continue
          end
          % Identify Group of activity
          if u_i == 1
            group_i = 1;
          elseif u_i == uCount
            group_i = 3;
          else
            group_i = 2;
          end
                    
          % Find NaN start
          nanStart = find(isnan(stimPSTH{stim_i, group_i}(:,1)), 1);
                   
          % Add to Larger Stack.
          stimPSTH{stim_i, group_i}(nanStart, :) = activityTrace;
        end
      end
      
      
    end

  end
  
  % For every stim, cycle through it, binning the spikes for each frame,
  % and saving this in a new structure
  
  % Once all the spike rasters have been collected 
      % Load the stimuli
    %     vidPath = fullfile(params.spikeEyeOverlay.stimDir, '**', uniqueStim{stim_i});
    %     vidObjPath = dir(vidPath);
    %     vidObj = VideoReader(fullfile(vidObjPath.folder, vidObjPath.name));    

  
  
end

end

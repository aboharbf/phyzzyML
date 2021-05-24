function [stimPSTH, dataType] = compilePsthPerStimArray(stimuliByRun, psthByImagePerRun, allStimuliVec)
% Go across run data, collecting the PSTHes for each stimuli.
% Output is stimPSTH{stimuli}{group}{dataType};


removeExcessNaNs = @(cellWNans) cellWNans(~isnan(cellWNans(:,1)), :);
dataType = {'PSTH', 'Run Ind'};

stimPSTH = cell(length(allStimuliVec), 3, length(dataType));
[stimPSTH{:, :, 1}] = deal(nan(1000, size(psthByImagePerRun{1}{1}{1}, 2)));
[stimPSTH{:, :, 2}] = deal(nan(1000, 1));

for stim_i = 1:length(allStimuliVec)
  % ID which runs have this stimulus
  runs2Check = find(cellfun(@(x) any(contains(x, allStimuliVec{stim_i})), stimuliByRun));
  
  for run_i = runs2Check'
    stimActInd = strcmp(stimuliByRun{run_i}, allStimuliVec(stim_i));
    
    for chan_i = 1:length(psthByImagePerRun{run_i})
      unitCount = length(psthByImagePerRun{run_i}{chan_i});
      for unit_i = 1:unitCount
        % ID The group
        if unit_i == 1
          group_i = 1;
        elseif unit_i == unitCount
          group_i = 3;
        else
          group_i = 2;
        end
        
        % Go through and stack it correctly.
        firstInd = find(isnan(stimPSTH{stim_i, group_i, strcmp(dataType, 'PSTH')}(:,1)), 1);
        act2Add = psthByImagePerRun{run_i}{chan_i}{unit_i}(stimActInd, :);
        
        stimPSTH{stim_i, group_i, strcmp(dataType, 'PSTH')}(firstInd:firstInd+size(act2Add,1)-1, :) = act2Add;
        stimPSTH{stim_i, group_i, strcmp(dataType, 'Run Ind')}(firstInd:firstInd+size(act2Add,1)-1, :) = run_i;
      end
    end
  end
  
end

% Post processing - remove rows which are just NaN
stimPSTH = cellfun(removeExcessNaNs, stimPSTH, 'UniformOutput', false);

end
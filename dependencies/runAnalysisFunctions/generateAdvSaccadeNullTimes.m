function [stimSaccadeArrays, stimSaccadeNullArrays] = generateAdvSaccadeNullTimes(eyeBehStatsByStim, psthParams)
% Goal of this function is to be fed in saccade times and to produce null
% time scrambles.
% Input 
% - eyeBehStatsByStim

preSaccWin = 300;
postSaccWin = 100;

psthPre = psthParams.psthPre;
psthImDur = psthParams.psthImDur;
psthPost = psthParams.psthPost;

[stimSaccadeArrays, stimSaccadeNullArrays] = deal(cell(size(eyeBehStatsByStim)));
trialBins = length(-psthPre:(psthImDur+psthPost));
for stim_i = 1:length(eyeBehStatsByStim)
  % Create a representation of all the saccades for particular trials of a
  % stimulus.
  [stimSaccadeTimes, stimSaccadeNullArrays{stim_i}] = deal(zeros(length(eyeBehStatsByStim{stim_i}), trialBins));
  for trial_i  = 1:length(eyeBehStatsByStim{stim_i})
    saccInd = eyeBehStatsByStim{stim_i}{trial_i}.saccadetimes + psthPre;
    if ~isempty(saccInd)
      stimSaccadeTimes(trial_i, saccInd(1,:)) = 1;
    end
  end
  stimSaccadeArrays{stim_i} = stimSaccadeTimes;
end

for stim_i = 1:length(stimSaccadeArrays)
  stimSaccadeTimes = stimSaccadeArrays{stim_i};
  trialIndLog = logical(eye(size(stimSaccadeTimes,1)));
  
  for trial_i = trialIndLog
    if any(stimSaccadeTimes(trial_i,:))
      % Check the area surrounding the saccade mark.
      saccBins = find(stimSaccadeTimes(trial_i, :));
      otherTrials = stimSaccadeTimes(~trial_i,:);
      
      for sacc_i = 1:length(saccBins)
        % Grab the saccade time
        saccWin = max(saccBins(sacc_i) - preSaccWin, 1):min(saccBins(sacc_i)+postSaccWin, size(otherTrials,2));
        saccWinOther = otherTrials(:, saccWin);
        
        % Identify other trials without saccades at the same spot
        nonSaccTrials = ~any(saccWinOther,2);
        
        if any(nonSaccTrials)
          % Add those spots to the non-saccade set
          stimSaccadeNullArrays{stim_i}(nonSaccTrials, saccBins(sacc_i)) = 1;
          
        else
          % Go through the other stimuli, until you find one where the same
          % window is unoccupied.
          nullSaccSearch = true;
          newStim_i = stim_i + 1;
          
          % Check until you run out of stim or you find stim.
          while nullSaccSearch && (newStim_i <= length(stimSaccadeArrays))
            nonSaccTrials = ~any(stimSaccadeArrays{newStim_i}(:,saccWin),2);
            
            if any(nonSaccTrials)
              % Once found, add them to the array, move on
              nullSaccSearch = false;
              stimSaccadeNullArrays{newStim_i}(nonSaccTrials, saccBins(sacc_i)) = 1;
            else
              % Otherwise, increment, check the next stim.
              newStim_i = newStim_i + 1;
            end
            
          end
        end
      end
    end
  end
  
end

end
      
      
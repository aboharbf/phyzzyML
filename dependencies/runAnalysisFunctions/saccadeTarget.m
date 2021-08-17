function eyeDataStruct = saccadeTarget(eyeDataStruct)
% For each saccade, label what the target of it is.
% Inputs:
% - eyeDataStruct.saccadeByStim from eyeStatsAnalysis.
% - .eyeBehStatsByStim - from calcEyeObjectTrace.
% Outputs:
% - eyeBehStatsByStim, with each trail struct having the new field
% 'objAttended', with an item for saccades which take place during stim.

% Note
% - Given the nature of CluterFix (identifying saccades as stretches of
% transitions, rather than a single one), it may identify the biggest jump
% in the saccade as being in the middle, and if it is followed by a  jump
% back, may look very small in distance. I've added the max distance jump
% info to the output of eyeStatsAnalysis to get around this.

eyeBehStatsByStim = eyeDataStruct.eyeBehStatsByStim;

for stim_i = 1:length(eyeBehStatsByStim)
  
  % Create a reference for the stimulus's frames
  frameTimes = [eyeDataStruct.attendedObjData.frameStartInd{stim_i}' ,eyeDataStruct.attendedObjData.frameEndInd{stim_i}'];
  stimEnd = max(frameTimes(:));
  attendedObjPerTrial = eyeDataStruct.attendedObjData.attendedObjVect{stim_i};
  
  % Identify all ths saccades
  for trial_i  = 1:length(eyeBehStatsByStim{stim_i})
    
    % Filter for saccades during stimulus presentation
    saccadetimes = eyeBehStatsByStim{stim_i}{trial_i}.saccadetimes;
    if ~isempty(saccadetimes)
      objAttended = cell(1, size(saccadetimes,2));
      saccKeep = saccadetimes(1,:) > 0 & saccadetimes(2,:) < stimEnd;
      
      % Look at the saccade end times, and see which frame that occurs
      % during, and use that to find the object.
      for sacc_i = 1:size(saccadetimes, 2)
        % If the saccade happens during the stim, ends during the stim.
        if saccKeep(sacc_i)
          
          saccadeEndTime = saccadetimes(2, sacc_i);
          
          frameInd = frameTimes(:,1) <= saccadeEndTime & frameTimes(:,2) >= saccadeEndTime;
          objAttended(sacc_i) = attendedObjPerTrial(trial_i, frameInd);
        end
        
      end
      
      % Save with the saccade data
      eyeBehStatsByStim{stim_i}{trial_i}.objAttended = objAttended;
    else
      eyeBehStatsByStim{stim_i}{trial_i}.objAttended = [];
    end
    
  end
  
end

eyeDataStruct.eyeBehStatsByStim = eyeBehStatsByStim;



end
function [saccadeMatArray, saccadeMatLabels] = saccadePerEvent(eyeBehStatsByStim, taskData, trialIDsByEvent, psthParams, params)
% A function which plots saccade frequency around major events.
% - fixationDotOnset
% - stimulusOnset
% - stimulusDuration for each stim.
% - reward

% Pull needed variables
taskEventStartTimes = taskData.taskEventStartTimes;
fixOnsetShift = round(taskData.taskEventStartTimes - taskData.taskEventFixDur);

rewardOnTimesPerTrial = round(taskData.juiceOnTimes - taskEventStartTimes);
eventIDs = taskData.taskEventList;
eventCount = length(eventIDs);

preEvent = params.preEventTime;
postEvent = params.postEventTime;
preStim = params.preStimTime;
stimDur = psthParams.psthImDur;
postStim = params.postStimTime;
smoothingWidth = params.smoothingWidth;

% Generate smoothing filter
filterPoints = -3*smoothingWidth:3*smoothingWidth;
smoothingFilter = exp(-1*filterPoints.^2/(2*smoothingWidth^2));
smoothingFilter = smoothingFilter/sum(smoothingFilter);

% Convert the saccadeByStim into just being about saccades (get rid of
% other values) label saccade times

[stimPresSacc, fixOnsetSacc, rewardSacc] = deal(cell(size(trialIDsByEvent)));
for event_i = 1:length(trialIDsByEvent)
  % Collect saccade times for this stimulus
  saccadeDataEvent = vertcat(eyeBehStatsByStim{event_i}{:});
  saccadeIndices = {saccadeDataEvent.saccadetimes};

  stimPresSacc{event_i} = indicies2Matrix(saccadeIndices, 0, preStim, stimDur, postStim);                                             % Times already aligned to stimulus onset.
  fixOnsetSacc{event_i} = indicies2Matrix(saccadeIndices, fixOnsetShift(trialIDsByEvent{event_i}), preEvent, 0, postEvent);           % Fix times are relative to stimOnset, so offset is postive.
  rewardSacc{event_i} = indicies2Matrix(saccadeIndices, -rewardOnTimesPerTrial(trialIDsByEvent{event_i}), preEvent, 0, postEvent);    % Reward times are post stimOnset, so offset is negative
      
end

% Set this into a a larger vector
saccadeMatArray = cell(eventCount + 3,1);
saccadeMatArray(1:eventCount) = stimPresSacc;

% Stim onset
stimOnsetSacc = cellfun(@(stimMat) stimMat(:, 1 : preStim + postEvent), stimPresSacc, 'UniformOutput', false);
saccadeMatArray{eventCount+1} = vertcat(stimOnsetSacc{:});

% Fix and reward
saccadeMatArray{eventCount+2} = vertcat(fixOnsetSacc{:});
saccadeMatArray{eventCount+3} = vertcat(rewardSacc{:});

% Create simplified traces - take mean and convolve
% saccadeMatArray = cellfun(@(x) conv(mean(x, 1), smoothingFilter, 'same'), saccadeMatArray, 'UniformOutput', false);

saccadeMatLabels = [eventIDs; 'Stimulus Onset'; 'Fixation Onset'; 'Reward Delivery'];

% figure()
% subplot(2,1,1)
% plot(vertcat(saccadeMatArray{1:eventCount})')
% legend(saccadeMatLabels(1:eventCount));
% 
% subplot(2,1,2)
% plot(saccadeMatArray{eventCount+1})
% hold on
% plot(saccadeMatArray{eventCount+2})
% plot(saccadeMatArray{eventCount+3})
% legend(saccadeMatLabels(eventCount+1:eventCount+3));

end
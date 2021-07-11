function [eyeDataStruct, eyeBehStatsByStim, eyeInByEventSmooth] = eyeStatsAnalysis(analogInByEvent, psthParams, eyeStatsParams, taskData, eyeDataStruct, figStruct)
% Function uses ClusterFix to generate a saccade map.

% Outputs
% - eyeDataStruct: Also an Input, has the saccadeByStim structure added,
% which create a color map per stimulus of eye behavior, as indicated in
% the eventInds and eventNames variables below.
% - eyeBehStatsByStim - a {stim}{trial} struct which details outputs of
% ClusterFix per trial.
% - eyeInByEventSmooth - the smoothed eye signal extracted from
% analogInByEvent and processed by ClusterFix, and returned to normal
% format.

% Relevant Variables for processing
blinkVals = taskData.eyeCal.blinkVals;
PixelsPerDegree = taskData.eyeCal.PixelsPerDegree;

% Code for the output saccadeByStim Image.
eventInds = [1, 2, 3, 4];      
eventNames = {'Fix', 'PreSaccade', 'Saccade','Blink'};
preSaccPeriod = 200;

% A few switches.
filterSaccades = 1;         % Performs a filter on saccades to exclude microsaccades, defined with respect to duration
saccadeDurThres = 35;       % Additional saccade filtering outside of clusterFix, removes saccades which don't last this long or more.
saccadeDistThrs = 1.5;      % As above, removes saccades with mean distances of less than this. 
plotPaths = 0;              % Code below which visualizes things trial by trial.

% visualized trial by trial trace parameters
saccadeColors = 'rbgcm';
psthPre = psthParams.psthPre;
psthImDur = psthParams.psthImDur;
lfpPaddedBy = (psthParams.movingWin(1)/2)+1;

% Remove Padding, as this messes with indexing and the code signal pads itself.
stimStartInd = lfpPaddedBy;
stimEndInd = size(analogInByEvent{1},4) - lfpPaddedBy;

% Reshape the analogInByEvent into eye data. Can not use cellfun w/ squeeze
% due to single trial stimuli paradigms.
eyeInByEventAll = cell(size(analogInByEvent));
for event_i = 1:length(analogInByEvent)
  eyeTmp = analogInByEvent{event_i}(:,1:2,:,stimStartInd:stimEndInd); %Reshapes analogInByEvent into eye signal.
  eyeSize = size(eyeTmp);
  eyeInByEventAll{event_i} = reshape(eyeTmp, [eyeSize(2:4)]); 
end

[fixStats, eyeInByEventTrialFiltered, blinkTimes] = deal(cell(size(analogInByEvent)));

% Remove Blinks, save their times in the same format clusterFix does.
blinkThreshold = round(mean(blinkVals)*0.8); % Removes top 20% of values.
for stim_i = 1:length(eyeInByEventAll)
  trialCount = size(eyeInByEventAll{stim_i}, 2);
  blinkTimes{stim_i} = cell(trialCount,1);
  for trial_i = 1:trialCount
    for eye_i = 1:size(eyeInByEventAll{stim_i}, 1)
      trace = squeeze(eyeInByEventAll{stim_i}(eye_i, trial_i, :));
      val2Replace = trace > blinkThreshold;
      if sum(val2Replace) > 0
        % Remove Blinks
        blinkInds = find(val2Replace);
        blinkTimesTmp = BehavioralIndex(blinkInds'); % Function from ClusterFix.
        for blink_i = 1:size(blinkTimesTmp,2)
          % find a start value - mean of 10 prior values
          startInds = (blinkTimesTmp(1,blink_i)-12 : blinkTimesTmp(1,blink_i) - 2);
          if any(startInds <= 0)
            startVal = trace(blinkTimesTmp(2,blink_i)+5);
            trace(blinkTimesTmp(1,blink_i):blinkTimesTmp(2,blink_i)+5) = deal(startVal);
          else
            % Replace values in the trace - easy method, may need to
            % complicate.
            startVal = mean(trace(startInds));
            trace(blinkTimesTmp(1,blink_i)-2:blinkTimesTmp(2,blink_i)) = deal(startVal);
          end
          
        end
        
        if eye_i == 1
          blinkTimes{stim_i}{trial_i} = blinkTimesTmp;
        end
        eyeInByEventAll{stim_i}(eye_i, trial_i, :) = trace;
      else
        blinkTimes{stim_i}{trial_i} = [];
      end        
    end
  end
end

% Turn each cell into the properly formated trial structure for clusterFix.
perTrialBreak = @(x) cellfun(@(y) squeeze(y), mat2cell(x, 2, ones(size(x, 2), 1), size(x, 3)), 'UniformOutput', 0);
eyeInByEventTrial = cellfun(perTrialBreak, eyeInByEventAll, 'UniformOutput', 0);

for stim_i = 1:length(eyeInByEventTrial)
  [fixStats{stim_i}, eyeInByEventTrialFiltered{stim_i}] = ClusterFix(eyeInByEventTrial{stim_i}, 1/1000);
end

% Filter saccades, make output image
% Convert the output of ClusterFix to an image which can be used behind
% Rasters.
[eyeBehImgByStim] = deal(cell(length(fixStats),1));
for stim_i = 1:length(fixStats)
  stimEye = fixStats{stim_i};
  [eyeBehImg] = deal(ones(length(stimEye), (stimEndInd-stimStartInd)+1));
  for trial_i = 1:length(stimEye)
    if ~isempty(stimEye{trial_i}.fixations)
      
      % Extract Fixations and mark their times
      fixationTimes = stimEye{trial_i}.fixationtimes;
      for fix_i = 1:size(fixationTimes,2)
        eyeBehImg(trial_i, fixationTimes(1, fix_i):fixationTimes(2, fix_i)) = deal(eventInds(strcmp(eventNames, 'Fix')));
      end
    end
    
    if ~isempty(stimEye{trial_i}.saccadetimes)
      
      % Extract Saccades and mark their times
      saccadeTimes = stimEye{trial_i}.saccadetimes;
      saccadeMaxDists = stimEye{trial_i}.SaccadeClusterValues(:,7)';
      
      % Filter if desired
      if filterSaccades
        % Filter based on duration (filterThres)
        sacDurTmp = saccadeTimes(2,:) - saccadeTimes(1,:);
        sacKeepInd = sacDurTmp <= saccadeDurThres;
        
        % Filter based on distance
        sacKeepInd2 = saccadeMaxDists < saccadeDistThrs;
        
        % Combine
        sacKeepInd = ~logical(sacKeepInd + sacKeepInd2);
        
        % Remove saccades and overwrite clusterFix output.
        [saccadeTimes, fixStats{stim_i}{trial_i}.saccadetimes] = deal(saccadeTimes(:, sacKeepInd));
        fixStats{stim_i}{trial_i}.SaccadeClusterValues = fixStats{stim_i}{trial_i}.SaccadeClusterValues(sacKeepInd,:);
      end
      
      % Count in reverse to avoid overwriting a saccade w/ pre saccade
      % period of a subsequent event.
      for sac_i = size(saccadeTimes,2):-1:1
        % Label Saccade times in the image
        eyeBehImg(trial_i, saccadeTimes(1, sac_i):saccadeTimes(2, sac_i)) = deal(eventInds(strcmp(eventNames, 'Saccade')));
        
        % If possible, label the period prior to the saccade
        preSaccWin = saccadeTimes(1, sac_i) - preSaccPeriod :saccadeTimes(1, sac_i);
        preSaccWin = preSaccWin(preSaccWin > 0);
        eyeBehImg(trial_i, preSaccWin) = deal(eventInds(strcmp(eventNames, 'PreSaccade')));
      end
      
      % 0s represents saccades which don't get over the threshold, for now
      % call them fixations.
      eyeBehImg(eyeBehImg == 0) = 1;
      
      % Add in Blinks
      blinksTrial = blinkTimes{stim_i}{trial_i};
      for blink_i = 1:size(blinksTrial,2)
        eyeBehImg(trial_i, blinksTrial(1, blink_i):blinksTrial(2, blink_i)) = deal(eventInds(strcmp(eventNames, 'Blink')));
      end
    else
      %disp('No saccades detected')
    end

  end
  eyeBehImgByStim{stim_i} = eyeBehImg;
end

if plotPaths
  eventIDs = eyeStatsParams.eventIDs;
  taskEventIDs = taskData.taskEventList;
  frameMotionInd = cellfun(@(x) find(strcmp(taskEventIDs, x)), eventIDs, 'UniformOutput', false);
  frameMotionInd = [frameMotionInd{:}];
  
  % Generate Images that illustrate Saccade Path
  for stim_i = 1:length(fixStats)
    % Find the right frame
    lowerLeft = -([taskData.frameMotionData(frameMotionInd(stim_i)).width/2 taskData.frameMotionData(frameMotionInd(stim_i)).height/2] ./ abs(taskData.eyeCal.PixelsPerDegree));
    
    h = figure('Units','normalized','Position', [0 0 1 0.9], 'Name', sprintf('Trial eye Traces - %s', eyeStatsParams.eventIDs{stim_i}));
    plotHandles = gobjects(length(fixStats{stim_i}),1);
    objTbGrp = uitabgroup('Parent', h);
    for trial_i = 1:length(fixStats{stim_i})
      objTab = uitab('Parent', objTbGrp, 'Title', sprintf('Trial %d', trial_i));
      plotHandles(trial_i) = axes('parent', objTab);
      xy = fixStats{stim_i}{trial_i}.XY;
      rectangle('Position', [[lowerLeft] abs([lowerLeft*2])], 'EdgeColor', 'b', 'LineWidth', 4);
      fixations = fixStats{stim_i}{trial_i}.fixations;
      fixationtimes = fixStats{stim_i}{trial_i}.fixationtimes;
      saccadetimes = fixStats{stim_i}{trial_i}.saccadetimes;
      saccMeanDist = round(fixStats{stim_i}{trial_i}.SaccadeClusterValues(:,3), 3);
      saccMaxDist = round(fixStats{stim_i}{trial_i}.SaccadeClusterValues(:,7), 3);
      
      hold on
      % Plot the entire trace in black for the sake of microsaccades which
      % were removed.
      a = plot(xy(1,:), xy(2,:), 'k');
      
      % Plot fixation means
      for ii = 1:size(fixationtimes, 2)
        plot(fixations(1,ii),fixations(2,ii),'y*'); % plot mean fixation location, yellow
      end
      
      % Plot saccades
      saccHandles = gobjects(size(saccadetimes, 2),1);
      saccLegend = cell(size(saccadetimes, 2),1);
      for ii = 1:size(saccadetimes,2)
        saccHandles(ii) = plot(xy(1,saccadetimes(1,ii):saccadetimes(2,ii)),...
          xy(2, saccadetimes(1,ii): saccadetimes(2,ii)),...
          saccadeColors(mod(ii - 1, length(saccadeColors)) + 1)); %   Each saccade gets a distinct color.
        saccLegend{ii} = sprintf('%s, %s', num2str(saccMaxDist(ii)), num2str(saccMeanDist(ii)));
      end
      
      %Plot the start and end
      text(xy(1,psthPre), xy(2,psthPre), 'S', 'Fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor','white')
      text(xy(1,psthPre+psthImDur), xy(2,psthPre+psthImDur), 'E', 'Fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor','white', 'Clipping', 'on')
      
      legend([a; saccHandles], [{'Fixations'}; saccLegend], 'location', 'NorthEastOutside')
      title(sprintf('Trial %d', trial_i));
      
    end
    % Link all the subplots
    linkprop(plotHandles, {'XLim', 'YLim', 'Position'});
    plotHandles(1).XLim = [-12 12];
    plotHandles(1).YLim = [-7 7];
    
    % Save the figure - Untested, just added to bring in line w/ other
    % code.
    saveFigure(eyeStatsParams.outDir, ['saccadeTraces_' eyeStatsParams.eventLabels{stim_i}], figData, figStruct, figStruct.figTag);
  end
  
end

% Reshape smoothed eye signal from clusterFix back into the appropriate
% form for the rest of phyzzy.
eyeInByEventSmooth = cell(size(eyeInByEventAll));
for stim_i = 1:length(eyeInByEventTrialFiltered)
  eyeSigTmp = eyeInByEventTrialFiltered{stim_i}';
  eyeSigReshaped = zeros(size(eyeSigTmp{1}, 1), length(eyeSigTmp), size(eyeSigTmp{1}, 2));
  for trial_i = 1:length(eyeSigTmp)
    eyeSigReshaped(1, trial_i, :) = eyeSigTmp{trial_i}(1,:);
    eyeSigReshaped(2, trial_i, :) = eyeSigTmp{trial_i}(2,:);
  end
  eyeInByEventSmooth{stim_i} = eyeSigReshaped;
end


% Shift saccade times in each stimuli, trial to put it with reference to
% stimulus onset, and add blinks.

for stim_i = 1:length(fixStats)
  for trial_i = 1:length(fixStats{stim_i})
    fixStats{stim_i}{trial_i}.fixationtimes = fixStats{stim_i}{trial_i}.fixationtimes - psthPre;
    fixStats{stim_i}{trial_i}.saccadetimes = fixStats{stim_i}{trial_i}.saccadetimes - psthPre;
    fixStats{stim_i}{trial_i}.blinktimes = blinkTimes{stim_i}{trial_i} - psthPre;
  end
end

eyeBehStatsByStim = fixStats;
eyeDataStruct.saccadeByStim = eyeBehImgByStim;

% Determine Saccade angles.
% padding = 20;
for stim_i = 1:length(eyeBehStatsByStim)
  for trial_i = 1:length(eyeBehStatsByStim{stim_i})
    % Pull and adjust saccade times
    saccadeTimes = eyeBehStatsByStim{stim_i}{trial_i}.saccadetimes + psthPre;
        
    [saccadeDir, saccadeAngle] = deal(zeros(1, size(saccadeTimes,2)));
    % Find the times just before it in the trace.
    trialEye = squeeze(eyeInByEventSmooth{stim_i}(:,trial_i, :));
    
    for sac_i = 1:size(saccadeTimes,2)
      saccadeEyePadded = trialEye(:, saccadeTimes(1, sac_i):saccadeTimes(2, sac_i));
      
      % Find the starting and stopping spots for the saccade
      startSpot = mean(saccadeEyePadded(:,1:3), 2);
      stopSpot = mean(saccadeEyePadded(:,end-3:end), 2);
      
      cardinalRanges = [0, 45, 90, 135, 180, 225, 270, 315, 360];
      delta_x = stopSpot(1) - startSpot(1);
      delta_y = stopSpot(2) - startSpot(2);
      angle = 180/pi * atan2(delta_y, delta_x); % Degrees initially.
      
      if angle < 0
        angle = 360 + angle;
      end
      
      % Find the angle for the saccade
      
      % Transform angles into Cardinal directions, and identify which each
      % belongs to.
      
      angleDiff = abs(angle - cardinalRanges);
      [~, dirInd] = min(angleDiff, [], 2);
      
      if 1 % Convert to Radians.
        angle = angle /(180/pi);
      end
      
      % Store
      saccadeDir(sac_i) = dirInd;
      saccadeAngle(sac_i) = angle;
      
    end
    
    % Store in output
    eyeBehStatsByStim{stim_i}{trial_i}.saccadeDirection = saccadeDir;
    eyeBehStatsByStim{stim_i}{trial_i}.saccadeAngle = saccadeAngle;
    
  end
end


end

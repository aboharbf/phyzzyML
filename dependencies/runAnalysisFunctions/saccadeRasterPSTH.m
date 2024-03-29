function saccadeRasterPSTH(eyeDataStruct, spikesByChannel, onsetsByEvent, rasterSaccParams, selTable, psthParams, figStruct)
% subEventAnalysis
% Description - looks through spikesByEvent, calculates PSTH for activity
% aligned to a specific event as well as a null distribution from the
% remaining stimuli. Produces statistical tests for the traces.
% Parameters
% eyeBehStatsbyStim - a structure generated by earlier eye analysis
% functions. {stim}{trial}
% spikesByChannel - a relatively raw stream of spikes, produced in
% processRun. time stamps are used to pull spikes from it.
% subEventParams - has path to eventData, stimDir, and variables used to
% plot PSTHes.

nonStimSaccadesOnly = true;
saccadeSplit = false;           % Split subEvent in the stimulus depending on whether they were followed by a saccade.
trialTemplate = [ones(1, psthParams.psthPre), ones(1, psthParams.psthImDur) * 2, ones(1, psthParams.psthPost) * 3];

eyeBehStatsByStim = eyeDataStruct.eyeBehStatsByStim;

% Highlight cells to which get to be plotted
saccadeSelCells = (selTable.subSel_saccades_pVal <= 0.05) | (selTable.subSel_saccadesNonStim_pVal  <= 0.05)...
  | (selTable.subSel_pre_saccades_pVal  <= 0.05) | selTable.saccDir_nonStim_pVal  <= 0.05;

% Highlight cells to which get to be plotted
% saccadeSelCells = (selTable.subSel_saccades_pVal <= 0.05) | (selTable.subSel_saccadesNonStim_pVal  <= 0.05)...
%   | (selTable.subSel_pre_saccades_pVal  <= 0.05) | (selTable.saccDir_all_pVal  <= 0.05);

% Generate/Organize Event times for eye movements here, concatonate onto
% the 'onsetsByEvent' cue.
trialsPerStim = cellfun('length', eyeBehStatsByStim);
circleInRadians = linspace(0, 2*pi, 9);
circleInRadians(end) = 0;
  
% Turn both the saccade times into rasters. In the same function, generate
% null rasters by sampling other trials of the same stimulus for
% non-saccade trials.
[stimSaccadeArrays, ~] = generateAdvSaccadeNullTimes(eyeBehStatsByStim, psthParams);

[blinkTimes, saccadeTimes, saccadeLabel, saccadeDir, saccadeImg] = deal([]);                          % Initialize vectors

for stim_i = 1:length(eyeBehStatsByStim)
  for trial_i = 1:trialsPerStim(stim_i)
    % Find Absolute eye event start times.
    eyeTrial = eyeBehStatsByStim{stim_i}{trial_i};
    
    % Grab blink times, add the appropriate offset to convert to absolute
    % time.
    if ~isempty(eyeTrial.blinktimes)
      blinkTimes = [blinkTimes; (eyeTrial.blinktimes(1,:)' + onsetsByEvent{stim_i}(trial_i))];
    end
    
    % Use the generated rasters to find saccade times, and shift by
    % psthPre to convert to indicies from times.
    saccadeBins = find(stimSaccadeArrays{stim_i}(trial_i,:)) - psthParams.psthPre;
    
    if ~isempty(eyeTrial.saccadetimes)
      % Add each saccade to the times
      saccadeTimes = [saccadeTimes; (saccadeBins' + onsetsByEvent{stim_i}(trial_i))];
      
      % for each added saccade, identify the appropriate label (preStim,
      % stim, postStim).
      preStim = saccadeBins <= 0;
      stimDur = (saccadeBins > 0 & saccadeBins <= psthParams.psthImDur)*2;
      postStim = (saccadeBins > psthParams.psthImDur)*3;
      saccadeLabelVec = preStim + stimDur + postStim;
      
      % contactonate to the vector.
      saccadeLabel = [saccadeLabel; saccadeLabelVec'];
      
      % Direction
      saccadeDir = [saccadeDir; eyeTrial.saccadeDirection'];
      
      % Image
      saccadeInds = saccadeBins + psthParams.psthPre;
      slicedData = extractSliceAndPad(trialTemplate, saccadeInds, rasterSaccParams.psthParams.psthPre, rasterSaccParams.psthParams.psthImDur, 'NaN');
      saccadeImg = [saccadeImg; slicedData];
    end
        
  end
end

% Filter based on saccades which took place outside of the stimulusPeriod
if nonStimSaccadesOnly
  saccKeepInd = saccadeLabel == 1 | saccadeLabel == 3;
  saccadeLabel = saccadeLabel(saccKeepInd);
  saccadeImg = saccadeImg(saccKeepInd, :);
  saccadeDir = saccadeDir(saccKeepInd);
  saccadeTimes = saccadeTimes(saccKeepInd);
end

eyeEventTimes = {blinkTimes; saccadeTimes}; 
% figure(); imagesc(saccadeImg);

if saccadeSplit
  % Event Type 2.5 - subEvents occuring in the stimuli, split by whether a
  % subsequent saccade took place.
  saccadeTimes = sort(saccadeTimes);
  subEventsNew = cell(length(subEventNames), 2);
  subEventNull = [onsetsByEventNull, onsetsByEventNull];
  subEventNamesNew = [strcat(subEventNames, '_sacc'), strcat(subEventNames, '_saccNone')];
  for event_i = 1:length(subEventNames)
    
    % Pull event times
    eventOnsets = onsetsByEvent{event_i};
    
    % See how they compare to previously collect saccade times
    eventOnsets2Check = repmat(saccadeTimes, [1, length(eventOnsets)]);
    eventOnsetsOffsetBySacc = eventOnsets2Check - eventOnsets';
    
    % If any number is between 0 and the previously defined window, it means
    % a saccade occured close enough after.
    eventsWSaccades = any(eventOnsetsOffsetBySacc > 0 & eventOnsetsOffsetBySacc < saccadeWindow, 1);
    
    % Store them as distinct events
    subEventsNew{event_i, 1} = eventOnsets(eventsWSaccades);
    subEventsNew{event_i, 2} = eventOnsets(~eventsWSaccades);
    
  end
  
  % Rearrange and overwrite originals
  onsetsByEvent = reshape(subEventsNew', [], 1);
  onsetsByEventNull = reshape(subEventNull', [], 1);
  subEventNames = reshape(subEventNamesNew', [], 1);
end

% Break up saccades based on their direction, sort by stim relation.
saccadeDir(saccadeDir == 9) = 1;
[saccadeTimesDir, saccadeLabelDir, saccImg] = deal(cell(8,1));

for dir_i = 1:8
  % Pull the saccades belonging to the direction
  saccInd = find(saccadeDir == dir_i);
  [sortedDir, sortInd] = sort(saccadeLabel(saccadeDir == dir_i));
  saccadeLabelDir{dir_i}  = sortedDir;
  tmp = saccadeTimes(saccadeDir == dir_i);
  saccadeTimesDir{dir_i} = tmp(sortInd);
  saccImg{dir_i} = saccadeImg(saccInd,:);
end

imgStruct.saccadeByStim = saccImg;
imgStruct.attendedObjData = eyeDataStruct.attendedObjData;

saccDirNames = string(circleInRadians(1:8) * (180/pi)); % Presaccade will have the same times as saccades, but the comparison window for the t test will stretch back.
onsetsByEvent = saccadeTimesDir;
rasterSaccParams.refOffset = 0;
[spikesBySubEvent, spikesEmptyBySubEvent] = alignSpikes(spikesByChannel, onsetsByEvent, ones(length(spikesByChannel),1), rasterSaccParams);

% Have an index for latter labeling
cardinalDirKeep = ~cellfun('isempty', onsetsByEvent);
saccDirNames = saccDirNames(cardinalDirKeep);
saccPerDir = cellfun('length', onsetsByEvent);
cardinalLabelInd = cumsum(saccPerDir);
cardinalLabelInd = cardinalLabelInd - floor(saccPerDir/2);
saccDirYLabel = repmat({''}, [sum(saccPerDir), 1]);
saccDirYLabel(cardinalLabelInd) = cellstr(saccDirNames);

% follow with putting spikesBySubEvent into calcPSTH.
if ~rasterSaccParams.spikeTimes
  spikesBySubEventBinned = calcSpikeTimes(spikesBySubEvent, rasterSaccParams.psthParams);  
  [psthBySubEvent, ~] = calcStimPSTH(spikesBySubEventBinned, spikesEmptyBySubEvent, rasterSaccParams.spikeTimes, rasterSaccParams.psthParams, rasterSaccParams);
else
  [psthBySubEvent, ~] = calcStimPSTH(spikesBySubEvent, spikesEmptyBySubEvent, rasterSaccParams.spikeTimes, rasterSaccParams.psthParams, rasterSaccParams);
end

% Plotting - Figures representing Rasters and PSTH with respect to
% direction of saccades.

for chan_i = 1:length(spikesBySubEventBinned{1})
  saccSelUnitInds = saccadeSelCells(strcmp(selTable.channel, figStruct.channelNames{chan_i}));
  
  % If channel has none selective, skip
  if ~any(saccSelUnitInds)
    continue
  end
  
  for unit_i = 2:length(spikesBySubEventBinned{1}{chan_i})-1 % Only do this for units.
    
    % If the unit appears selective to saccades, plot the normal lines, otherwise skip. 
    if saccSelUnitInds(unit_i) == 1
      spikePattern = 1;
    else
      %spikePattern = 3;
      continue
    end
          
    if length(spikesBySubEventBinned{1}{chan_i}) == 2 && unit_i == 1
      continue
    end
    
    % Raster, color coded
    chanUnitTag = sprintf('%s%s', figStruct.channelNames{chan_i}, figStruct.channelUnitNames{chan_i}{unit_i});
    if 0
      figTitle = sprintf('SaccDir_Raster - %s', chanUnitTag);
      fh = figure('Name', figTitle, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0.4 0.05 0.37 0.95]);
      
      %     subAx = subplot(3, 1, 1:2);
      subAx = axes();
      subAx.YTick = []; subAx.XTick = [];
      rasterColorCoded(subAx, spikesBySubEvent, saccDirNames, rasterSaccParams.psthParams, 500, chan_i, unit_i, struct(), 0, spikePattern);
      xlabel('Time relative to Saccade onset (ms)')
      yticklabels(saccDirYLabel);
      ylabel('Saccade Direction (Degrees)');
      title(strrep(figTitle, '_', ' '));
      delete(findobj(fh, 'Type', 'Legend'));
      saveFigure(figStruct.figDir, figTitle, [], figStruct, figStruct.figTag);
    end

    % PSTH, on the right
    figTitle = sprintf('SaccDir_PSTH - %s', chanUnitTag);
    fh = figure('Name', figTitle, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0.4047 0.0583 0.36 0.38]);
    subAx = axes();
    
    % subAx2 = subplot(3, 1, 3);
    plotPSTH(psthBySubEvent{chan_i}{unit_i}, [], [], rasterSaccParams.psthParams, 'color', 'PSTH', saccDirNames);
    title(strrep(figTitle, '_', ' '));
    xlabel('Time from saccade onset (ms)');
    ylabel('Direction (Deg)');
    subAx.YLabel.FontSize = 14;
    % Save
    saveFigure(figStruct.figDir, figTitle, [], figStruct, figStruct.figTag);
    
  end
end

end

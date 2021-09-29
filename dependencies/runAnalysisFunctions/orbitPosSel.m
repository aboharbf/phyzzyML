function selTable = orbitPosSel(selTable, spikesByChannel, analogIn, psthParams, taskDataAll, figStruct)
% a function which determiens whether a unit is showing selectivity for the
% eyes being in a particular position
% Input:
% - spikeByChannel
% - eye data
% - params for proper alignment
% - figStruct

% Method:
% - Put every spike down on a grid (according to eye position during firing) and then compare the segments as though
% they are independent.

% Identify times of the trial when stimuli is on.
stimStarts = round(taskDataAll.taskEventStartTimes);
stimEnds = round(taskDataAll.taskEventEndTimes);
stimOnVec = zeros(1, length(analogIn));

% Remove NaNs, as these are stim which never turned on (fail during fix)
stimEnds = stimEnds(~isnan(stimStarts));
fixTime = taskDataAll.fixTime(~isnan(stimStarts));
stimStarts = stimStarts(~isnan(stimStarts));

excludeFix = true;
if excludeFix
  stimStarts = stimStarts - fixTime;
end

for ii = 1:length(stimStarts)
  stimOnVec(stimStarts(ii):stimEnds(ii)) = deal(true);
end

excludeStim = true;
if excludeStim
  stimOnTimesAll = find(stimOnVec);
  stimOffTimesAll = find(~stimOnVec);
else
  % If we're not excluding stim, pretend its never on.
  stimOnTimesAll = []; %find(true(size(stimOnVec)));
  stimOffTimesAll = find(true(size(stimOnVec))); % find(true(size(stimOnVec)));
end

% - Grab eye position data and spike data.
% - for every spike a unit discharges, place onto a grid (maybe 20 by 20
% grid).
% - Smooth the grid. 
% - see if any region rises above the rest somehow...
gaussFilt = gausswin(5)*gausswin(5)';
gaussFilt = gaussFilt./sum(gaussFilt(:));
unitIndex = 1;
[spikeEyeCorr, spikeEyeCorrNull] = deal(nan(size(selTable,1),1));

for chan_i = 1:length(figStruct.channelNames)
  
  % Identify units to use
  unitNamesChan = selTable.unitType(strcmp(selTable.channel, figStruct.channelNames{chan_i}));
  units2Check = unique(spikesByChannel(chan_i).units);
  unitCheckMat = logical([eye(length(units2Check)), true(length(units2Check),1)]); % check each unit, then all of them for MU.
  
  for unit_i = 1:size(unitCheckMat, 2)
    
    %  Grab the correct units
    spikeIndex = ismember(spikesByChannel(chan_i).units, units2Check(unitCheckMat(:, unit_i)));
    spikeTimes = round(spikesByChannel(chan_i).times(spikeIndex));
    
    % Remove any spikes from when the stim is on
    spikeTimes = spikeTimes(~ismember(spikeTimes, stimOnTimesAll));
    
    % If there is an odd number, remove a random spike.
    if mod(length(spikeTimes),2) ~= 0
      spikeTimes = spikeTimes(randperm(length(spikeTimes), length(spikeTimes)-1));
    end
    
    % Take halves iteratively and compare them - for real spikes
    corrDist = splitHalvesCorrelation(spikeTimes, analogIn, 1000, gaussFilt);
          
    % Create some random distributions of times when the stimulus was not
    % on.
    corrDistNull = nan(1000,1);
    for rand_i = 1:1000
      spikeTimesNull = stimOffTimesAll(randperm(length(stimOffTimesAll), length(spikeTimes)));
      corrDistNull(rand_i) = splitHalvesCorrelation(spikeTimesNull, analogIn, 1, gaussFilt);
    end
    
    % figure(); histogram(corrDist,100); hold on; histogram(corrDistNull,100); legend('Real', 'Scramble');
    spikeEyeCorr(unitIndex) = mean(corrDist);
    spikeEyeCorrNull(unitIndex) = mean(corrDistNull);
    unitIndex = unitIndex + 1;
    
  end
end

selTable.spikeEyeCorr = spikeEyeCorr;
selTable.spikeEyeCorrNull = spikeEyeCorrNull;

end
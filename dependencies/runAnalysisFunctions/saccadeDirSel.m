function selTable = saccadeDirSel(spikesByEventBinned, eyeBehStatsByStim, psthParams, taskData, selTable)
% A function which uses the circStats toolbox to implement saccade
% selectivity.

psthPre = psthParams.psthPre;
psthImDur = psthParams.psthImDur;

% Inputs:
% - spikesByEventBinned: {stim}{chan}{unit}(trial*bin) structure of binned
% spike counts.
% - eyeBehStatsByStim: a {stim}{trial} structure generated by ClusterFix,
% reporting saccade parameters, including angle in SaccadeClusterValues(:,5))
% - selTable: also an Output, used to report results.

% Important note - spikesByUnitBinned is {chan}{unit}(trial*bin), and
% saccadeTimesStats is saccadeTimesStats.field(trial). This works because
% both the inputs originate from a {event}{trial} structure, and that is
% preserved through out. If this changed, make sure changes are made below

% Variables
scrambleCount = 2000;
preSacc = 200;
postSacc = 0;
cardinalDir = true;       % A switch for analyzing things in terms of cardinal directions vs raw angles.

% Reorganize spikesByEventBinned{event}{chan}{unit} into
% spikesByUnitBinned{chan}{unit}
chanCount = length(spikesByEventBinned{1});
unitCounts = cellfun('length', spikesByEventBinned{1});
trialCounts = cellfun('length', eyeBehStatsByStim);
binCount = size(spikesByEventBinned{1}{1}{1}, 2);

spikesByUnitBinned = cell(chanCount,1);
binnedSpikePerUnitTemplate = nan(sum(trialCounts), binCount);
endInds = cumsum(trialCounts);
startInds = [1; endInds(1:end-1) + 1];
indTab = [startInds, endInds];

for chan_i = 1:chanCount
  spikesByUnitBinned{chan_i} = cell(unitCounts(chan_i),1);
  for unit_i = 1:unitCounts(chan_i)
    
    % Initialize array
    spikesByUnitBinnedTmp = binnedSpikePerUnitTemplate;  
    
    % Loop through each stim, pulling the activity block
    for stim_i = 1:length(spikesByEventBinned)
      spikesByUnitBinnedTmp(indTab(stim_i,1):indTab(stim_i,2),:) = spikesByEventBinned{stim_i}{chan_i}{unit_i};
    end
    
    % Store
    spikesByUnitBinned{chan_i}{unit_i} = spikesByUnitBinnedTmp;

  end
end

% simplify eyeBehStatsByStim into an easily referenced structure, removing
% anything not needed.
saccadeTimesTmp = [eyeBehStatsByStim{:}]';
saccadeTimesTmp = [saccadeTimesTmp{:}];
fields2Remove = fields(saccadeTimesTmp);
fields2Remove = setdiff(fields2Remove, {'saccadetimes', 'SaccadeClusterValues', 'saccadeDirection', 'saccadeAngle'});
saccadeTimesStats = rmfield(saccadeTimesTmp, fields2Remove);

trialCount = 1:length(saccadeTimesStats);
saccPerTrial = cellfun('length', {saccadeTimesStats.saccadeDirection});
trialIndsAll = arrayfun(@(ind) repmat(trialCount(ind), [saccPerTrial(ind), 1]), trialCount, 'UniformOutput', false);
trialIndsAll = vertcat(trialIndsAll{:});

if cardinalDir
  % Cardinal Direction based analysis - convert the labeled direction to a radian value
  saccadeDir = [saccadeTimesStats.saccadeDirection]';
  circleInRadians = linspace(0, 2*pi, 9);
  saccadeDirAll = circleInRadians(saccadeDir)';
  binSpacing = circleInRadians(2);
else
  % Raw angles
  saccadeDirAll = [saccadeTimesStats.saccadeAngle]';
  binSpacing = [];
end

saccadetimesAll = [saccadeTimesStats.saccadetimes];

% Collect the offsets to shift times from stimulus onset aligned to
% fixation dot onset aligned.
fixDotDuration = round(taskData.taskEventStartTimes - taskData.taskEventFixDur)';
trialTimeOffset = fixDotDuration(trialIndsAll);

comparisonType = {'fixOnset', 'stimOnset', 'stim', 'all'};

for comparison_i = 1:length(comparisonType)
  % Filter out saccades which don't happen during the correct period. 
  switch comparisonType{comparison_i}
    case 'fixOnset'
      saccadetimesFix = saccadetimesAll + trialTimeOffset;    % Shift times to the appearance of the fixation dot.
      saccKeep = saccadetimesFix(1,:) > 0 & saccadetimesFix(1,:) < 300;
    case 'stimOnset'
      saccKeep = saccadetimesAll(1,:) > 0 & saccadetimesAll(1,:) < 300;
    case 'stim'
      saccKeep = saccadetimesAll(1,:) > 0 & saccadetimesAll(1,:) < psthImDur;
    case 'all'
      saccKeep = true(1, size(saccadetimesAll,2));
  end
   
  saccadeCount = sum(saccKeep);
  saccadeDir = saccadeDirAll(saccKeep);
  saccadetimes = saccadetimesAll(:, saccKeep) + psthPre;
  trialInds = trialIndsAll(saccKeep);
  
  % Create index for table
  unitTabInd = 1;
  unitSelVec = zeros(size(selTable,1), 1);
  saccadeSpikes = zeros(saccadeCount, 1);

  for chan_i = 1:chanCount
    for unit_i = 1:length(spikesByUnitBinned{chan_i})
      
      % For each saccade, gather spikes which occured during that period.
      for sacc_i = 1:saccadeCount
        %       spikeBin = saccadetimes(1, sacc_i):saccadetimes(2, sacc_i); % Use actual sacc times
        spikeBin = saccadetimes(1, sacc_i)-preSacc:saccadetimes(1, sacc_i)+postSacc; % use method in paper
        spikeBin = spikeBin(spikeBin > 0);
        saccadeSpikes(sacc_i) = sum(spikesByUnitBinned{chan_i}{unit_i}(trialInds(sacc_i), spikeBin));
      end
      
      % Each Saccade has a direction and a count associated. for all the
      % directions, cycle through
      [trueMu, ~, ~] = circ_mean(saccadeDir, saccadeSpikes);
      trueR = circ_r(saccadeDir, saccadeSpikes, binSpacing);
      
      % Scramble the labels 2000 times,
      [scrambMu, scrambR] = deal(nan(scrambleCount,1));
      for scramb_i = 1:scrambleCount
        scrambDir = saccadeDir(randperm(saccadeCount, saccadeCount));
        [scrambMu(scramb_i), ~, ~] = circ_mean(scrambDir, saccadeSpikes);
        scrambR(scramb_i) = circ_r(scrambDir, saccadeSpikes, binSpacing);
      end
      
      % plots
      if 0
        figure()
        % Plot the Mu for the mean of the direction.
        subplot(1,2,1)
        h = histogram(scrambMu);
        hold on
        title(sprintf('true Mu = %s', num2str(trueMu, 2)));
        ylim(h.Parent.YLim);
        plot([trueMu trueMu], ylim(), 'linewidth', 3, 'color', 'r');
        
        subplot(1,2,2)
        h = histogram(scrambR);
        hold on
        title(sprintf('true R = %s', num2str(trueR, 2)));
        ylim(h.Parent.YLim);
        plot([trueR trueR], ylim(), 'linewidth', 3, 'color', 'r');
      end
      
      % Preform a t test, store
      [~, pVal, ~, ~] = ttest2(scrambR, trueR);
      unitSelVec(unitTabInd) = pVal;
      
      % Increment
      unitTabInd = unitTabInd + 1;
      
    end
  end
  
  % store in the selTable
  selTable.(sprintf('saccDir_%s_pVal', comparisonType{comparison_i})) = unitSelVec;
  
end
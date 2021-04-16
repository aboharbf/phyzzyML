function selTable = saccadeSel(spikesByEventBinned, eyeBehStatsByStim, psthPre, selTable)
% A function which uses the circStats toolbox to implement saccade
% selectivity.

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
preSacc = 250;
postSacc = 50;

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
trialInds = arrayfun(@(ind) repmat(trialCount(ind), [saccPerTrial(ind), 1]), trialCount, 'UniformOutput', false);
trialInds = vertcat(trialInds{:});

saccadeDir = [saccadeTimesStats.saccadeDirection]';
saccadeCount = length(saccadeDir);
saccadetimes = [saccadeTimesStats.saccadetimes];
saccadeSpikes = zeros(saccadeCount, 1);

% Filter out saccades which don't happen during stim pres. Update reference
% variables accordingly.
saccKeep = saccadetimes(1,:) > 0;
saccadeDir = saccadeDir(saccKeep);
circleInRadians = linspace(0, 2*pi, 9);
saccadeDir = circleInRadians(saccadeDir)';

saccadetimes = saccadetimes(:, saccKeep) + psthPre;
saccadeSpikes = saccadeSpikes(saccKeep);
trialInds = trialInds(saccKeep);
saccadeCount = sum(saccKeep);

% Create index for table
unitTabInd = 1;
unitSelVec = nan(size(selTable,1), 1);

for chan_i = 1:chanCount
  for unit_i = 1:length(spikesByUnitBinned{chan_i});
    
    % For each saccade, gather spikes which occured during that period.
    for sacc_i = 1:size(saccadetimes,2)
%       spikeBin = saccadetimes(1, sacc_i):saccadetimes(2, sacc_i); % Use actual sacc times
      spikeBin = saccadetimes(1, sacc_i)-preSacc:saccadetimes(1, sacc_i)+postSacc; % use method in paper
      spikeBin = spikeBin(spikeBin > 0);
      saccadeSpikes(sacc_i) = sum(spikesByUnitBinned{chan_i}{unit_i}(trialInds(sacc_i), spikeBin));
    end
    
    % Each Saccade has a direction and a count associated. for all the
    % directions, cycle through
%     [trueMu, ul, ll] = circ_mean(saccadeDir, saccadeSpikes);
    trueR = circ_r(saccadeDir, saccadeSpikes);
    
    % Scramble the labels 2000 times,
    [~, scrambR] = deal(nan(scrambleCount,1));
%     [scrambMu, scrambR] = deal(nan(scrambleCount,1));
    for scramb_i = 1:scrambleCount
      scrambDir = saccadeDir(randperm(saccadeCount, saccadeCount));
%       [scrambMu(scramb_i), ~, ~] = circ_mean(scrambDir, saccadeSpikes);
      scrambR(scramb_i) = circ_r(scrambDir, saccadeSpikes);
    end
    
    % plots
%     figure()
    % Plot the Mu for the mean of the direction.
%     subplot(1,2,1)
%     h = histogram(scrambMu);
%     hold on
%     title(sprintf('true Mu = %d', trueMu));
%     ylim(h.Parent.YLim);
%     plot([trueMu trueMu], ylim(), 'linewidth', 3, 'color', 'r');
    
%     subplot(1,2,2)
    if 0
      figure()
      h = histogram(scrambR);
      hold on
      title(sprintf('true R = %s', num2str(trueR, 2)));
      ylim(h.Parent.YLim);
      plot([trueR trueR], ylim(), 'linewidth', 3, 'color', 'r');
    end
    
    % Preform a t test
    [~, pVal, ~, ~] = ttest2(scrambR, trueR);
    if pVal < 0.05
      unitSelVec(unitTabInd) = 1;
    else
      unitSelVec(unitTabInd) = 0;
    end
    %     [A, B, C, D] = ttest2(scrambMu, trueMu);
    
    unitTabInd = unitTabInd + 1;
    
  end
end

selTable.saccSel = unitSelVec;
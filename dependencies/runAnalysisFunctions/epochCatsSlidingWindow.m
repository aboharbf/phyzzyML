function [selTable] = epochCatsSlidingWindow(spikesByEvent, saccadeByStim, selTable, eventIDs, paradigm, psthParams, params)
% This function performs a sliding window ANOVA, taking into account
% Social Label, Category Labe, and eye behavior label. It adds columns to
% selTable corresponding to the p-Val for each factor in each bin, as well
% as the fraction of explained variability in that bin.

stimParamsFilename = params.stimParamsFilename;
preStimSpikePad = (psthParams.movingWin(1)/2) + psthParams.psthPre;
preStimEyePad = psthParams.psthPre;
totalTime = psthParams.psthImDur + psthParams.psthPost;

binSize = params.binSize;
binStep = params.binStep;
alpha = params.alpha;                           % alpha value to set while looking for units.

% Simple variables necessary for later
chanCount = length(spikesByEvent{1});
unitCounts = cellfun('length', spikesByEvent{1});
trialCounts = cellfun(@(x) size(x{1}{1}, 1), spikesByEvent);
groupsAndLabels = params.(paradigm);
comparisonLabel = groupsAndLabels.comparisonLabel;
comparisonCategoryLabel = groupsAndLabels.comparisonCategoryLabels;

%% Binning for regression

binStarts = 0:binStep:totalTime-binSize;
binStarts(1) = 1;
binEnds = binStarts + binSize;
timeBinSpikes = [binStarts', binEnds'] + preStimSpikePad;
timeBinEyes = [binStarts', binEnds'] + preStimEyePad;
binnedSpikes = initNestedCellArray(spikesByEvent);
binnedEye = cell(length(eventIDs),1);

% Collect spikes for every Bin.
for event_i = 1:length(eventIDs)
  
  % Bin Spikes
  for chan_i = 1:chanCount
    for unit_i = 1:unitCounts(chan_i)
      
      % Initialize structure for storage
      binnedCountsPerUnit = nan(trialCounts(event_i), size(timeBinSpikes,1));
      
      for epoch_i = 1:size(timeBinSpikes,1)
        binnedCountsPerUnit(:, epoch_i) = sum(spikesByEvent{event_i}{chan_i}{unit_i}(:, timeBinSpikes(epoch_i, 1): timeBinSpikes(epoch_i, 2)), 2);
        %spikeCountsByImageByEpoch{epoch}{channel}{unit}{event}.rates = trials*1
      end
      
      % Store
      binnedSpikes{event_i}{chan_i}{unit_i} = binnedCountsPerUnit;
      
    end
  end
  
  binnedEyeByStim = nan(trialCounts(event_i), size(timeBinSpikes,1));
  
  % Bin eye behavior
  for epoch_i = 1:size(timeBinEyes,1)
    for trial_i = 1:trialCounts(event_i)
      binnedEyeByStim(:, epoch_i) = round(mean(saccadeByStim{event_i}(:, timeBinEyes(epoch_i, 1): timeBinEyes(epoch_i, 2)), 2));
    end
  end
  
  binnedEye{event_i} = binnedEyeByStim;
  
end

% Find unique eye signal
eyeLevels = unique(vertcat(binnedEye{:}));
eyeLabels = {'Fix', 'preSacc', 'Sacc', 'Blink'};

%% For every group

for group_i = 1:length(comparisonLabel)
  
  % Extract labels + group labels
  %   catLabel = comparisonLabel{group_i};
  %   stimCat = contains(comparisonCategoryLabel{group_i}, 'ing') + 1;
  stimLabels = comparisonCategoryLabel{group_i};
  stimLabelCount = length(stimLabels);
  stimCatLabels = {'Social', 'nonSocial'};
  
  % Use plotInd
  plotParams.stimParamsFilename = stimParamsFilename;
  plotParams.plotLabels = comparisonCategoryLabel{group_i};
  plotParams.outLogic = 1;
  stimCatInd = plotIndex(eventIDs, plotParams);
  
  % Rates for tests
  spikesPerLabel = cell(stimLabelCount, size(timeBinSpikes,1), sum(unitCounts));
  eyePerLabel = cell(stimLabelCount, size(timeBinSpikes,1));
  %   [spikesPerLabel{:}] = deal(nan(150,1));
  
  % Initialize matrix of entries for results.
  allUnitInd = 1;
  
  % Restructure the spike data for use
  for chan_i = 1:chanCount
    for unit_i = 1:unitCounts(chan_i)
      for cat_i = 1:stimLabelCount
        catInd = stimCatInd(:, cat_i);
        
        % Extract the Rates across trials
        stim2Combine = horzcat(binnedSpikes{catInd});
        stim2Combine = stim2Combine(chan_i, :);
        
        % Extract the Eye behavior
        stimEye2Combine = vertcat(binnedEye{catInd});
        
        for ep_i = 1:size(timeBinSpikes,1)
          
          % Pull spike data across this epoch.
          stimRates = cellfun(@(x) x{unit_i}(:, ep_i), stim2Combine, 'UniformOutput', 0)';
          stimRates = vertcat(stimRates{:});
          
          % Collapse eye data across this epoch.
          stimEye = stimEye2Combine(:, ep_i);
          
          % Store in the larger array
          spikesPerLabel{cat_i, ep_i, allUnitInd}= stimRates;
          eyePerLabel{cat_i, ep_i} = stimEye;
        end
        
      end
      
      % Increment index
      allUnitInd = allUnitInd + 1;
      
    end
  end
  
  % For every unit, do the comparison
  trialsPerCat = cellfun('length', spikesPerLabel(:, 1, 1));
  trialsTotal = sum(trialsPerCat);
  
  % Labels for category
  labelNum = cell(length(trialsPerCat),1);
  for lab_i = 1:stimLabelCount
    labelNum{lab_i} = repmat(lab_i, [trialsPerCat(lab_i), 1]);
  end
  labelNum = vertcat(labelNum{:});
   
  % Run the comparison
  pMat = nan(size(spikesPerLabel,3), size(spikesPerLabel,2), 5);
  explainedVar = nan(size(spikesPerLabel,3), size(spikesPerLabel,2), 5);
  totalError = nan(size(spikesPerLabel,3), size(spikesPerLabel,2), 3);

  % Add in Social vs non-social to make a nested model
  socLabelNum = (labelNum >= 5) + 1;
 
  for ep_i = 1:size(spikesPerLabel,2);
    for unit_i = 1:size(spikesPerLabel,3)
      % Find spikes for this epoch and unit
      spikeCounts = spikesPerLabel(:, ep_i, unit_i);
      spikeCounts = vertcat(spikeCounts{:});
      
      % Find eye movements for this epoch.
      labelEye = vertcat(eyePerLabel{:, ep_i});
      
      % Perform Test
      [pMat(unit_i, ep_i, :), B, ~, ~] = anovan(spikeCounts, [{stimCatLabels(socLabelNum)}; {stimLabels(labelNum)}; {eyeLabels(labelEye)}], 'model', 'interaction', 'nested', ...
        [0 0 0; 1 0 0; 0 0 0;], 'display', 'off', 'varnames', {'SocialLabel', 'Category', 'Eye'}, 'alpha', alpha);
      
      % labelsUsed
      labelsForTable = B(2:end-2, 1);
      labelsForTable = strrep(labelsForTable, '(', '_');
      labelsForTable = strrep(labelsForTable, ')', '_');
      labelsForTable = strrep(labelsForTable, '*', '_');
      
      SS_PerFactor = [B{2:end, 2}];
      
      % Store the explained variance of the total variance for each factor.
      for SS_i = 1:length(labelsForTable)
        explainedVar(unit_i, ep_i, SS_i) = SS_PerFactor(SS_i)/SS_PerFactor(end);
      end
      
      totalError(unit_i, ep_i, 1) = sum(SS_PerFactor(1:end-2));
      totalError(unit_i, ep_i, 2) = SS_PerFactor(end-1);
      totalError(unit_i, ep_i, 3) = SS_PerFactor(end);

    end
  end
  
  % Plotting
  if 0
    pMatSqueezed = pMat;
    pMatSqueezed = pMatSqueezed < 0.05
%     pMatSqueezed = -log10(pMat);
%     pMatSqueezed(pMatSqueezed > 5) = 5;
%     pMatSqueezed(pMatSqueezed < 2) = 2;
    for p_i = 1:size(pMatSqueezed, 3)
      subplot(1, size(pMatSqueezed, 3), p_i)
      imagesc(pMatSqueezed(:,:,p_i))
      colorbar()
    end
    
  end
  
  % Max stretch saving
  if 0
    sigMat = (pMat < alpha);
    maxStetchLength = findStretches(sigMat);
  end
  
  % Add to the output table
  pValPerBin = mat2cell(pMat, ones(size(pMat,1), 1), size(pMat,2), ones(size(pMat,3), 1));
  ExpVarPerBin = mat2cell(explainedVar, ones(size(explainedVar,1), 1), size(explainedVar,2), ones(size(explainedVar,3), 1));
  
  for tab_i = 1:length(labelsForTable)
    selTable.([labelsForTable{tab_i} 'pVal']) = pValPerBin(:, :, tab_i);
    selTable.([labelsForTable{tab_i} 'VarExp']) = ExpVarPerBin(:, :, tab_i);
  end
    
end

end

% MDS Notepad

% Load the binned data
binnedData = load('D:\DataAnalysis\batchAnalysis\NeuralDecodingTB\NaturalSocial\rasterData\rasterData_binned_150ms_bins_50ms_sampled.mat');
plotIndParams.stimParamsFilename = 'C:\Users\aboha\Onedrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\stimParamFileLib\StimParamFileSocialVids_Full.mat';

% binned_data - N Site cell array, where each entry is trials*bins
% binned_labels - a struct with some labels, each of which is a N site cell
% array, where each entry is trials*1.
% binned_site_info -  a struct with some labels, each of which is a N * 1
% structure with different values each can take on (many are logical).

unitType = binnedData.binned_site_info.UnitType';
muaInd = contains(unitType, 'M') & binnedData.binned_site_info.socIntSel_any;

trialsPerSite = cellfun('length', binnedData.binned_labels.socialCat)';
binCount = size(binnedData.binned_data{1}, 2);
binnedDataPerGroup = nan(sum(trialsPerSite), binCount);

% Just grab MUA related things
stimPerSite = binnedData.binned_labels.stimuli;
labels2Include = {{'objects', 'goalDirected', 'idle', 'socialInteraction'}};

for group_i = 1:length(labels2Include)
  
  % Look for the indices to include
  plotIndParams.plotLabels = labels2Include(group_i);
  indices2Grab = cell(sum(muaInd), length(labels2Include{group_i}));
  
  % Initialize array for storing data.
  binnedDataArray = cell(length(labels2Include{group_i}),1);
  [binnedDataArray{:}] = deal(binnedDataPerGroup);
  
  % Collect the indices which will be used to collect data.
  for site_i = find(muaInd)'
    
    % Identify trials to be collected from the binned_data.
    stimuliIndPerLabel = plotIndex(stimPerSite{site_i}, plotIndParams);
    [A, ~, C] = unique(stimuliIndPerLabel);
    a_counts = accumarray(C, 1); % 1 is 0, and shouldn't be counted.
    counts2Keep = find(A ~= 0)';
    a_counts = a_counts(counts2Keep);
    
    % Grab the number of labels which are 
    trialsPerLabel = min(a_counts);
    
    for label_i = counts2Keep
      % Identify trials to collect
      trialInd = C == label_i;
      
      % Collect their indicies
      label_inds = find(trialInd);
      
      % If there are too many trials, grab a random set of trials to reach the minimum
      if length(label_inds) > trialsPerLabel
        label_inds = label_inds(randperm(length(label_inds), trialsPerLabel));
      end
      
      % Store in the larger matrix
      indices2Grab{site_i, A(label_i)} = label_inds;
            
    end
    
  end
  
  % Find the sites
  trialsPerSite2Gather = cellfun('length', indices2Grab(:,1));
  
  % Cycle through each site, and label
  labelCountTotal = length(labels2Include{group_i});
  for site_i = 1:size(indices2Grab,1)
    
    % Initialize vectors for saving
    tmpArray = nan(labelCountTotal, binCount);
    
    for label_i = 1:labelCountTotal
      
      % Find the indices to index into the data
      indices2Use = indices2Grab{site_i, label_i};
      
      % Collect the data using the indicies
      data2Store = mean(binnedData.binned_data{site_i}(indices2Use,:));
      
      % Store it in this spot
%       binnedDataArray{label_i}(nanInd:nanInd+size(data2Store, 1)-1, :) = data2Store;
      tmpArray(label_i,:) = data2Store;
      
    end
    
    % Normalize all the activity from this specific site, across
    % conditions.
    tmpArray = normalize(tmpArray);
    
    for label_i = 1:labelCountTotal
      % Find the next slot to put data into.
      nanInd = find(isnan(binnedDataArray{label_i}(:,1)),1);
      binnedDataArray{label_i}(nanInd:nanInd+size(data2Store, 1)-1, :) = tmpArray(label_i,:);
    end
    
  end
  
  % Remove the NaNs
  for grp_i = 1:length(binnedDataArray)
    nanInd = find(isnan(binnedDataArray{grp_i}(:,1)),1);
    binnedDataArray{grp_i} = binnedDataArray{grp_i}(1:nanInd-1, :);
  end
  
  % Combine across the arrays
  labelCount = size(binnedDataArray{1}, 1);
  
  binnedDataStack = vertcat(binnedDataArray{:});
  binnedLabelStack = repmat(1:length(binnedDataArray), [labelCount, 1]);
  binnedLabelStack = binnedLabelStack(:);
  
  % pca time
  dataIn = binnedDataStack';
  [coeff,score,latent] = pca(binnedDataStack);
  eigenvalues = diag(latent);
  
  plot(cumsum(latent)/sum(latent))
  
  disp('done w/ binning data');
  
  
end
function popDistPlot(params)
% Function which looks in the appropriate directory for pcaData, and plots
% it.

% Points to label in trajectory
tmpColorMat = [0 0 0.7; 0.8 0 0]; % Segment colors;
% string2Plot = {'Stim Onset', 'Stim End', 'Reward'};
string2Plot = {'Stim Onset', 'Stim End'};
pointsForLabels = [0 2400];
pcaInputType = {'meanResponse', 'singleTrial'};

pcCountTotal = 12;
recordSpinningVid = false;              
removeRewardPeriod = true;
removeFixPeriod = true;
pcaDataPlots = false;
excludeIndividualTimeCourse = true;

% Define variables
if removeRewardPeriod && removeFixPeriod
  outputDir = fullfile(params.dimRedParams.outputDir, 'PCA_justStim');
elseif removeRewardPeriod
  outputDir = fullfile(params.dimRedParams.outputDir, 'PCA_noReward');
else
  outputDir = fullfile(params.dimRedParams.outputDir, 'PCA');
end

if ~exist(outputDir, 'dir')
  mkdir(outputDir)
end

% colors = 'rbgkcmyk';
colors = [1 0 0; 0.8 0 0; 0.6 0 0; 0.4 0 0;...
          0 1 0; 0 0.8 0; 0 0 1; 0 0 0.6;]; % Segment colors;

markerSize = 30;

% dimRedLabel = {'categories_Combo_NaturalSocial_MUA_SNT',  'categories_Combo_NaturalSocial_Units_SNT', ...
%   'categories_Combo_headTurnCon_MUA_SNT', 'categories_Combo_headTurnCon_Units_SNT'};

if removeRewardPeriod && removeFixPeriod
  dimRedPCs = {[1 2 4]; [1 2 3]};
  cameraPos = {[2458.88 -1063.7 1628.3]; [2458.88 -1063.7 1628.3]};
else
  dimRedPCs = {[3 4 5]; [3 4 5]};
  cameraPos = {[-33.2969 -140.4141   92.1205]; [-33.2969 -140.4141   92.1205]};
end

% Identify files of interest
dimRedParams = dir(fullfile(params.dimRedParams.preprocDir, '*.mat'));
dimRedParams = fullfile({dimRedParams.folder}, {dimRedParams.name})';

fileNames = extractBetween(dimRedParams, 'dimRedData_', '.mat');
titles = strrep(fileNames, '_', ' ');
pcVar = nan(12, length(dimRedParams));

activityProcessTag = 'SNT';
dimRedLabel = {'categories_Combo_NaturalSocial_Units', 'categories_Combo_headTurnCon_Units'};
dimRedLabel = strcat(dimRedLabel, '_', activityProcessTag);

dimRedFiltInd = contains(titles, 'meanResponse') & contains(titles, '200') & contains(titles, activityProcessTag) & contains(titles, 'categories') & contains(titles, 'Combo');
dimRedParams = dimRedParams(dimRedFiltInd);
titles = titles(dimRedFiltInd);

for data_i = 1:length(dimRedParams)
  
  % Load the data
  load(dimRedParams{data_i})
  
  binStarts_shifted = binLabelInfo.the_bin_start_times_shift;
  binLabels = binLabelInfo.bins_to_label;
  
  if removeRewardPeriod || removeFixPeriod
    
    % if you need to remove some fraction of the activity, do so below.
    % First, remove the binStarts_shifted.
    keepInd = true(size(binStarts_shifted));
    if removeRewardPeriod
      keepInd = keepInd & binStarts_shifted < 2800;      
    end
    
    if removeFixPeriod
      keepInd = keepInd & binStarts_shifted >= 0;      
    end
        
    label_keep = binLabels >= find(keepInd, 1, 'first') & binLabels <= find(keepInd, 1, 'last');
    binStarts_shifted = binStarts_shifted(keepInd);
    
    binLabelInfo.points_to_label = binLabelInfo.points_to_label(label_keep);  
    binLabelInfo.bins_to_label = interp1(binStarts_shifted, 1:length(binStarts_shifted), binLabelInfo.points_to_label);
    binLabelInfo.x_for_lines = interp1(binStarts_shifted, 1:length(binStarts_shifted), binLabelInfo.points_for_lines);
    
    binnedDataArray = cellfun(@(x) x(:, keepInd), binnedDataArray, 'UniformOutput', false);
    binCount = sum(keepInd);
       
  end
  
  binLabelInfo.the_bin_start_times_shift = binStarts_shifted;
  binnedDataArrayMean = cellfun(@(x) mean(x, 1), binnedDataArray, 'UniformOutput', false);
  binnedDataArrayMean = vertcat(binnedDataArrayMean{:});
  
  correlationCoefMat = corr(binnedDataArrayMean');
  
  tmpBlock = ones(4,4);
  comparisonBlock = [[tmpBlock, tmpBlock.*2]; [tmpBlock.*3, tmpBlock.*4]];
  compMean = nan(2,2);
  for ii = 1:4
    ind2Samp = find(comparisonBlock == ii);
    compMean(ii) = mean(correlationCoefMat(ind2Samp));
  end
  
  % Plot it up
  figure();
  axesH = axes();
  imagesc(correlationCoefMat);
  axesH.XTickLabel = labels;
  axesH.YTickLabel = labels;
  numString = sprintf('S-S: %s NS-NS: %s S-NS: %s', num2str(compMean(1),3), num2str(compMean(4),3), num2str(compMean(2),3));
  title(sprintf('Mean Population Distances (via Correlation Coefficient) - %s', numString))
  axesH.FontSize = 14;
  
  % Perform PCA on the data
  [coeff, score, ~, ~, explained, ~] = pca(binnedDataStack');

  
end


end
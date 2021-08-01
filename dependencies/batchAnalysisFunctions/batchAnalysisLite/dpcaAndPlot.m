function dpcaAndPlot(params)
% Function which looks in the appropriate directory for pcaData, and plots
% it.

path2dPCA = 'C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\dimRed\dPCA';
addpath(genpath(path2dPCA))

% Identify files of interest
dimRedFiles = dir(fullfile(params.dimRedParams.preprocDir, '*.mat'));
dimRedFiles = fullfile({dimRedFiles.folder}, {dimRedFiles.name})';
pcCount = 12;

% Define variables
outputDir = fullfile(params.dimRedParams.outputDir, 'figures');
fileNames = extractBetween(dimRedFiles, 'dimRedData_', '.mat');
testName = extractBefore(fileNames, '_');
titles = strrep(fileNames, '_', ' ');
colors = 'rbgkcmyk';
markerSize = 30;

for data_i = 1:length(dimRedFiles)
  
  % Load the data
  load(dimRedFiles{data_i}, 'binnedDataArray', 'binCount', 'paradigmLabel', 'unitType', 'labels', 'binLabelInfo')
  
  % Reshape appropriately
  binnedDataStack = cat(3, binnedDataArray{:});
  %   binnedDataStack = permute(binnedDataStack, [1, 3, 2]);
  binnedDataStack(isnan(binnedDataStack)) = 0;
  margNames = {'Time', 'Category', 'Mystery'};
  
  flip = true;
  if flip
    binnedDataStack = permute(binnedDataStack, [1, 3, 2]);
%     margNames = margNames([2, 1, 3]);
  end
  
  combinedParams = {{1, [1 2]}, {2, [2 1]}};
  margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
  time = (1:size(binnedDataStack,2)) / 10;
  timeEvents = [1.2 4];
  
  %% Step 1: PCA of the dataset

if 0
  X = binnedDataStack(:,:);
  X = bsxfun(@minus, X, mean(X,2));
  
  [W,~,~] = svd(X, 'econ');
  W = W(:,1:20);
  
  % minimal plotting
  dpca_plot(binnedDataStack, W, W, @dpca_plot_default);
  
  % computing explained variance
  explVar = dpca_explainedVariance(binnedDataStack, W, W, ...
    'combinedParams', combinedParams);
  
  % a bit more informative plotting
  dpca_plot(binnedDataStack, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);
end

%%
  
  % run dPCA
  [W, V, whichMarg] = dpca(binnedDataStack, 15);
  
  explVar = dpca_explainedVariance(binnedDataStack, W, V, 'combinedParams', combinedParams);
  
  dpca_plot(binnedDataStack, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);

end
end
function pcaAndPlot(params)
% Function which looks in the appropriate directory for pcaData, and plots
% it.

% Points to label in trajectory
tmpColorMat = [0 0 0.7; 0.8 0 0]; % Segment colors;
% string2Plot = {'Stim Onset', 'Stim End', 'Reward'};
string2Plot = {'Stim Onset', 'Stim End'};
pointsForLabels = [0 2400];
pcaInputType = {'meanResponse', 'singleTrial'};

pcCountTotal = 12;
recordSpinningVid = true;              
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
  dimRedPCs = {[1 2 4]; [2 3 5]};
  cameraPos = {[2458.88 -1063.7 1628.3]; [2458.88 -1063.7 1628.3]; [2458.88 -1063.7 1628.3]; [2458.88 -1063.7 1628.3];...
    [2458.88 -1063.7 1628.3]; [2458.88 -1063.7 1628.3]; [2458.88 -1063.7 1628.3]; [2458.88 -1063.7 1628.3];...
    [2458.88 -1063.7 1628.3]; [2458.88 -1063.7 1628.3]; [2458.88 -1063.7 1628.3]; [2458.88 -1063.7 1628.3]};
else
  dimRedPCs = {[3 4 5]; [3 4 5]};
  cameraPos = {[-33.2969 -140.4141   92.1205]; [-33.2969 -140.4141   92.1205]; [-33.2969 -140.4141   92.1205]; [-33.2969 -140.4141   92.1205]};
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

keepInd = contains(dimRedParams, 'taskModStim');
dimRedParams = dimRedParams(keepInd);
titles = titles(keepInd);

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
  binnedDataStack = horzcat(binnedDataArray{:});
  binnedDataStack(isnan(binnedDataStack)) = 0;
  
  % Perform PCA on the data
  [coeff, score, ~, ~, explained, ~] = pca(binnedDataStack');
    
  % Plots
  % Plot 1 - Data Input
  if pcaDataPlots
    figTitle = sprintf('Data input to PCA, blocked by type - %s - %s', paradigmLabel, unitType);
    figure('name', figTitle)
    imagesc(binnedDataStack)
    hold on
    ylabel('Unit')
    xlabel('Bin')
    start = 1:binCount:size(score, 1);
    yLims = ylim();
    ylim(yLims);
    for ii = 2:length(start)
      plot([start(ii), start(ii)], yLims, 'color', 'k', 'linewidth', 2)
    end
    title(figTitle);
    saveFigure(outputDir, ['1. ' figTitle], [], params.figStruct, [])
    
    % Plot 2 - Score per unit
    figTitle = sprintf('Score per Unit, blocked by type - %s - %s', paradigmLabel, unitType);
    figure('name', figTitle)
    imagesc(score)
    hold on
    ylabel('Unit')
    xlabel('Bin')
    start = 1:binCount:size(score, 1);
    yLims = ylim();
    ylim(yLims);
    for ii = 2:length(start)
      plot([start(ii), start(ii)], yLims, 'color', 'k', 'linewidth', 2)
    end
    title(figTitle);
    saveFigure(outputDir, ['2. ' figTitle], [], params.figStruct, [])
    
    % Plot 3 - Coefficients by unit.
    figTitle = sprintf('PC Coefficients for each Unit - %s - %s', paradigmLabel, unitType);
    figure('name', figTitle)
    imagesc(coeff)
    hold on
    ylabel('Eigenvector')
    xlabel('Bin start time')
    xticks(binLabelInfo.bins_to_label);
    xticklabels(binLabelInfo.points_to_label)
    xlim([1, length(binLabelInfo.the_bin_start_times_shift)]);
    for ii = 1:length(binLabelInfo.x_for_lines)
      plot([binLabelInfo.x_for_lines(ii), binLabelInfo.x_for_lines(ii)], ylim(), 'linewidth', 4, 'color', [0.2, 0.2, 0.2])
    end
    title(figTitle);
    saveFigure(outputDir, ['3. ' figTitle], [], params.figStruct, [])
    
    % Plot 4 - cumulative sum across PCs.
    figTitle = sprintf('Latent variables cumulative sum - %s - %s', paradigmLabel, unitType);
    figure('name', figTitle)
    plot(cumsum(explained))
    hold on
    ylim([0 100])
    xlim([1 length(explained)])
    xlabel('Principal Component')
    ylabel('Cumulative Explained Variance')
    title(figTitle)
    for ii = 2:length(start)
      plot([start(ii), start(ii)], yLims, 'color', 'k', 'linewidth', 2)
    end
    saveFigure(outputDir, ['4. ' figTitle], [], params.figStruct, [])
  end
  
  % Additional variables for PC trajectory plotting.
  timeAxis = 1:binCount;
  start = 1:binCount:size(score, 1);
  stops = [start(2:end) - 1 size(score, 1)];
  newScore = nan(length(labels), binCount, size(score,2));
  for lab_i = 1:length(labels)
    newScore(lab_i, :, :) = score(start(lab_i):stops(lab_i), :);
  end
  PCind = cellfun(@(x) contains(dimRedParams{data_i}, x), dimRedLabel);
  pc2Plot = dimRedPCs{PCind};
  scoreToPlot = newScore(:, :, pc2Plot);
  
  % Trajectories per Label + Time course by per Label
  figTitle = sprintf('PC 1-3 Trajectories %s', titles{data_i});
  axHands = gobjects(length(labels)+1,1);
  
  if excludeIndividualTimeCourse
    figH = figure('name', figTitle, 'units','normalized','position', [0.0700    0.2300    0.3878    0.5900]);
    axHands(1) = subplot(1, 1, 1, 'CameraPosition', cameraPos{data_i});
  else
    figH = figure('name', figTitle, 'units','normalized','position', [0.07 0.12 0.85 0.7]);
    possiblePositions = true(12, 1);
    bigPlotInds = [3 4 7 8];
    axHands(1) = subplot(3, 4, bigPlotInds, 'CameraPosition', cameraPos{data_i});
    possiblePositions(bigPlotInds) = false;
  end
  
  hold on
    
  for i = 1:length(labels)
    scatter3(scoreToPlot(i,:,1), scoreToPlot(i, :, 2), scoreToPlot(i, :, 3), markerSize, colors(i, :), 'filled');
  end
  
  xlabel(sprintf('PC %d', pc2Plot(1))); ylabel(sprintf('PC %d', pc2Plot(2))); zlabel(sprintf('PC %d', pc2Plot(3))); 
  title(figTitle); legH = legend(labels); grid;
  legH.Position = [0.7480    0.7117    0.1997    0.2655];
  
  % Find the points between that mark transitions
  [pointsForLabelsPlot, binsForLabels, binLabelsNew, colorAxis] = bin2Label(binLabelInfo, pointsForLabels, string2Plot, tmpColorMat);
  if ~excludeIndividualTimeCourse
    for i = 1:length(labels)
      % find the next slot
      subplotInd = find(possiblePositions,1);
      possiblePositions(subplotInd) = false;
      
      % Create the plot
      axHands(i+1) = subplot(3, 4, subplotInd);
      scatter3(scoreToPlot(i, :, 1), scoreToPlot(i, :, 2), scoreToPlot(i, :, 3), markerSize, colorAxis)
      
      % Add Labels
      for string_i =1:length(string2Plot)
        text(scoreToPlot(i, pointsForLabelsPlot(string_i), 1), scoreToPlot(i, pointsForLabelsPlot(string_i), 2), scoreToPlot(i, pointsForLabelsPlot(string_i), 3), binLabelsNew{string_i}, ...
          'HorizontalAlignment', 'left', 'FontSize', 10);
      end
      title(labels{i});
      
    end
    linkprop(axHands, {'XLim', 'YLim', 'ZLim', 'CameraPosition', 'CameraUpVector'});
  end
  % Create video rotating about Axes
  if recordSpinningVid
    % Remove the Title
    figH.Children(2).Title.String = '';
    figH.Color = [1 1 1];

    vidObj = VideoWriter(fullfile(outputDir, strcat(figTitle, '.avi')));
    vidObj.FrameRate = vidObj.FrameRate/2;
    open(vidObj);
    
    % Set initial position
    campos([98.8655  -41.5652   30.3368]);
    view([51.6794 16.4829]);
    
    % Rotation to add
    rotArray = zeros(120,2);
    rotArray(1:120,1) = 2;
    
    for frame_i = 1:120
      % Rotate slightly
      camorbit(rotArray(frame_i,1), rotArray(frame_i,2))
      
      % Get frame, and save
      frame = getframe(gcf);
      writeVideo(vidObj, frame);
    end
    close(vidObj)
  end
  
  % saveFigure(outputDir, ['1. ' figTitle], [], params.figStruct, [])
  continue
  
  % Plot each PC
  figTitle = sprintf('Traj in Pop %s - %d PCs, %s var', titles{data_i}, pcCountTotal, num2str(round(sum(explained(1:pcCountTotal)), 2)));
  figH2 = figure('name', figTitle, 'units', 'normalized', 'position', [0 0 0.74 0.97]);
  sgtitle(figTitle);
  subplotHands = gobjects(pcCountTotal, 1);
  ylimPlots = [round(min(score(:)), 2), ceil(max(score(:)))];
  perRow = 3;
  perCol = ceil(pcCountTotal/perRow);

  for PC = 1:pcCountTotal
    subplotHands(PC) = subplot(perCol, perRow, PC);
    hold on
    for i = 1:length(labels)
      plot(1:binCount, newScore(i, :, PC), 'color', colors(i, :), 'LineWidth', 2);
    end
    
    % Label the plot
    xlim([binLabelInfo.bins_to_label(1), size(newScore, 2)]);
    ylim(ylimPlots)
    
    for line_i = 1:length(binLabelInfo.x_for_lines)
      plot([binLabelInfo.x_for_lines(line_i), binLabelInfo.x_for_lines(line_i)], [ylimPlots(2), ylimPlots(1)], 'color', 'k', 'lineWidth', 3)
    end
    
    if PC >= (pcCountTotal - perRow ) + 1
      xticks(binLabelInfo.bins_to_label);
      xticklabels(binLabelInfo.points_to_label)
      xlabel('Time');
    end
    
    ylabel('a.u.');
    title(sprintf('PC %d - %s%%', PC, num2str(round(explained(PC), 2))))
    hold on
  end
  
  linkprop(subplotHands, {'XLim', 'YLim', 'ZLim', 'CameraPosition', 'CameraUpVector'});
  saveFigure(outputDir, ['2. ' figTitle], [], params.figStruct, [])
  
end


end
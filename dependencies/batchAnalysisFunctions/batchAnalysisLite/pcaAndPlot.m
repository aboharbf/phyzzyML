function pcaAndPlot(params)
% Function which looks in the appropriate directory for pcaData, and plots
% it.

% Identify files of interest
dimRedParams = dir(fullfile(params.dimRedParams.preprocDir, '*.mat'));
dimRedParams = fullfile({dimRedParams.folder}, {dimRedParams.name})';
pcCountTotal = 12;
pcTransition2Label = [0 2800 3000];
recordSpinningVid = true;              

% Only certain files
dimRedParams = dimRedParams(contains(dimRedParams, 'SNT'));

% Define variables
outputDir = fullfile(params.dimRedParams.outputDir, 'PCA');

fileNames = extractBetween(dimRedParams, 'dimRedData_', '.mat');
testName = extractBefore(fileNames, '_');
titles = strrep(fileNames, '_', ' ');
colors = 'rbgkcmyk';
markerSize = 30;
pcVar = nan(12, length(dimRedParams));

dimRedParams = dimRedParams(contains(dimRedParams, 'categories') & contains(dimRedParams, 'Combo'));

dimRedLabel = {'categories_Combo_NaturalSocial_MUA_SNT',  'categories_Combo_NaturalSocial_Units_SNT', ...
               'categories_Combo_headTurnCon_MUA_SNT', 'categories_Combo_headTurnCon_Units_SNT'};
dimRedPCs = {[3 4 9]; [3 4 5];...
             [4 5 6]; [5 6 8];};

for data_i = 1:length(dimRedParams)
  
  % Load the data
  load(dimRedParams{data_i})
  
  binnedDataStack = horzcat(binnedDataArray{:});
  binnedDataStack(isnan(binnedDataStack)) = 0;
  
  % Perform PCA on the data
  [coeff, score, ~, ~, explained, ~] = pca(binnedDataStack');
  
  % Plots
  % Plot 1 - Data Input
  pcaDataPlots = false;
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
  
  % Trajectories per Label + Time course by per Label
  figTitle = sprintf('PC 1-3 Trajectories %s', titles{data_i});
  figH = figure('name', figTitle, 'units','normalized','position', [0.07 0.12 0.85 0.7]);
  axHands = gobjects(length(labels)+1,1);
  
  possiblePositions = true(12, 1);
  bigPlotInds = [3 4 7 8];
  axHands(1) = subplot(3, 4, bigPlotInds, 'CameraPosition', [-33.2969 -140.4141   92.1205]);
  possiblePositions(bigPlotInds) = false;
  hold on
    
  % 
  PCind = cellfun(@(x) contains(dimRedParams{data_i}, x), dimRedLabel);
  pc2Plot = dimRedPCs{PCind};
  
  for i = 1:length(labels)
    scatter3(score(start(i):stops(i), pc2Plot(1)), score(start(i):stops(i), pc2Plot(2)), score(start(i):stops(i), pc2Plot(3)), markerSize, colors(i));
  end
  
  xlabel(sprintf('Principal Component %d', pc2Plot(1))); 
  ylabel(sprintf('Principal Component %d', pc2Plot(2))); 
  zlabel(sprintf('Principal Component %d', pc2Plot(3))); 
  title(figTitle);
  legend(labels);
  grid;
  
  % Find the points between that mark transitions
  bins_to_label = interp1(binLabelInfo.the_bin_start_times_shift, 1:length(binLabelInfo.the_bin_start_times_shift), pcTransition2Label);
  tmpColorMat = [0 0 0.6; 0.8 0 0; 0 0.6 0; 0 0.4 0]; % Segment colors;
  bins_to_label = [1, bins_to_label, stops(1)];
  
  colorAxis = [];
  for ii = 1:length(bins_to_label)-1
    colorAxis = [colorAxis; repmat(tmpColorMat(ii,:), length(bins_to_label(ii):bins_to_label(ii+1)-1), 1)];
%     colorAxis = [; ...
%       repmat(tmpColorMat(2,:), length(bins_to_label(1):bins_to_label(2)), 1); ...
%       repmat(tmpColorMat(3,:), length(bins_to_label(2):stops(1)), 1)];
  end
  colorAxis = [colorAxis; colorAxis(end,:)];
  
  for i = 1:length(labels)
    % find the next slot
    subplotInd = find(possiblePositions,1);
    possiblePositions(subplotInd) = false;

    % Create the plot
    axHands(i+1) = subplot(3, 4, subplotInd);
    scatter3(score(start(i):stops(i), pc2Plot(1)), score(start(i):stops(i),pc2Plot(2)), score(start(i):stops(i),pc2Plot(3)), markerSize, colorAxis)
    text(score(start(i)                    ,pc2Plot(1)), score(start(i)                    ,pc2Plot(2)), score(start(i)                    ,pc2Plot(3)), 'Start', 'HorizontalAlignment', 'left', 'FontSize', 10);
    text(score(start(i) + bins_to_label(2) ,pc2Plot(1)), score(start(i) + bins_to_label(2) ,pc2Plot(2)), score(start(i) + bins_to_label(2) ,pc2Plot(3)), 'Stim Onset', 'HorizontalAlignment', 'left', 'FontSize', 10);
    text(score(start(i) + bins_to_label(3) ,pc2Plot(1)), score(start(i) + bins_to_label(3) ,pc2Plot(2)), score(start(i) + bins_to_label(3) ,pc2Plot(3)), 'Stim End', 'HorizontalAlignment', 'left', 'FontSize', 10);
    text(score(start(i) + bins_to_label(4) ,pc2Plot(1)), score(start(i) + bins_to_label(4) ,pc2Plot(2)), score(start(i) + bins_to_label(4) ,pc2Plot(3)), 'Reward', 'HorizontalAlignment', 'left', 'FontSize', 10);
    title(labels{i});
        
  end
  
  linkprop(axHands, {'XLim', 'YLim', 'ZLim', 'CameraPosition', 'CameraUpVector'});
  saveFigure(outputDir, ['1. ' figTitle], [], params.figStruct, [])

  % Create video rotating about Axes
  if recordSpinningVid
    vidObj = VideoWriter(fullfile(outputDir, strcat(figTitle, '.avi')));
    vidObj.FrameRate = vidObj.FrameRate/2;
    open(vidObj);
    
    % Set initial position
    campos([98.8655  -41.5652   30.3368]);
    view([51.6794 16.4829]);
    
    % Rotation to add
    rotArray = zeros(120,2);
    rotArray(1:120,1) = 3;
    
    for frame_i = 1:120
      % Rotate slightly
      camorbit(rotArray(frame_i,1), rotArray(frame_i,2))
      
      % Get frame, and save
      frame = getframe(gcf);
      writeVideo(vidObj, frame);
    end
    close(vidObj)
  end
  
  %     axHands(1).CameraPosition = [-1 0 180];
  %     axHands(1).CameraUpVector = [0 1 0];
  %     saveFigure(outputDir, ['1a. ' figTitle], [], params.figStruct, [])
  
  % Plot per PC.
  figTitle = sprintf('Paradigm %s - Trajectories in Population %s - top %d PCs, %s var', paradigmLabel, titles{data_i}, pcCountTotal, num2str(round(sum(explained(1:pcCountTotal)), 2)));
  figH2 = figure('name', figTitle, 'units', 'normalized', 'position', [0 0 0.74 0.97]);
  sgtitle(figTitle);
  subplotHands = gobjects(pcCountTotal, 1);
  ylimPlots = [round(min(score(:)), 2), ceil(max(score(:)))];
  
  for PC = 1:pcCountTotal
    subplotHands(PC) = subplot(ceil(pcCountTotal/3),3, PC);
    hold on
    for i = 1:length(labels)
      PrincipleAxis = score(start(i):stops(i), PC);
      plot(1:binCount, PrincipleAxis, colors(i), 'LineWidth', 2);
    end
    
    % Label the plot
    xticks(binLabelInfo.bins_to_label)
    xticklabels(binLabelInfo.points_to_label)
    ylim([ylimPlots])
    
    for line_i = 1:length(binLabelInfo.x_for_lines)
      plot([binLabelInfo.x_for_lines(line_i), binLabelInfo.x_for_lines(line_i)], [ylimPlots(2), ylimPlots(1)], 'color', 'k', 'lineWidth', 3)
    end
    
    xlabel('Time'); ylabel('a.u.');
    title(sprintf('PC %d - %s%%', PC, num2str(round(explained(PC), 2))))
    hold on
  end
  
  linkprop(subplotHands, {'XLim', 'YLim', 'ZLim', 'CameraPosition', 'CameraUpVector'});

  
  saveFigure(outputDir, ['2. ' figTitle], [], params.figStruct, [])
  
end


end
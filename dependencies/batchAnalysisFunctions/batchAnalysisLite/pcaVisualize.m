function pcaVisualize(params)
% Function which looks in the appropriate directory for pcaData, and plots
% it.

% Identify files of interest
pcaDataFiles = dir(fullfile(params.pcaParams.preprocDir, '*.mat'));
pcaDataFiles = fullfile({pcaDataFiles.folder}, {pcaDataFiles.name})';
pcCount = 12;
    
% Define variables
outputDir = fullfile(params.pcaParams.outputDir, 'figures');
fileNames = extractBetween(pcaDataFiles, 'pcaData_', '.mat');
testName = extractBefore(fileNames, '_');
titles = strrep(fileNames, '_', ' ');

for data_i = 1:length(pcaDataFiles)
   
  colors = 'rbgkcmyk';
  markerSize = 30;  
  load(pcaDataFiles{data_i});
  
  if data_i == 1
    timeAxis = 1:binCount;
    start = 1:binCount:size(score, 1);
    stops = [start(2:end) - 1 size(score, 1)];
  end
  
  % Color scatter plot distinguishing clusters
  figTitle = sprintf('PC 1-3 Trajectories %s', titles{data_i});
  figH = figure('name', figTitle, 'units','normalized','position', [0.4214    0.1528    0.5568    0.6917]);
  axHands = gobjects(length(bhvLabels)+1,1);
  
  axHands(1) = subplot(3,1, 1:2, 'CameraPosition', [-33.2969 -140.4141   92.1205]);
  hold on
  
%   bhvLabels = bhvLabels(1:size(score,1)/binCount);
  
  for i = 1:length(bhvLabels)
    scatter3(score(start(i):stops(i), 1), score(start(i):stops(i), 2), score(start(i):stops(i), 3), markerSize, colors(i));
  end
  
  xlabel('1st Principal Component')
  ylabel('2nd Principal Component')
  zlabel('3rd Principal Component')
  title(figTitle);
  legend(bhvLabels);
  grid;
  
  colorAxis = repmat(timeAxis,numel(score(:,1)),1);
  colorAxis = colorAxis(1,:);
  for i = 1:length(bhvLabels)
    axHands(i+1) = subplot(3, length(bhvLabels), sub2ind([length(bhvLabels), 3], i, 3));
    scatter3(score(start(i):stops(i),1),score(start(i):stops(i),2),score(start(i):stops(i),3), markerSize, colorAxis)
    text(score(start(i),1), score(start(i),2),0,['Start'],'HorizontalAlignment','left','FontSize',10);
    title(bhvLabels{i});
    
    if length(bhvLabels) == i
      c = colorbar;
      c.Location = 'manual';
      c.Position = [0.9144    0.1086    0.0106    0.2153];
      ylabel(c, 'Bin (time)')
    end
    
  end
  
  linkprop(axHands, {'XLim', 'YLim', 'ZLim', 'CameraPosition','CameraUpVector'});
  saveFigure(outputDir, ['1. ' figTitle], [], params.figStruct, [])
  
%     axHands(1).CameraPosition = [-1 0 180];
%     axHands(1).CameraUpVector = [0 1 0];
%     saveFigure(outputDir, ['1a. ' figTitle], [], params.figStruct, [])
  
%   Trajectories
    figTitle = sprintf('Paradigm %s - Trajectories in Population %s', paradigmLabel, titles{data_i});
    figH2 = figure('name', figTitle, 'units', 'normalized', 'position', [0 0 0.74 0.97]);
    sgtitle(figTitle);
    for PC = 1:pcCount
      subplot(ceil(pcCount/3),3, PC);
      hold on
      for i = 1:length(bhvLabels)
        PrincipleAxis = score(start(i):stops(i), PC);
        plot(1:binCount, PrincipleAxis, colors(i), 'LineWidth', 2);
      end
  
      % Label the plot
      xticks(binLabelInfo.bins_to_label)
      xticklabels(binLabelInfo.points_to_label)
      yLims = ylim();
      ylim(yLims)
      for line_i = 1:length(binLabelInfo.x_for_lines)
        plot([binLabelInfo.x_for_lines(line_i), binLabelInfo.x_for_lines(line_i)], [yLims(1), yLims(2)], 'color', 'k', 'lineWidth', 3)
      end
  
      xlabel('Time'); ylabel('a.u.');
      title(sprintf('PC %d', PC))
      hold on
    end
  
    saveFigure(outputDir, ['2. ' figTitle], [], params.figStruct, [])
  
end
end
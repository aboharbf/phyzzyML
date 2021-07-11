% perTimePt data projected onto PCs 1-3; labeled by bhvtype

monkeyName = {'Mo', 'Sam'};
dataDir = {'D:\DataAnalysisMo\batchAnalysis\PCA\preProc', 'D:\DataAnalysisSam\batchAnalysis\PCA\preProc'};
% dataFiles = {'PCA_data_allMUA_3smooth', 'PCA_data_SocMUA_3smooth', 'PCA_data_allUUS_3smooth', 'PCA_data_SocUUS_3smooth'};
% dataFiles = {'pcaData_HTC_MUA', 'pcaData_HTC_SocMUA', 'pcaData_HTC_US&U', 'pcaData_HTC_SocUS&U', 'pcaData_NS_MUA', 'pcaData_NS_SocMUA', 'pcaData_NS_US&U', 'pcaData_NS_SocUS&U'};
dataFiles = {'pcaData_NS_MUA', 'pcaData_HTC_MUA'};
titles = {'MUA', 'Soc MUA', 'US&U', 'Soc US&U', 'MUA', 'Soc MUA', 'US&U', 'Soc US&U'};
% dataFiles = {'pcaData_Combo_MUA', 'pcaData_Combo_US&U'};
% titles = {'MUA', 'US&U'};

colors = 'rbgkcmyk';

outputDir = 'D:\DataAnalysis\batchAnalysis\PCA';

if ~exist(outputDir, 'dir')
  mkdir(outputDir);
end

markerSize = 30;
% bhvLabels = {'socialInteraction', 'goalDirected', 'idle', 'objects'};
pcCount = 3;

figStruct.saveFig = 1;      % save the figure in its output directory.           
figStruct.closeFig = 0;     % close the figure once it is saved
figStruct.exportFig = 1;    % export figure using export_fig.
figStruct.saveFigData = 0;  % save data with the figure.
figStruct.noOverWrite = 1;  % If a figure is already there, don't make it again.
verbosity = 'INFO';         %other options, 'DEBUG', 'VERBOSE';

for dir_i = 1:length(dataDir)
for data_i = 1:length(dataFiles)
  
  load(fullfile(dataDir{dir_i}, dataFiles{data_i}));
  
  if data_i == 1
    timeAxis = 1:binCount;
    start = 1:binCount:size(score, 1);
    stops = [start(2:end) - 1 size(score, 1)];
  end
  
  % Color scatter plot distinguishing clusters
  figTitle = sprintf('%s - %s - %s PC 1-3 time evolution bhvtype - %s', monkeyName{dir_i}, paradigmLabel, titles{data_i}, strjoin(bhvLabels, '-'));
  figH = figure('name', figTitle, 'units','normalized','position', [0.4214    0.1528    0.5568    0.6917]);
  axHands = gobjects(length(bhvLabels)+1,1);

  axHands(1) = subplot(3,1, 1:2, 'CameraPosition', [-33.2969 -140.4141   92.1205]);
  hold on
  
  bhvLabels = bhvLabels(1:size(score,1)/binCount);

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
  saveFigure(outputDir, ['1. ' figTitle], [], figStruct, [])
  
  
%   axHands(1).CameraPosition = [-1 0 180];
%   axHands(1).CameraUpVector = [0 1 0];
%   saveFigure(outputDir, ['1a. ' figTitle], [], figStruct, [])  
%   
  % Trajectories
%   figTitle = sprintf('Paradigm %s - Trajectories in Population %s', paradigmLabel, titles{data_i});
%   figH2 = figure('name', figTitle, 'units', 'normalized', 'position', [0.25    0.3    0.86    0.45]);
%   sgtitle(figTitle);
%   for PC = 1:pcCount
%     subplot(1,pcCount,PC);
%     hold on
%     for i = 1:length(bhvLabels)
%       PrincipleAxis = score(start(i):stops(i), PC);
%       plot(1:binCount,PrincipleAxis, colors(i), 'LineWidth', 2);
%     end
%     
%     % Label the plot
%     xticks(binLabelInfo.bins_to_label)
%     xticklabels(binLabelInfo.points_to_label)
%     yLims = ylim();
%     ylim(yLims)
%     for line_i = 1:length(binLabelInfo.x_for_lines)
%       plot([binLabelInfo.x_for_lines(line_i), binLabelInfo.x_for_lines(line_i)], [yLims(1), yLims(2)], 'color', 'k', 'lineWidth', 3)
%     end
%     
%     xlabel('Time'); ylabel('a.u.');
%     title(sprintf('PC %d : bhvtype %s', PC, bhvLabels{i}))
%     hold on
%   end
%   
%   saveFigure(outputDir, ['2. ' figTitle], [], figStruct, [])
  
  
end
end
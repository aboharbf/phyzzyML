function selectivityPerGridHole(spikePathBankParadigm, selTable, params)


[selIndGridMat, labelArray] = selInd2GridMat(selTable);

barPlotParams = params.selParam;

% Analyze the batch sheet
batchSheet = 'D:\DataAnalysis\BatchRunResults.xlsx';
phyzzyDir = 'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML';
runList = strcat(selTable.dateSubj, selTable.runNum);
unitTypes = barPlotParams.UnitTypes;
unitTypePlot = barPlotParams.UnitTypePlot;
outDir = fullfile(barPlotParams.outputDir, 'selectivityPerGridHole');

addpath(genpath(phyzzyDir));

% Will be used to normalize counts, need to match unitTypes.
gridUnitCounts = cat(3, gridChannels, gridUnits);

% Recover the events represented
uniqueEvents = selTable.Properties.VariableNames(contains(selTable.Properties.VariableNames, '_selInd'));

gridHoleSelectivity = zeros([size(gridRuns) length(uniqueEvents) length(unitTypes)]);

for ii = 1:length(gridHoles)
  % Identify spot to add
  aInd = strcmp(gridHoles{ii}(1), gridHoleA);
  bInd = strcmp(gridHoles{ii}(2:end), gridHoleB);
  
  %   sprintf('Paradigm %s has %d runs', paradigms{ii}, gridRuns);                     % Runs
  fprintf('Grid Hole %s has %d channels \n', gridHoles{ii}, gridChannels(aInd, bInd));                         % Channels
  fprintf('Grid Hole %s has %d units \n', gridHoles{ii}, gridUnits(aInd, bInd));     % Get a unit count per paradigm
  
  % find the slice of the selTable which represents this gridHole.
  selTableHole = selTable(strcmp(selTable.gridHole, gridHoles{ii}), :);
  
  % Find headTurn per grid hole.
  for event_i = 1:length(uniqueEvents)
    for unit_i = 1:length(unitTypes)
      gridHoleSelectivity(aInd, bInd, event_i, unit_i) = sum(selTableHole.(uniqueEvents{event_i}) ~= 0 & contains(selTableHole.unitType, unitTypes{unit_i}));
      fprintf('Grid Hole %s has %d units selective for %s', gridHoles{ii}, gridHoleSelectivity(aInd, bInd, event_i, unit_i), uniqueEvents{event_i});
      fprintf('\n')
    end
  end
  
end

% Get normalized units for each paradigm, in each hole
gridHoleSelectivityNormalized = gridHoleSelectivity;
for unit_i = 1:length(unitTypes)
  for event_i = 1:length(uniqueEvents)
    gridHoleSelectivityNormalized(:,:, event_i, unit_i) = round(gridHoleSelectivityNormalized(:,:, event_i, unit_i) ./ gridUnitCounts(:,:, unit_i) * 100, 1);
  end
end

%% Figure time
% Reference gridHoleSelectivity, gridHoleSelectivityNormalized

% Runs per hole
figTitle = 'Runs per Grid Hole, Total';
h = figure();
himg = imagesc(gridRuns);
colorbar()

axesH = findobj(h, 'Type', 'Axes');
axesH.YTick = 1:3;
axesH.YTickLabel = gridHoleA;
axesH.XTick = 1:4;
axesH.XTickLabel = gridHoleB;
title(figTitle)

saveFigure(outDir, figTitle, [], params.figStruct, []);

% Units per hole
paradigmPlotName = sprintf('Units per Grid Hole, %s Paradigm', barPlotParams.paradigm);
h = figure('Name',paradigmPlotName,'NumberTitle','off','units','normalized','outerposition', params.figStruct.figPos);
himg = imagesc(round(mean(gridUnits, 3)));
colorbar()

h.Children(2).YTick = 1:length(gridHoleA);
h.Children(2).YTickLabel = gridHoleA;
h.Children(2).XTick = 1:length(gridHoleB);
h.Children(2).XTickLabel = gridHoleB;
title(paradigmPlotName)
saveFigure(outDir, paradigmPlotName, [], params.figStruct, []);

% Runs/Units per hole per paradigm (Figure w/ 4 subplots)

% Unit selectivity in each Paradigm
dataArray = [{gridHoleSelectivity}, {gridHoleSelectivityNormalized}];
normTag = [{''}, {' Normalized'}];
uniqueEventsPlot = strrep(uniqueEvents, '_', ' ');

for data_i = 1:length(dataArray)
  for unit_i = 1:length(unitTypes)
    %     TopPlotName = sprintf("Units selective for events per Grid Hole, %s %s", unitTypePlot{unit_i}, normTag{data_i});
    %     h = figure('Name', TopPlotName, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', figStruct.figPos);
    %     sgtitle(TopPlotName)
    
    for event_i = 16:length(uniqueEvents)
      %       subplot(4, 5, event_i)
      TopPlotName = sprintf("Units selective for event %s per Grid Hole, %s %s", uniqueEventsPlot{event_i}, unitTypePlot{unit_i}, normTag{data_i});
      h = figure('Name', TopPlotName, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', params.figStruct.figPos);
%       sgtitle(TopPlotName)
      plotName = sprintf("%s", strrep(uniqueEvents{event_i}, '_', ' '));
      imagesc(dataArray{data_i}(:, :, event_i, unit_i));
      colorbar()
      
      h.Children(2).YTick = 1:length(gridHoleA);
      h.Children(2).YTickLabel = gridHoleA;
      h.Children(2).XTick = 1:length(gridHoleB);
      h.Children(2).XTickLabel = gridHoleB;
      h.Children(2).FontSize = 16;
      title(TopPlotName)
      
      saveFigure(outDir, TopPlotName, [], params.figStruct, []);
      
    end

  end
end

disp('y')
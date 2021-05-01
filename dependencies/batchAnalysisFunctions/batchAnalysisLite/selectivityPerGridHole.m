function selectivityPerGridHole(spikePathBankParadigm, barPlotParams, selTable)
% Accepts a selTable structure with the 

% Analyze the batch sheet
batchSheet = 'D:\DataAnalysis\BatchRunResults.xlsx';
phyzzyDir = 'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML';
runList = strcat(selTable.dateSubj, selTable.runNum);
unitTypes = barPlotParams.UnitTypes;
unitTypePlot = barPlotParams.UnitTypePlot;
outDir = fullfile(params.outputDir, 'selectivityPerGridHole');

if ~exist(outDir, 'dir')
  mkdir(outDir);
end

% Add combo Events
selTable = expandSelTableComboEvents(selTable, barPlotParams);

addpath(genpath(phyzzyDir));

figStruct = struct();
figStruct.saveFig = 1;                
figStruct.closeFig = 0;               
figStruct.exportFig = 1;             
figStruct.saveFigData = 0;            
figStruct.figPos = [0 0 .6 0.7];      % Normalized units for figure position

if ~exist(outDir, 'dir')
  mkdir(outDir);
end

% Add in combo events for later counting

% Make Grid hole arrays
gridHoles = unique(selTable.gridHole);
gridHoles = gridHoles(2:end);
gridHoleA = unique(string(extractBefore(gridHoles, 2)));
gridHoleB = string(unique(double(string(extractAfter(gridHoles, 1)))));
[gridRuns, gridChannels, gridUnits]  = deal(zeros(length(gridHoleA), length(gridHoleB)));

% Get a unit count per grid hole
for grid_i = 1:length(gridHoles)
  % Identify spot to add
  aInd = strcmp(gridHoles{grid_i}(1), gridHoleA);
  bInd = strcmp(gridHoles{grid_i}(2:end), gridHoleB);
  
  % Identify runs
  runInds = strcmp(selTable.gridHole, gridHoles{grid_i});
  runsInGridHole = unique(runList(runInds));
  gridRuns(aInd, bInd) = gridRuns(aInd, bInd) + length(runsInGridHole);
  
  % Identify channel count
  for run_i = 1:length(runsInGridHole)
    gridChannels(aInd, bInd) = gridChannels(aInd, bInd) + spikePathBankParadigm{strcat('S', runsInGridHole{run_i}), "chanCount"};
  end
  
  % Identify Isolated Units
  for run_i = 1:length(runsInGridHole)
    unitInd = strcmp(runList, runsInGridHole{run_i}) & contains(selTable.unitType, digitsPattern);
    gridUnits(aInd, bInd) = gridUnits(aInd, bInd) + sum(unitInd);
  end
  
end

% Will be used to normalize counts, need to match unitTypes.
gridUnitCounts = cat(3, gridChannels, gridUnits);

% Recover the events represented
uniqueEvents = selTable.Properties.VariableNames(contains(selTable.Properties.VariableNames, 'Sel'));

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
  fprintf('\n')
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

h.Children(2).YTick = 1:3;
h.Children(2).YTickLabel = gridHoleA;
h.Children(2).XTick = 1:4;
h.Children(2).XTickLabel = gridHoleB;
title(figTitle)

saveFigure(outDir, figTitle, [], figStruct, []);

% Units per hole
paradigmPlotName = sprintf('Units per Grid Hole, %s Paradigm', barPlotParams.paradigm);
h = figure('Name',paradigmPlotName,'NumberTitle','off','units','normalized','outerposition',figStruct.figPos);
himg = imagesc(round(mean(gridUnits, 3)));
colorbar()

h.Children(2).YTick = 1:length(gridHoleA);
h.Children(2).YTickLabel = gridHoleA;
h.Children(2).XTick = 1:length(gridHoleB);
h.Children(2).XTickLabel = gridHoleB;
title(paradigmPlotName)
saveFigure(outDir, paradigmPlotName, [], figStruct, []);

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
      h = figure('Name', TopPlotName, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', figStruct.figPos);
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
      
      saveFigure(outDir, TopPlotName, [], figStruct, []);
      
    end

  end
end

disp('y')
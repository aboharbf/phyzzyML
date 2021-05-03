function selectivityCurveCountBinned(spikePathBank, batchAnalysisParams)
% a function which pulls encoding curves from the selectivity table for
% each paradigm.

runNames = extractAfter(spikePathBank.Properties.RowNames, 'S');
stretchAllSame = false;

% Collect unit selectivity
[anovaTablePerRun] = spikePathLoad(spikePathBank, {'anovaTable'}, batchAnalysisParams.spikePathLoadParams);

selTablePerRun = spikePathBank.selTable;

outputDir = fullfile(batchAnalysisParams.selParam.outputDir, 'selectivityCurveCount');
if ~exist(outputDir, 'dir')
  mkdir(outputDir);
end

paradigmList = unique(spikePathBank.paradigmName);
unitType = {'MUA'};
unitTypePlot = {'MUA', 'U+US'};

% Bins for the x-axis
binStep = 25;
binSize = 150;
stimSize = 2800 + 500;

the_bin_start_times = 0:binStep:stimSize-binSize;
points_to_label = [0  1000 2000 2800];
points_for_lines = [0 2800];

the_bin_start_times_shift = the_bin_start_times;
bins_to_label = interp1(the_bin_start_times, 1:length(the_bin_start_times), points_to_label);
x_for_lines = interp1(the_bin_start_times, 1:length(the_bin_start_times), points_for_lines);
labelOrder = {'chasing', 'fighting' 'grooming', 'mounting', 'goalDirected', 'idle', 'objects'};
alpha = 0.05;

for par_i = 1:length(paradigmList)
  
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  anovaTableParadigmPerRun = anovaTablePerRun(pInd);
  selTablePerRunParadigm = selTablePerRun(pInd);
  
  anovaTableParadigm = vertcat(anovaTableParadigmPerRun{:});
  uniqueObjects = unique(vertcat(anovaTableParadigm.slidingWin_broadCategoryTest_prefStim{:}));
  uniqueObjects = uniqueObjects(~contains(uniqueObjects, 'None'));
  selTablePerRunParadigm = vertcat(selTablePerRunParadigm{:});
  tableVars = selTablePerRunParadigm.Properties.VariableNames;
  
  % Find rows w/ Var
  swVars = tableVars(contains(tableVars, 'SW') & contains(tableVars, 'sigObj'))';
  swVarsParts = cellfun(@(x) strsplit(x, '_'), swVars, 'UniformOutput', false);
  testsPerformed = unique(cellfun(@(x) x(2), swVarsParts));
  epochSequence = {'stimOnset (0-500 ms)', 'stimPres (500-2800 ms)', 'reward (2800-3050 ms)'};
  
  for test_i = 1:length(testsPerformed)
    swVars = tableVars(contains(tableVars, 'SW') & contains(tableVars, testsPerformed{test_i}) & contains(tableVars, 'sigObj'))';
    swVarsParts = cellfun(@(x) strsplit(x, '_'), swVars, 'UniformOutput', false);
    epochs = cellfun(@(x) x(3), swVarsParts);
    
    % Turn the pVal vector into counts.
    for unitType_i = 1:length(unitTypePlot)
      
      objectArrays = nan(length(swVars), length(uniqueObjects));      
      for epoch_i = 1:length(epochs)
        
        % Grab one of the unit types
        if unitType_i == 1
          unitInd = contains(selTablePerRunParadigm.unitType, 'MUA');
        else
          unitInd = ~contains(selTablePerRunParadigm.unitType, 'MUA');
        end
        
        unitCount = sum(unitInd);
        
        % Extract and contactonate the indicies pointing to objects.
        epochObjData = selTablePerRunParadigm.(swVars{epoch_i})(unitInd);
        epochObjStack = vertcat(epochObjData{:});
        
        objectArrays(epoch_i, :) = sum(epochObjStack);
        
        if epoch_i == 1
          atleastOne = any(epochObjStack,2);
        else
          atleastOne = atleastOne | any(epochObjStack,2);
        end

      end
      
      % Plot as Bar graph.
      figTitle = sprintf('Significant %s counts per Category, Sliding Window Test, %s (%d of %d unique)', unitTypePlot{unitType_i}, paradigmList{par_i}, sum(atleastOne), unitCount);
      figH = createBarPlotWithChanceLine(epochSequence, objectArrays', 0, unitCount, figTitle, uniqueObjects);
      figH.Position = [0    0.0370    1.0000    0.8917];
      saveFigure(outputDir, figTitle, [], batchAnalysisParams.figStruct, [])

    end
        
  end
end

end
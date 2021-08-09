function selectivityCurveCountBinned(spikePathBank, batchAnalysisParams)
% a function which pulls encoding curves from the selectivity table for
% each paradigm.

barplotParams = struct();
barplotParams.labels = {'socialInteraction', 'chasing', 'fighting', 'mounting', 'grooming',...
                        'goalDirected', 'idle', 'objects', 'scene'};
barplotParams.colors = [1 0 0; 1 0 0; .95 0 0; .9 0 0; .85 0 0;...
                        .4 0 0; 0.9 0.6 0.1; 0 0.5 0; 0.5 0.5 0.5];

% Collect unit selectivity
[anovaTablePerRun] = spikePathLoad(spikePathBank, {'anovaTable'}, batchAnalysisParams.spikePathLoadParams);

selTablePerRun = spikePathBank.selTable;

outputDir = fullfile(batchAnalysisParams.selParam.outputDir, 'selectivityCurveCount');
if ~exist(outputDir, 'dir')
  mkdir(outputDir);
end

paradigmList = unique(spikePathBank.paradigmName);
unitLabel = {'MUA', digitsPattern};
unitTypePlot = {'MUA', 'Units'};

% Bins for the x-axis
binStep = 25;
binSize = 150;
stimSize = 2800 + 500;

% Resort for later
BaseArrayResort = {{'chasing', 'fighting', 'mounting', 'grooming', 'goalDirected', 'idle', 'objects', 'scene'}, {'socialInteraction', 'goalDirected', 'idle', 'objects', 'scene'}};   
BaseArrayTestNames = {'categoriesTest', 'broadCatTest'};

for par_i = 1:length(paradigmList)
  
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  anovaTableParadigmPerRun = anovaTablePerRun(pInd);
  anovaTableParadigm = vertcat(anovaTableParadigmPerRun{:});
  anovaTableVars = anovaTableParadigm.Properties.VariableNames';

  selTablePerRunParadigm = selTablePerRun(pInd);
  selTablePerRunParadigm = vertcat(selTablePerRunParadigm{:});
  selTableVars = selTablePerRunParadigm.Properties.VariableNames';
  
  tableVarsPrefStim = anovaTableVars(contains(anovaTableVars, '_prefStim'));
  
  for pref_i = 1:length(tableVarsPrefStim)
    
    testName = extractBetween(tableVarsPrefStim{pref_i}, '_', '_');
    uniqueObjects = unique(vertcat(anovaTableParadigm.(tableVarsPrefStim{pref_i}){:}));
    uniqueObjects = uniqueObjects(~contains(uniqueObjects, 'None'));

    swVars = selTableVars(contains(selTableVars, 'SW') & contains(selTableVars, testName{1}) & contains(selTableVars, 'sigObj'));
    swVarsParts = cellfun(@(x) strsplit(x, '_'), swVars, 'UniformOutput', false);
    epochs = cellfun(@(x) x(3), swVarsParts);
%     epochSequence = {'preFix (ITI)', 'Fixation (-800 - 0)', 'stimEarly (0-500 ms)', 'stimLate (500-2800 ms)', 'reward (2800-3050 ms)'};
%     epochSequence = {'preFix (ITI)', 'Fixation (-800 - 0)', 'Stimulus (0-2800 ms)', 'Reward (2800-3050 ms)'};
    epochSequence = {'Stimulus (0-2800 ms)'};
    
    % Turn the pVal vector into counts.
    for unitType_i = 1:length(unitLabel)
      
      objectArrays = nan(length(swVars), length(uniqueObjects));
      
      for epoch_i = 1:length(epochs)
        
        % Grab one of the unit types
        unitInd = contains(selTablePerRunParadigm.unitType, unitLabel{unitType_i});
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
      
      % Prior to plotting, resort
      testInd = strcmp(BaseArrayTestNames, testName);
      [~, resortInd] = ismember(BaseArrayResort{testInd}, uniqueObjects);
      objectArrays = objectArrays(:, resortInd);
      uniqueObjects = uniqueObjects(resortInd);
      
      % Plot as Bar graph.
      figTitle = sprintf('%s Significant %s counts per Category, Sliding Window Test (%d/%d unique)', testName{1}, unitTypePlot{unitType_i}, sum(atleastOne), unitCount);
      figH = createBarPlotWithChanceLine(epochSequence, objectArrays', 0, unitCount, figTitle, uniqueObjects, barplotParams);
      figH.Position = [.3    0.185    0.62    0.5];
      saveFigure(outputDir, sprintf('1.%s - %s', paradigmList{par_i}, strrep(figTitle, '/', ' of ')), [], batchAnalysisParams.figStruct, [])
      
    end
    
  end
end


end
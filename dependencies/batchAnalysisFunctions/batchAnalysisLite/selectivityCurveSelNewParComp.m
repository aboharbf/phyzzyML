function selectivityCurveSelNewParComp(spikePathBank, batchAnalysisParams)
% a function which pulls encoding curves from the selectivity table for
% each paradigm.

% Collect unit selectivity
[anovaTablePerRun, epochSWparamsPerRun] = spikePathLoad(spikePathBank, {'anovaTable', 'epochSWparams'}, batchAnalysisParams.spikePathLoadParams);
selTablePerRun = spikePathBank.selTable;

paradigmList = unique(spikePathBank.paradigmName);
unitLabel = {digitsPattern};
unitTypePlot = {'Units'};

% Bins for the x-axis
binStep = 25;
binSize = 150;
stimSize = 2800 + 500;
the_bin_start_times = 0:binStep:stimSize-binSize;
points_to_label = [0  1000 2000 2800];
points_for_lines = [0 2800];
barGraphData = nan(2,8);

for par_i = 1:length(paradigmList)
  
  outputDir = fullfile(batchAnalysisParams.selParam.outputDir, paradigmList{par_i}, 'selectivityCurveSel');
  
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});  
  epochSWParams = epochSWparamsPerRun{find(pInd, 1)};
  paradigmName = [lower(paradigmList{par_i}(1)) paradigmList{par_i}(2:end)];
  paradigmParams = epochSWParams.(paradigmName);
  testsPerformed = paradigmParams.testLabel(1:2);

  anovaTableParadigmPerRun = anovaTablePerRun(pInd);
  anovaTableParadigm = vertcat(anovaTableParadigmPerRun{:});
  tableVars = anovaTableParadigm.Properties.VariableNames';
  
  selTableParadigmPerRun = selTablePerRun(pInd);
  selTableParadigmPerRun = vertcat(selTableParadigmPerRun{:});
  tableVarsSel = selTableParadigmPerRun.Properties.VariableNames';
  
  % Cycle through tests
  for test_i = 1%:length(testsPerformed)
    
    testInd = contains(tableVars, testsPerformed{test_i}) & ~contains(tableVars, 'EyesTest');
    testObjs = paradigmParams.testCategoryLabels{test_i};
    
    % Find rows w/ Var - ANOVA Table
    pValLabels = tableVars(contains(tableVars, 'pVal') & testInd);
    varLabels = tableVars(contains(tableVars, 'VarExp') & testInd);
    prefSelLabels = tableVars(contains(tableVars, '_prefStim') & testInd);
    
    testInd = contains(tableVarsSel, testsPerformed{test_i});

    % Find rows w/ Var - sel Table
    selTableVars = tableVarsSel(contains(tableVarsSel, 'SW') & testInd);
    sigStretchInd = selTableVars(contains(selTableVars, 'sigStretch'));
    prefOjbInd = selTableVars(contains(selTableVars, 'sigObj'));
    
    dataArray = {pValLabels; varLabels; prefSelLabels};
    
    % Cycle through analyses, combining results and visualizing them
    objCount = length(selTableParadigmPerRun.(prefOjbInd{1}){1});
    %barGraphData = nan(length(sigStretchInd), objCount);
    
    % Find epoch labels
    tmp = cellfun(@(x) strsplit(x, '_'), sigStretchInd, 'UniformOutput', false);
    epochNames = cellfun(@(x) x(3), tmp);
    
    for unit_i = 1:length(unitLabel)
      
      unit_index = contains(selTableParadigmPerRun.unitType, unitLabel{unit_i});
      
      % Store data per epoch
      vennDiagramVectors = nan(sum(unit_index), length(epochNames));
      
      selStretchIndMat = nan(sum(unit_index), length(epochNames));
      for epoch_i = 1:length(epochNames)
        % Identify runs with significant stretches, and retrieve the labels
        % associated with them
        selStretchIndMat(:,epoch_i) = selTableParadigmPerRun.(sigStretchInd{epoch_i})(unit_index);
        
        % Significant objects
        sigObjs = selTableParadigmPerRun.(prefOjbInd{epoch_i})(unit_index);
        sigObjs = vertcat(sigObjs{:});
        barGraphData(par_i,:) = sum(sigObjs);
        
        % Store for VD
        vennDiagramVectors(:, epoch_i) = any(sigObjs,2);
        
      end
      
      atleastOne = sum(any(selStretchIndMat,2));
      unitCount(par_i) = sum(unit_index);

    end
    
  end
end

for i=1:size(testObjs,2)
    testObjs{i}(1) = upper(testObjs{i}(1));
end

% Plotting
legend2Add = strcat({'Natural', 'Animated'}, ' (', string(unitCount), ')');
figTitle = sprintf('%s selectivity during %s test, %s (%d/%d Unique)', unitTypePlot{unit_i}, testsPerformed{test_i}, paradigmList{par_i}, atleastOne, unitCount);
plotTitle = sprintf('Category Selectivity');
createBarPlotWithChanceLine(testObjs, barGraphData, 0, unitCount, figTitle, legend2Add);
title(plotTitle);

axesH = gca;
axesH.FontSize = 16;
axesH.Title.FontSize = 16;
axesH.XTickLabelRotation = 0;
axesH.XTickLabel{6} = 'Goal\newlineDirected';

figH = gcf;
figH.Position = [0.2531    0.2333    0.6035    0.4618];

saveFigure(fullfile(outputDir, 'BarGraph'), sprintf('%s - Paradigm Compare', plotTitle), [], batchAnalysisParams.figStruct, [])


end
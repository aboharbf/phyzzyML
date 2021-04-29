function selectivityCurveCount(spikePathBank, batchAnalysisParams)
% a function which pulls encoding curves from the selectivity table for
% each paradigm.

runNames = extractAfter(spikePathBank.Properties.RowNames, 'S');
stretchAllSame = false;

% Collect unit selectivity
[selTablePerRun] = spikePathLoad(spikePathBank, {'anovaTable'}, batchAnalysisParams.spikePathLoadParams);
outputDir = batchAnalysisParams.selParam.outputDir;

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
alpha = 0.01;

for par_i = 1:length(paradigmList)
  
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  selTableParadigmPerRun = selTablePerRun(pInd);
  
  selTableParadigm = vertcat(selTableParadigmPerRun{:});
  
  % Find rows w/ Var
  pValLabels = selTableParadigm.Properties.VariableNames(contains(selTableParadigm.Properties.VariableNames, 'pVal'));
  pValLabels = pValLabels(1);
  prefSelLabels = selTableParadigm.Properties.VariableNames(contains(selTableParadigm.Properties.VariableNames, 'prefSel'));
  
  % Turn the pVal vector into counts.
  for p_i = 1:length(pValLabels)
    for unitType_i = 1:length(unitTypePlot)
      
      % Grab one of the unit types
      if unitType_i == 1
        unitInd = contains(selTableParadigm.unitType, 'MUA');
      else
        unitInd = ~contains(selTableParadigm.unitType, 'MUA');
      end
      
      unitCount = sum(unitInd);
      
      unitDatapVal = selTableParadigm.(pValLabels{p_i})(unitInd,:);
      unitDataPrefStim = selTableParadigm.(prefSelLabels{p_i})(unitInd,:); 
      unitDatapVal = vertcat(unitDatapVal{:});
      unitDataPrefStim = vertcat(unitDataPrefStim{:});
      
      % Preferred stimuli need to be turned into numbers
      uniqueStimLabels = unique(unitDataPrefStim(:));
      [~, B, ~] = intersect(uniqueStimLabels, labelOrder);
      uniqueStimLabels = uniqueStimLabels(B);
      dataStackFinalTmp = zeros([size(unitDataPrefStim), length(uniqueStimLabels)]);
      for ii = 1:length(uniqueStimLabels)
        layerInd = strcmp(uniqueStimLabels{ii}, labelOrder);
        dataStackFinalTmp(:, :, layerInd) = strcmp(unitDataPrefStim, uniqueStimLabels{ii});
      end
      
%       unitDataPrefStim = sum(dataStackFinalTmp, 3);
      unitDataPrefStim = dataStackFinalTmp;
      
      % Find units with significant runs
      stretchRuns = findStretchRuns(unitDataPrefStim);
      
      % Sig count per label
      countPerLabel = sum(squeeze(stretchRuns) >= 5);
      
      % Plot as Bar graph.
      figTitle = sprintf('Significant %s counts per Category, %s', unitTypePlot{unitType_i}, paradigmList{par_i});
      createBarPlotWithChanceLine(labelOrder, countPerLabel, alpha^5, unitCount, figTitle)
      
    end
  end
  
end

end
function selectivityCurveSel(spikePathBank, batchAnalysisParams)
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

for par_i = 1:length(paradigmList)
  
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  selTableParadigmPerRun = selTablePerRun(pInd);
  
  selTableParadigm = vertcat(selTableParadigmPerRun{:});
  
  % Find rows w/ Var
  pValLabels = selTableParadigm.Properties.VariableNames(contains(selTableParadigm.Properties.VariableNames, 'pVal'));
  varLabels = selTableParadigm.Properties.VariableNames(contains(selTableParadigm.Properties.VariableNames, 'VarExp'));
  prefSelLabels = selTableParadigm.Properties.VariableNames(contains(selTableParadigm.Properties.VariableNames, 'prefSel'));
  
  dataArray = {pValLabels; varLabels; prefSelLabels};
  dataLabels = {'-log10(p) Encoding Strength', 'Variance Explained', 'Preferred Stimuli'};
  stretchIndMat = cell(2, length(dataArray{1}));
  
  % Cycle through analyses, combining results and visualizing them
  for label_i = 1:length(dataLabels)
    % Grab the columns which represent this analysis
    analysisLabels = strrep(dataArray{label_i}, '_', ' '); %selTableParadigm.Properties.VariableNames(contains(selTableParadigm.Properties.VariableNames, [anovaFactors{label_i} 'Var']));
    
    % For each factor of the analysis, plot a
    factorCount = length(dataArray{label_i});
    figTitle = sprintf('ANOVA on Label %s - Paradigm %s', dataLabels{label_i}, paradigmList{par_i});
    figH = figure('Name', figTitle, 'NumberTitle', 'off', 'units', 'normalized', 'position', [0.0214    0.2287    0.9615    0.5213]);
    sgtitle(figTitle);
    
    for factor_i = 1:factorCount
      for unit_i = 1:2
        
        % Grab one of the unit types
        if unit_i == 1
          unitInd = contains(selTableParadigm.unitType, 'MUA');
        else
          unitInd = ~contains(selTableParadigm.unitType, 'MUA');
        end
        
        % Extract values for unit
        dataStack = vertcat(selTableParadigm.(dataArray{label_i}{factor_i}){unitInd,:});
        
        if contains(dataLabels{label_i}, 'Encoding')
          
          pValInd2Keep = dataStack < 0.05;
          stretchLength = findStretchRuns(pValInd2Keep);
          
          % For units with the right length of stretches, visualize them
          [stretchIndMat{unit_i, factor_i}, ind2Keep] = deal(stretchLength >= 5);
          dataStackKeep = dataStack(ind2Keep,:);
          stretchLengthKeep = stretchLength(ind2Keep,:);
          [~, ai] = sort(stretchLengthKeep);
          dataStackKeep = dataStackKeep(ai,:);
          
          % convert to negative log 2 to 5
          dataStackFinal = -log10(dataStackKeep);
          dataStackFinal(dataStackFinal > 5) = 5;
          dataStackFinal(dataStackFinal < 2) = 2;
          
        elseif contains(dataLabels{label_i}, 'Preferred')
          
          % Preferred stimuli need to be turned into numbers
          uniqueStimLabels = unique(dataStack(:));
          dataStackFinalTmp = zeros([size(dataStack), length(uniqueStimLabels)]);
          for ii = 1:length(uniqueStimLabels)
            dataStackFinalTmp(:,:,ii) = strcmp(dataStack, uniqueStimLabels{ii}) * (ii - 1);
          end
          dataStack = sum(dataStackFinalTmp, 3);
          
          if stretchAllSame
            stretchLengths = findStretchRuns(dataStackFinalTmp);
            [ind2Keep, stretchIndMat{unit_i, factor_i}] = deal(logical(sum(squeeze(stretchLengths) >= 5, 2)));
          else
            ind2Keep = stretchIndMat{unit_i, factor_i};
          end
          
          dataStackFinal = dataStack(ind2Keep,:);
          
        else
          
          % Sort the explained variance in the same way the p value stacks
          % were sorted.
          ind2Keep = stretchIndMat{unit_i, factor_i};
          dataStackFinal = dataStack(ind2Keep,:);
          
        end
                
        % Plotting
        ax = subplot(2, factorCount, sub2ind([factorCount 2], factor_i, unit_i));
        imagesc(dataStackFinal)
        title(sprintf('%s, %s - %d/%d', unitTypePlot{unit_i}, analysisLabels{factor_i}, sum(ind2Keep), length(ind2Keep)))
        
        if label_i == 2 || label_i == 3
          cbHand = colorbar();
        end
        
        % Convert Bins to times, add to X axis.
        ax.XTick = bins_to_label;
        ax.XTickLabel = points_to_label;
        ylabel('Unit #');
        xlabel('Bin (ms)');
        
      end
    end
    
    % Save the Figure
    saveFigure(batchAnalysisParams.selParam.outputDir, figTitle, [], batchAnalysisParams.selParam.figStruct, [])
    
    % See if Units are together or not.
    figTitleSum = sprintf('Counts of Traces with Significant Stretches %s', paradigmList{par_i});
    figure('Name', figTitleSum, 'units', 'normalized', 'position', [0.0635    0.0380    0.6339    0.8833])
    unitLabel = {'Units', 'MUA'};     % Units
    for ii = 1:length(unitLabel)
      ax = subplot(1,length(unitLabel),ii);
      anyUnit = [stretchIndMat{ii, :}];
      unitCountPerLabel = sum(anyUnit);
      anyUnitCount = sum(anyUnit, 2) ~= 0;
      allUnit = [anyUnit, anyUnitCount];
      imagesc(allUnit);
      title(sprintf('%s (Total %d/%d)', unitLabel{ii}, sum(anyUnitCount), length(anyUnitCount)))
      
      % Generate Labels
      analysisLabelsPlot = analysisLabels;
      for jj = 1:length(analysisLabels)
        analysisLabelsPlot{jj} = sprintf('%s (%d)', analysisLabels{jj}, unitCountPerLabel(jj));
      end
      xticks(1:length(analysisLabelsPlot)+1)
      xticklabels([analysisLabelsPlot, {'All Units'}])
      ax.XTickLabelRotation = 45;
    end
    sgtitle(figTitleSum);
    ylabel('Unit #')
    xlabel('ANOVA Factor');
    
    % Save the paneled image.
    saveFigure(batchAnalysisParams.selParam.outputDir, figTitleSum, [], batchAnalysisParams.selParam.figStruct, [])
    
  end
  
end

end
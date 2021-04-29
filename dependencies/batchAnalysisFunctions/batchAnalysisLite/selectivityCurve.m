function selectivityCurve(spikePathBank, batchAnalysisParams)
% a function which pulls encoding curves from the selectivity table for
% each paradigm. 

runNames = extractAfter(spikePathBank.Properties.RowNames, 'S');

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

for par_i = 2:3 %length(paradigmList)
  
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  spikePathBankParadigm = spikePathBank(pInd,:);
  runNamesParadigm = runNames(pInd,:);
  selTableParadigmPerRun = selTablePerRun(pInd);
  
  selTableParadigm = vertcat(selTableParadigmPerRun{:});

  % Find rows w/ Var
  pValLabels = selTableParadigm.Properties.VariableNames(contains(selTableParadigm.Properties.VariableNames, 'pVal'));
  varLabels = selTableParadigm.Properties.VariableNames(contains(selTableParadigm.Properties.VariableNames, 'VarExp'));
  dataArray = {pValLabels; varLabels};
  dataLabels = {'-log10(p) Encoding Strength', 'Variance Explained'};
  stretchIndMat = cell(2, length(dataArray{1}));
  
  % Cycle through analyses, combining results and visualizing them
  for label_i = 1:2
    % Grab the columns which represent this analysis
    analysisLabels = strrep(dataArray{label_i}, '_', ' '); %selTableParadigm.Properties.VariableNames(contains(selTableParadigm.Properties.VariableNames, [anovaFactors{label_i} 'Var']));
    analysisLabels(4:5) = {'Interaction 1 & 3', 'Interaction 2 & 3'};
    
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
          
          % Clean up the pValStack by adding in a stretch of significance
          % filter. Do this by adding a pad of zeros, then running diff, and
          % removing the pads
          
          zPad = zeros(size(dataStack,1),1);
          pValInd2Keep = dataStack < 0.05;
          pValStackDiff = diff([zPad, pValInd2Keep, zPad], 1, 2);
          stretchLength = zeros(size(pValStackDiff,1),1);
          for row_i = 1:length(stretchLength)
            starts = find(pValStackDiff(row_i, :) == 1);
            stops =  find(pValStackDiff(row_i, :) == -1);
            lengths = stops - starts;
            
            if ~isempty(lengths)
              stretchLength(row_i) = max(lengths);
            end
            
          end
          
          % For units with the right length of stretches, visualize them
          [stretchIndMat{unit_i, factor_i}, ind2Keep] = deal(stretchLength >= 6);
          dataStackKeep = dataStack(ind2Keep,:);
          stretchLengthKeep = stretchLength(ind2Keep,:);
          [~, ai] = sort(stretchLengthKeep);
          dataStackKeep = dataStackKeep(ai,:);
          
          % convert to negative log 2 to 5
          dataStackFinal = -log10(dataStackKeep);
          dataStackFinal(dataStackFinal > 5) = 5;
          dataStackFinal(dataStackFinal < 2) = 2;
          
        else
          
          % Sort the explained variance in the same way the p value stacks
          % were sorted.
          ind2Keep = stretchIndMat{unit_i, factor_i};
          dataStackFinal = dataStack(ind2Keep,:);
          
        end
        
        sigCount = sum(ind2Keep);
        
        % Plotting
        ax = subplot(2, factorCount, sub2ind([factorCount 2], factor_i, unit_i));
        imagesc(dataStackFinal)
        title(sprintf('%s, %s - %d/%d', unitTypePlot{unit_i}, analysisLabels{factor_i}, sigCount, length(ind2Keep)))
        if label_i == 2
          colorbar()
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
      xticklabels([analysisLabelsPlot, {'All Units'}])
      ax.XTickLabelRotation = 45;
    end
    sgtitle(figTitleSum);
    ylabel('Unit #')
    xlabel('ANOVA Factor');
    
    saveFigure(batchAnalysisParams.selParam.outputDir, figTitleSum, [], batchAnalysisParams.selParam.figStruct, [])  
    
  end
  
end



end
function selCount(spikePathBank, batchAnalysisParams)
% A Function to tally across processed runs and the 'selTable' produced.

[selTablePerRun] = spikePathLoad(spikePathBank, {'selTable'}, batchAnalysisParams.spikePathLoadParams);
outputDir = batchAnalysisParams.selParam.outputDir;

if ~exist(outputDir, 'dir')
  mkdir(outputDir);
end

paradigmList = unique(spikePathBank.paradigmName);
selCheck = {'socInt', 'headTurn', 'fullModel', 'turnToward'};

UnitTypes = {'MUA', digitsPattern};
UnitTypePlot = {'MUA', 'Units'};
colNamePoss = {'stimOnset','stimPres','reward'};
colNamePlotAll = {'stim Onset','stim Presentation','Reward'};

for par_i = 2%:length(paradigmList)
  
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  spikePathBankParadigm = spikePathBank(pInd,:);
  selTableParadigmPerRun = selTablePerRun(pInd);
     
  % Combine across tables
  selTableParadigm = vertcat(selTableParadigmPerRun{:});
  
  for unitType_i = 1:length(UnitTypes)
    % Filter the table per unit
    unitInd = contains(selTableParadigm.unitType, UnitTypes{unitType_i});
    selTableParadigmUnit = selTableParadigm(unitInd, :);
    
    % Plots 1 - fixation/saccade
    unitCount = sum(unitInd);
    saccSelCount = sum(selTableParadigmUnit.saccSel);
    fixSelCount = sum(selTableParadigmUnit.fixationSel ~= 0);
    
    for sel_i = 1:length(selCheck)
      % Identify which columns were part of this analysis
      columns2Combine = contains(selTableParadigmUnit.Properties.VariableNames, [selCheck{sel_i} 'Sel_']);
      
      if any(columns2Combine)
        % Assuming each only has 3 possible comparisons (onset, pres,
        % reward at most).
        
        % Create the groups + their intersections
        selTableParadigmUnitCat = selTableParadigmUnit{:, columns2Combine};
        selTableParadigmVarNames = selTableParadigmUnit.Properties.VariableNames(columns2Combine);
        colNames = extractAfter(selTableParadigmVarNames, [selCheck{sel_i} 'Sel_']);
        [~, A] = intersect(colNamePoss, colNames);
        A = sort(A);
        colNamesPlot = colNamePlotAll(A);
        
        % make the names easier
        for col_i = 1:length(colNames)
          eval(sprintf('%s = selTableParadigmUnit.%s;', colNames{col_i}, selTableParadigmVarNames{col_i}));
        end

        % Generate Data structure
        switch length(colNames)
          case 3
            dataMat = [[sum(stimOnset > 0), sum(stimOnset < 0)];...
              [sum(stimPres > 0), sum(stimPres < 0)];...
              [sum(reward > 0), sum(reward < 0)]];
            
            atleastOne = sum(stimOnset | stimPres | reward);
          case 2
            dataMat = [[sum(stimPres > 0), sum(stimPres < 0)];...
              [sum(reward > 0), sum(reward < 0)]];
            
            atleastOne = sum(stimPres | reward);

        end
        
        % Plot
        figTitle = sprintf('%s activity selective for %s during %s (%d Unique)', UnitTypePlot{unitType_i}, selCheck{sel_i}, paradigmList{par_i}, atleastOne);
        figH = figure('Name', figTitle, 'NumberTitle','off','units','normalized', 'outerposition', [.3 .3 .4 .6]);
        X = categorical(colNamesPlot);
        X = reordercats(X,colNamesPlot);
        barh = bar(X, dataMat, 'stacked');
        for bar_i = 1:length(barh)
          xtips1 = barh(bar_i).XEndPoints;
          ytips1 = barh(bar_i).YEndPoints;
          labels1 = string(barh(bar_i).YData);
          labels2 = string(round(barh(bar_i).YData/unitCount,3));
          labels3 = strcat(labels1, '(', labels2, ')');
          
          text(xtips1, ytips1, labels3, 'HorizontalAlignment', 'center', 'VerticalAlignment','top')
        end
        YlimN = figH.Children.YLim;
        figH.Children.YLim(2) = YlimN(2) * 1.2;
        legend('Increased firing', 'decreased firing');
        title(figTitle);
        
        batchAnalysisParams.figStruct
        saveFigure(outputDir, figTitle, [], batchAnalysisParams.figStruct, [])
        
      end
    end
    
  end
  
end


end
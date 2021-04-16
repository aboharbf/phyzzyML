function selectivityPerEpochBarGraphs(selTableParadigm, paradigm, barPlotParams)

selCheck = barPlotParams.selCheck;
UnitTypes = barPlotParams.UnitTypes;
UnitTypePlot = barPlotParams.UnitTypePlot;
colNamePoss = barPlotParams.colNamePoss;
colNamePlotAll = barPlotParams.colNamePlotAll;

% Generate bar plots showing total counts in each category
for unitType_i = 1:length(UnitTypes)
  % Filter the table per unit
  unitInd = contains(selTableParadigm.unitType, UnitTypes{unitType_i});
  selTableParadigmUnit = selTableParadigm(unitInd, :);
  
  if ~isempty(selTableParadigmUnit)
    % Plots 1 - fixation/saccade
    unitCount = sum(unitInd);
    chanceUnitCount = round(sum(unitCount)*0.025);
    
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
        colNames = colNames(A);
        colNamesPlot = colNamePlotAll(A);
        
        colNamesPlot = colNamesPlot(~contains(colNames' ,'stimWhole'));
        colNames = colNames(~contains(colNames' ,'stimWhole'));
        
        % make the names easier
        for col_i = 1:length(colNames)
          eval(sprintf('%s = selTableParadigmUnit.%s;', colNames{col_i}, selTableParadigmVarNames{col_i}));
          eval(sprintf('%s(isnan(%s)) = 0;', colNames{col_i}, colNames{col_i}));
        end
        
        % Get rid of redundant elements
        
        
        % Generate Data structure
        switch length(colNames)
          case 3
            dataMat = [[sum(stimOnset > 0), sum(stimOnset < 0)];...
              [sum(stimPres > 0), sum(stimPres < 0)];...
              [sum(reward > 0), sum(reward < 0)]];
            
            atleastOne = sum(stimOnset | stimPres | reward);
            
%             dataMat = [[sum(stimOnset > 0), sum(stimOnset < 0)];...
%               [sum(stimPres > 0), sum(stimPres < 0)]];
%             colNamesPlot = colNamesPlot(1:2);
%             
%             atleastOne = sum(stimOnset | stimPres);
          case 2
            dataMat = [[sum(stimPres > 0), sum(stimPres < 0)];...
              [sum(reward > 0), sum(reward < 0)]];
            
            atleastOne = sum(stimPres | reward);
            
        end
        
        % Plot
        figTitle = sprintf('%s activity selective for %s during %s (%d Unique)', UnitTypePlot{unitType_i}, selCheck{sel_i}, paradigm, atleastOne);
        figH = figure('Name', figTitle, 'NumberTitle','off','units','normalized', 'outerposition', [.3 .3 .4 .6]);
        
        %       dataMat = dataMat(~strcmp(colNamesPlot, 'stim Whole Presentation'), :);
        %       colNamesPlot = colNamesPlot(~strcmp(colNamesPlot, 'stim Whole Presentation'));
        
        X = categorical(colNamesPlot);
        X = reordercats(X,colNamesPlot);
        barh = bar(X, dataMat);
        for bar_i = 1:length(barh)
          xtips1 = barh(bar_i).XEndPoints;
          ytips1 = barh(bar_i).YEndPoints;
          labels1 = string(barh(bar_i).YData);
          labels2 = string(round(barh(bar_i).YData/unitCount,3));
          labels3 = strcat(labels1, '(', labels2, ')');
          
          text(xtips1, ytips1, labels3, 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom')
        end
        YlimN = figH.Children.YLim;
        figH.Children.YLim(2) = YlimN(2) * 1.2;
        legend('Increased firing', 'decreased firing', 'Location','northeastoutside');
        xlim(xlim());
        hold on
        chanceLine = plot([X(1) X(end)], [chanceUnitCount chanceUnitCount], 'Color', 'red', 'LineStyle', '--', 'linewidth', 3);
        chanceLine.DisplayName = 'Chance (2.5%)';
        figH.Children(2).FontSize = 16;
        xlabel('Epoch Compared')
        ylabel('Unit Count')
        
        %       title(figTitle);
        
        saveFigure(barPlotParams.outputDir, figTitle, [], barPlotParams.figStruct, [])
        
      end
    end
  end
  
end
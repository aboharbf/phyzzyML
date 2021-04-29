function selectivityPerEpochVennDiagram(selTableParadigm, paradigm, barPlotParams)

selCheck = barPlotParams.(paradigm).selCheck;
UnitTypes = barPlotParams.UnitTypes;
UnitTypePlot = barPlotParams.UnitTypePlot;
colNamePoss = barPlotParams.colNamePoss;
alpha = 0.01;

% Generate bar plots showing total counts in each category
for unitType_i = 1:length(UnitTypes)
  % Filter the table per unit
  unitInd = contains(selTableParadigm.unitType, UnitTypes{unitType_i});
  selTableParadigmUnit = selTableParadigm(unitInd, :);
  
  if ~isempty(selTableParadigmUnit)    
    for sel_i = 1:length(selCheck)
      % Identify which columns were part of this analysis
      columns2Combine = contains(selTableParadigmUnit.Properties.VariableNames, [selCheck{sel_i} 'Sel_']);
      
      if any(columns2Combine)
        
        % Create the groups + their intersections
        selTableParadigmVarNames = selTableParadigmUnit.Properties.VariableNames(columns2Combine);
        colNames = extractAfter(selTableParadigmVarNames, [selCheck{sel_i} 'Sel_']);
        %         colNames = colNames(~(strcmp(colNames, 'stimWhole') | strcmp(colNames, 'any'))); % Don't include these.
        [~, A] = intersect(colNames, colNamePoss);
        A = sort(A);
        colNames = colNames(A);
        selTableParadigmVarNames = selTableParadigmVarNames(A);
        sigInd = selTableParadigmUnit{:, selTableParadigmVarNames};
        % Get rid of NaNs
        if iscell(sigInd)
          
          % Find unique labels and create a 3rd dimension for each.
          uniqueLabels = unique(sigInd(:));
          sigIndTmp = zeros([size(sigInd), length(uniqueLabels)]);
          for lab_i = 1:length(uniqueLabels)
            sigIndTmp(:, :, lab_i) = strcmp(sigInd, uniqueLabels(lab_i));
          end
          
          % Get rid of 'None'.
          sigInd = sigIndTmp(:, :, ~strcmp(uniqueLabels, 'None'));
          uniqueLabels = uniqueLabels(~strcmp(uniqueLabels, 'None'));
          
          % Collapse across objects
          sigInd = any(sigInd, 3);
                    
        else
          sigInd(isnan(sigInd)) = 0;
          sigInd = logical(sigInd);
        end
        
        atleastOne = sum(any(sigInd, 2));               
        figTitle = sprintf('%s activity selective for %s during %s (%d Unique)', UnitTypePlot{unitType_i}, selCheck{sel_i}, paradigm, atleastOne);
        vennXExpanded(sigInd, figTitle, colNames)
        saveFigure(barPlotParams.outputDir, figTitle, [], barPlotParams.figStruct, [])
      end
      
    end
  end
  
end
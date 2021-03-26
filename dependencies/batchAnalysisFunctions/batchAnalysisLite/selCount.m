function selCount(spikePathBank, batchAnalysisParams)
% A Function to tally across processed runs and the 'selTable' produced.

[selTablePerRun] = spikePathLoad(spikePathBank, {'selTable'}, batchAnalysisParams.spikePathLoadParams);

paradigmList = unique(spikePathBank.paradigmName);
selCheck = {'socialInteraction', 'headTurn', 'fullModel', 'turnToward'};
UnitTypes = {'MUA', 'US', digitsPattern};
UnitTypePlot = {'MUA', 'Unsorted', 'Units'};

for par_i = 2:length(paradigmList)
  
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  spikePathBankParadigm = spikePathBank(pInd,:);
  selTableParadigmPerRun = selTablePerRun(pInd);
    
  % Add a runID column to each table
  runListParadigm = spikePathBankParadigm.Properties.RowNames;
  for run_i = 1:length(runListParadigm)
    rowCount = size(selTableParadigmPerRun{run_i}, 1);
    selTableParadigmPerRun{run_i}.runID = repmat(string(runListParadigm{run_i}), [rowCount 1]);
  end
  
  % Combine across tables
  selTableParadigm = vertcat(selTableParadigmPerRun{:});
  
  for unitType_i = 1:length(UnitTypes)
    % Filter the table per unit
    unitInd = contains(selTableParadigm.unitTypeVec, UnitTypes{unitType_i});
    selTableParadigmUnit = selTableParadigm(unitInd, :);
    unitCount = sum(unitInd);
    
    % Plots 1 - fixation/saccade
    saccSelCount = sum(selTableParadigmUnit.saccSel);
    fixSelCount = sum(selTableParadigmUnit.fixationSel ~= 0);
    
    for sel_i = 1:length(selCheck)
      % Identify which columns were part of this analysis
      columns2Combine = contains(selTableParadigmUnit.Properties.VariableNames, selCheck{sel_i});
      
      if any(columns2Combine)
        % Assuming each only has 3 possible comparisons (onset, pres, rew at
        % most).
        
        % Create the groups + their intersections
        selTableParadigmUnitCat = selTableParadigmUnit{:, columns2Combine};
        selTableParadigmVarNames = selTableParadigmUnit.Properties.VariableNames(columns2Combine);
        colNames = extractAfter(selTableParadigmVarNames, [selCheck{sel_i} 'Sel_']);
        
        % make the names easier
        for col_i = 1:length(colNames)
          eval(sprintf('%s = selTableParadigmUnit.%s ~= 0;', colNames{col_i}, selTableParadigmVarNames{col_i}));
        end
        
        switch length(colNames)
          case 3
            A = stimOnset;
            B = stimPres;
            C = reward;
            
            ABC = sum(A & B & C);
            AB = sum(A & B) - ABC;
            AC = sum(A & C) - ABC;
            BC = sum(B & C) - ABC;
            
            A = sum(A) - AB - AC - ABC;
            B = sum(B) - AB - BC - ABC;
            C = sum(C) - AC - BC - ABC;
            
            % Plot for counts
            vennMat = [A AB B BC C AC ABC];            
          case 2
            A = stimPres;
            B = reward;
            
            AB = sum(A & B);
            
            A = sum(A) - AB;
            B = sum(B) - AB;
            vennMat = [A AB B];            
        end
        
        % Plotting
        vennX(vennMat, 0.01)
        title(sprintf('Selectivity for %s in Paradigm %s, Counts of %s', selCheck{sel_i}, paradigmList{par_i}, UnitTypePlot{unitType_i}))
        
        % Plot for fractions
        vennX(vennMat/sum(unitInd), 0.01)
        title(sprintf('Selectivity for %s in Paradigm %s, Fraction of %s', selCheck{sel_i}, paradigmList{par_i}, UnitTypePlot{unitType_i}))
        
        
      end
    end
    
  end
  
end


end
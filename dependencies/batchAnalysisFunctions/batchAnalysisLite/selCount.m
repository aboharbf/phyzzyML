function selCount(spikePathBank, batchAnalysisParams)
% A Function to tally across processed runs and the 'selTable' produced.
runNames = extractAfter(spikePathBank.Properties.RowNames, 'S');

% Collect unit selectivity
[selTablePerRun] = spikePathLoad(spikePathBank, {'selTable'}, batchAnalysisParams.spikePathLoadParams);
outputDir = batchAnalysisParams.selParam.outputDir;

if ~exist(outputDir, 'dir')
  mkdir(outputDir);
end

paradigmList = unique(spikePathBank.paradigmName);

for par_i = 2%:length(paradigmList)
  
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  spikePathBankParadigm = spikePathBank(pInd,:);
  runNamesParadigm = runNames(pInd,:);
  selTableParadigmPerRun = selTablePerRun(pInd);
  
  % Combine across tables
  selTableParadigm = vertcat(selTableParadigmPerRun{:});
  
  % Replace the channel names for 20201123Mo
  selTableParadigm = replaceChanNum_1123(selTableParadigm);
  

  % Run the function
%   selectivityPerEpochBarGraphs(selTableParadigm, paradigmList{par_i}, batchAnalysisParams.selParam)
  
  % Generate an array of gridHoles for each unit
  paradigmRuns = extractAfter(spikePathBankParadigm.Properties.RowNames, 'S');
  
  [chan, runPerChan] = deal([]); %cell(sum(chanPerRun), 1);
  gridHoleAll = [];
  for run_i = 1:length(paradigmRuns)
    % Determine shape and distribution of gridHoles
    gridHoleRun = strsplit(spikePathBankParadigm.GridHole{run_i}, ';');
    chNum = spikePathBankParadigm.chanCount(run_i);
    
    % Create reference structures
    channels2Name = strcat("Ch", string(1:chNum))';
    runPerChan = [runPerChan; repmat(string(paradigmRuns{run_i}), [chNum, 1])];
    chan = [chan; channels2Name];
    
    try
      switch chNum
        case 32 % One gridHole for all of them
          gridHoleChan = repmat(gridHoleRun, [chNum, 1]);
        case 64
          gridHoleChan = [repmat(gridHoleRun(1), [chNum/2, 1]); repmat(gridHoleRun(2), [chNum/2, 1])];
        case 33
          gridHoleChan = [repmat(gridHoleRun(1), [32, 1]); gridHoleRun(2)];
        otherwise
          if chNum == 1 && isempty(gridHoleRun{1})
            continue
          end
          
          gridHoleChan = [];
          for ii = 1:chNum
            gridHoleChan  = [gridHoleChan; gridHoleRun(ii)];
          end
          
      end
    end
    
    gridHoleAll = [gridHoleAll; gridHoleChan];
    
  end
  
  % Now use that matrix to find the gridHolePerUnit
  gridHolePerUnit = cell(size(selTableParadigm,1),1);
  fullRunName = strcat(selTableParadigm.dateSubj, selTableParadigm.runNum);
  for runChan_i = 1:length(runPerChan)
    popInd = strcmp(fullRunName, runPerChan(runChan_i)) & strcmp(selTableParadigm.channel, chan(runChan_i));
    gridHolePerUnit(popInd) = deal(gridHoleAll(runChan_i));
  end
  
  selTableParadigm.gridHole = gridHolePerUnit;
  batchAnalysisParams.selParam.paradigm = paradigmList{par_i};
  
  % Now, for each unitType and selectivity, map out 
  selectivityPerGridHole(spikePathBankParadigm, batchAnalysisParams.selParam, selTableParadigm)
  
end


end
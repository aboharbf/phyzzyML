function spikeEyeCorrHist(spikePathBank, batchAnalysisParams)
% A Function to tally across processed runs and the 'selTable' produced.

% Now, selectivity across paradigms...
% selCountCrossParadigm(spikePathBank, selTablePerRun, batchAnalysisParams);

paradigmList = unique(spikePathBank.paradigmName);
monkeyList = {'Sam', 'Mo', 'Combo'};
unitList = {'MU', 'U'};
unitTags = {'MUA', digitsPattern};
selParams = batchAnalysisParams.selParam;
mainOutDir = selParams.outputDir;
selParams.unitList = unitList;
selParams.unitTags = unitTags;

for m_i = 3%1:length(monkeyList)
  for par_i = 1:length(paradigmList)
    
    selParams.paradigmTag = paradigmList{par_i};
    selParams.monkeyTag = monkeyList{m_i};
    
    % Filter at the stage of spikePathBank - monkey and paradigm
    pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
    if ~strcmp(monkeyList{m_i}, 'Combo')
      mInd = contains(spikePathBank.Row, monkeyList{m_i});
    else
      mInd = true(size(pInd));
    end
    
    selTableParadigmPerRun = spikePathBank.selTable(pInd & mInd);
    
    % Combine across tables
    selTableParadigm = vertcat(selTableParadigmPerRun{:});
    
    % Set output directory
    selParams.outputDir = fullfile(mainOutDir, paradigmList{par_i}, monkeyList{m_i});
    
    % Cycle through cell types, create the plots
    for unit_i = 2%1:length(unitList)
      
      % Set output directory
      selParams.outputDir = fullfile(mainOutDir, paradigmList{par_i}, monkeyList{m_i}, unitList{unit_i});
      
      unitIndex = contains(selTableParadigm.unitType, unitTags{unit_i});
      selTableParadigmUnit = selTableParadigm(unitIndex, :);
      selParams.unitTag = unitList{unit_i};
            
      % Find units which specifically are highly orbit position dependent
      sigSpikeEyeCorrInd = find(selTableParadigmUnit.spikeEyeCorr >= prctile(selTableParadigmUnit.spikeEyeCorrNull, 95));
      spikeEyeCorrDiffSig = selTableParadigmUnit(sigSpikeEyeCorrInd, :);
      [~, sortInd] = sort(spikeEyeCorrDiffSig, 'descend');
      
      % Simple plot
%       sigUnits = selTableParadigmUnit(sigSpikeEyeCorrInd(sortInd),:);
      sigUnits = selTableParadigmUnit;
      figH = figure('units', 'normalized', 'position', [0.1922 0.3481 0.5995 0.3278]); 
      subplot(1,2,1)
      realHist = histogram(sigUnits.spikeEyeCorr, 20); 
      hold on; 
      nullHist = histogram(sigUnits.spikeEyeCorrNull, 20);
      nullHist.BinEdges = realHist.BinEdges;
      title('Split Halves Correlation, Real & Scramble');
      xlabel('Split halves correlation')
      ylabel('Count')
      legend({'Real Spike Times', 'Scrambled Spike Times'});
      figH.Children(2).FontSize = 12;
      
      subplot(1,2,2)
      compHist = histogram(sigUnits.spikeEyeCorr - sigUnits.spikeEyeCorrNull, 50);
      title('Split Halves Correlation, Real - Scramble');
      figH.Children(1).FontSize = 12;
      xlabel('Correlation Difference')
      ylabel('Count');
      
      
    end
  end
end

end
function socSaccAnalysis(spikePathBank, batchAnalysisParams)
% A Function to tally across processed runs and the 'selTable' produced.

% Now, selectivity across paradigms...
% selCountCrossParadigm(spikePathBank, selTablePerRun, batchAnalysisParams);

paradigmList = unique(spikePathBank.paradigmName);
monkeyList = {'Sam', 'Mo', 'Combo'};
unitList = {'U'};
unitTags = {digitsPattern};
selParams = batchAnalysisParams.selParam;
mainOutDir = selParams.outputDir;
selParams.unitList = unitList;
selParams.unitTags = unitTags;

for m_i = 3%1:length(monkeyList)
  for par_i = 1%:length(paradigmList)
    
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
    for unitType_i = 1:length(unitList)
      
      % Set output directory
      selParams.outputDir = fullfile(mainOutDir, paradigmList{par_i}, monkeyList{m_i}, unitList{unitType_i});
      
      unitIndex = contains(selTableParadigm.unitType, unitTags{unitType_i});
      selTableParadigmUnit = selTableParadigm(unitIndex, :);
      selParams.unitTag = unitList{unitType_i};
      
      selTableVars = selTableParadigmUnit.Properties.VariableNames';
      saccSelInd = contains(selTableVars, 'saccSel') & contains(selTableVars, 'selInd');
      subEventVars = selTableVars(saccSelInd);
      subEventName= extractBetween(subEventVars, 'saccSel_', '_selInd');
      subEventNamePlot = strrep(subEventName, '_', ' ');
       
      saccSelUnits = selTableParadigmUnit(selTableParadigmUnit{:, subEventVars(1)}, :);
      
      % Resort according to cohensD
      [~, resortInd] = sort(saccSelUnits.saccSel_socVNonSoc_diff, 'descend');
      saccSelUnits = saccSelUnits(resortInd, :);
      
      %
      openUnitFigs(saccSelUnits(1:10, :), 'imPSTH*', [], batchAnalysisParams.analysisDirectory)

      
      % Pre-saccade overlap vs Unit
      subEventMat = any(selTableParadigmUnit{:, subEventVars(1:4)}, 2);
      preSaccMat = any(selTableParadigmUnit{:, 'subSel_pre_saccades_selInd'}, 2);
      vennXExpanded([subEventMat, preSaccMat], 'subEvent vs Pre-saccade', {'subEvent Selective', 'Pre-Saccade Selective'})
      
    end
  end

end
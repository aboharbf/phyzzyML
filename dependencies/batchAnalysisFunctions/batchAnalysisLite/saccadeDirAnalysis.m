function saccadeDirAnalysis(spikePathBank, batchAnalysisParams)
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
    for unit_i = 1:length(unitList)
      
      % Set output directory
      selParams.outputDir = fullfile(mainOutDir, paradigmList{par_i}, monkeyList{m_i}, unitList{unit_i});
      
      unitIndex = contains(selTableParadigm.unitType, unitTags{unit_i});
      selTableParadigmUnit = selTableParadigm(unitIndex, :);
      selParams.unitTag = unitList{unit_i};
     
      % Pull 
      units2Find = 'saccGen';
      switch units2Find
        case 'saccDir'
          selTableSaccDir = selTableParadigmUnit(selTableParadigmUnit.saccDir_nonStim_selInd, :);
          unitScore = selTableSaccDir.saccDir_nonStim_pVal;
          figTag2Find = 'SaccDir_PSTH*';
          filtTag = [];
        case 'saccGen'
          %selTableSaccDir = selTableParadigmUnit(selTableParadigmUnit.subSel_pre_saccades_selInd, :);
          selTableSaccDir = selTableParadigmUnit(selTableParadigmUnit.subSel_saccades_selInd, :);
          unitScore = selTableSaccDir.subSel_pre_saccades_pVal;
          figTag2Find = 'subEvent_*';
          filtTag = 'saccade';
      end
      
      otherTag = 'subEventPSTH_*';
      
      unitLabels = strcat(selTableSaccDir.dateSubj, selTableSaccDir.runNum, '_', selTableSaccDir.channel, selTableSaccDir.unitType);

      % Sort
      [~, sortInd] = sort(unitScore);
      unitLabelsSorted = unitLabels(sortInd);
      
      % Presaccade
%       presaccadeUnits = {'20210625Sam005_Ch10U1', '20210629Sam003_Ch45U1', '202001204Mo003_Ch21U1','20201120Mo003_Ch2U1'};
      unitLabelsSortedTmp = unitLabelsSorted(contains(unitLabelsSorted, 'Ch8U1'));

      for ii = 1:20%length(unitLabelsSorted)
        %openUnitFigs(unitLabelsSorted(ii), figTag2Find, [], batchAnalysisParams.analysisDirectory)
        openUnitFigs(unitLabelsSorted(ii), otherTag, filtTag, batchAnalysisParams.analysisDirectory)
      end
      
      % Post saccade units
      for ii = 11:20%length(unitLabelsSorted)
        openUnitFigs(unitLabelsSortedTmp, figTag2Find, [], batchAnalysisParams.analysisDirectory)
        openUnitFigs(unitLabelsSortedTmp, otherTag, 'saccade', batchAnalysisParams.analysisDirectory)
      end
      
    end
  end
end

end
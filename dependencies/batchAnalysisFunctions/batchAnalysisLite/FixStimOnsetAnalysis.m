function FixStimOnsetAnalysis(spikePathBank, batchAnalysisParams)
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
      
      vennXExpanded([selTableParadigmUnit.subSel_fixation_selInd, selTableParadigmUnit.subSel_saccades_selInd], 'Stimulus Onset vs Pre-Saccade Selectivity', {'Stimulus Onset', 'Pre-Saccade Period'})
     
      % Read out Fix event vs Fix Epoch
      fixEvent = selTableParadigmUnit.subSel_fixation_selInd;
      fixEpoch = selTableParadigmUnit.baseV_Fix_selInd;
      vennXExpanded([fixEvent, fixEpoch], 'Title', {'Event', 'Epoch'});
      
      % Event not Epoch
      eventNotEpoch = find(fixEvent & ~fixEpoch);
      openUnitFigs(selTableParadigmUnit(eventNotEpoch(1:10),:), 'subEventPSTH_*', 'fix', [])
      openUnitFigs(selTableParadigmUnit(eventNotEpoch(1:10),:), 'imPSTH_*', [], [])
      
      % Pull 
      units2Find = 'stimOnset';
      switch units2Find
        case 'fixOnset'
          selTableSaccDir = selTableParadigmUnit(selTableParadigmUnit.subSel_fixation_selInd, :);
          unitScore = selTableSaccDir.subSel_fixation_pVal;
          figTag2Find = 'subEvent_*';
          filtTag = 'fix';
        case 'fixRamp'
          selTableSaccDir = selTableParadigmUnit(selTableParadigmUnit.subSel_fixation_selInd, :);
          unitScore = selTableSaccDir.subSel_fixation_pVal;
          figTag2Find = 'subEvent_*';
          filtTag = 'fix';
        case 'stimOnset'
          %selTableSaccDir = selTableParadigmUnit(selTableParadigmUnit.subSel_pre_saccades_selInd, :);
          selTableSaccDir = selTableParadigmUnit(selTableParadigmUnit.subSel_stimOnset_selInd, :);
          unitScore = selTableSaccDir.subSel_stimOnset_pVal;
          figTag2Find = 'subEvent_*';
          filtTag = 'stimOnset';
      end
      
      otherTag = 'subEventPSTH_*';
      unitLabels = strcat(selTableSaccDir.dateSubj, selTableSaccDir.runNum, '_', selTableSaccDir.channel, selTableSaccDir.unitType);
      [~, sortInd] = sort(unitScore);
      unitLabelsSorted = unitLabels(sortInd);
      
      % Fixation
      for ii = 1:10%length(unitLabelsSorted)
        openUnitFigs(unitLabelsSorted(ii), otherTag, filtTag, batchAnalysisParams.analysisDirectory)
        openUnitFigs(unitLabelsSorted(ii), 'imPSTH_*', [], batchAnalysisParams.analysisDirectory)
      end

      stimOnsetUnit = selTableParadigmUnit.subSel_stimOnset_selInd & selTableParadigmUnit.subSel_stimOffset_selInd;
      saccadeUnit = selTableParadigmUnit.subSel_pre_saccadesNonStim_selInd | selTableParadigmUnit.subSel_saccadesNonStim_selInd;
      fixOnsetUnit = selTableParadigmUnit.subSel_fixation_selInd;
      %vennXExpanded([selTableParadigmUnit.subSel_stimOnset_selInd, selTableParadigmUnit.subSel_stimOffset_selInd], 'Title', {'Stim Onset', 'Stim Offset'});
      vennXExpanded([saccadeUnit, fixOnsetUnit], 'Title', {'Pre/Post Saccade', 'Fix Onset'});
      
      selTableSaccDir = selTableParadigmUnit(stimOnsetUnit, :);
      unitLabels = strcat(selTableSaccDir.dateSubj, selTableSaccDir.runNum, '_', selTableSaccDir.channel, selTableSaccDir.unitType);

      openUnitFigs(unitLabels, 'imPSTH_*', [], batchAnalysisParams.analysisDirectory)

    end
  end
end

end
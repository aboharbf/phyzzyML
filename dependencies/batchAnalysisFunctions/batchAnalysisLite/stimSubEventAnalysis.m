function stimSubEventAnalysis(spikePathBank, batchAnalysisParams)
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
      subEventInd = contains(selTableVars, 'subSel') & contains(selTableVars, 'selInd');
      subEventVars = selTableVars(subEventInd);
      subEventName= extractBetween(subEventVars, 'subSel_', '_selInd');
      subEventNamePlot = strrep(subEventName, '_', ' ');
      
      subEventVarsCore = subEventVars(~contains(subEventName, 'sacc'));
      filtNames = subEventName(~contains(subEventName, 'sacc'));
      
      for event_i = 1:length(subEventVarsCore)
        
        % First, make a Venn Diagram comparing Sacc and Non-Sacc versions
        if event_i < 5
          eventAllInd = contains(subEventVars, filtNames{event_i});
          eventNamePlotAll = subEventNamePlot(eventAllInd)';
          eventNamesAll = subEventVars(eventAllInd);
          eventDataAll = selTableParadigmUnit{:, eventNamesAll};
          vennXExpanded(eventDataAll, strcat(strrep(filtNames{event_i}, '_', ' '), ' All'), {'All', 'Saccades', 'No Saccade'})      
        end
        
        if 0
        evTable = selTableParadigmUnit(selTableParadigmUnit.(subEventVarsCore{event_i}), :);
        unitScore = evTable.(strrep(subEventVarsCore{event_i}, 'selInd', 'pVal'));
        [~, sortInd] = sort(unitScore);
        unitTag = strcat(evTable.dateSubj, evTable.runNum, '_', evTable.channel, evTable.unitType);
        unitTagSorted = unitTag(sortInd);
        openUnitFigs(unitTagSorted(1:50), 'subEventPSTH*', subEventNamePlot{event_i}, batchAnalysisParams.analysisDirectory)
        end
        
      end
      
      % Pre-saccade overlap vs Unit
      subEventMat = any(selTableParadigmUnit{:, subEventVars(1:4)}, 2);
      preSaccMat = any(selTableParadigmUnit{:, 'subSel_pre_saccades_selInd'}, 2);
      vennXExpanded([subEventMat, preSaccMat], 'subEvent vs Pre-saccade', {'subEvent Selective', 'Pre-Saccade Selective'})

      
      
    end
  end

end
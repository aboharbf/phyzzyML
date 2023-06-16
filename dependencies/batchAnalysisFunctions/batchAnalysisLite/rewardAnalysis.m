function rewardAnalysis(spikePathBank, batchAnalysisParams)
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
      
      % Identify the 3 types of units we're interesting
      % First are classical reward responsive
      rewardCellUnits = selTableParadigmUnit(selTableParadigmUnit.subSel_rewardCombo_selInd, :);
      
      % Sort them on effect size
      rankingUnit = rewardCellUnits.subSel_rewardAbsent_pVal;
      [~, sortInd] = sort(rankingUnit);
      rankingUnit = rankingUnit(sortInd);
      rewardCellUnits = rewardCellUnits(sortInd, :);
      unitLabels = strcat(rewardCellUnits.dateSubj, rewardCellUnits.runNum, '_', rewardCellUnits.channel, rewardCellUnits.unitType);

      num2Check = 1:10;
      openUnitFigs(unitLabels(num2Check), 'subEventPSTH_*', '_reward',  batchAnalysisParams.analysisDirectory)
      % Unit - "20201209Mo005_Ch45U1"
      
      % Second, Reward Anticipation
      rewardAntCellUnits = selTableParadigmUnit(selTableParadigmUnit.subSel_rewardAnt_selInd, :);
      [~, sortInd] = sort(rewardAntCellUnits.subSel_rewardAnt_cohensD, 'descend');
      rewardAntCellUnits = rewardAntCellUnits(sortInd,:);
      unitLabels = strcat(rewardAntCellUnits.dateSubj, rewardAntCellUnits.runNum, '_', rewardAntCellUnits.channel, rewardAntCellUnits.unitType);
      
      num2Check = 1:10;
      openUnitFigs(unitLabels(num2Check), 'subEventPSTH_*', '_rewardAnt',  batchAnalysisParams.analysisDirectory)
      % Unit IDed - "20201204Mo003_Ch16U1"
      
      % Finally, Reward Anticipation coming from the stimulus all the way
      % to reward. trickier. 
      openUnitFigs(unitLabels(1:50), 'subEventPSTH_*', '_rewardAnt',  batchAnalysisParams.analysisDirectory)
      h =  findobj('type','figure');
      [figString, ChUString] = deal(cell(size(h)));
      for fig_i = 1:length(h)
        figString{fig_i} = h(fig_i).Children(1).Children.String;
        ChUString{fig_i} = h(fig_i).Name;
      end
      
      figString = strrep(figString, ',Run', '');
      ChUString = extractBefore(ChUString, ',');
      
      fig2CheckList = strcat(figString, '_', ChUString);
      openUnitFigs(fig2CheckList(11:22), 'subEventPSTH_*', '_rewardAnt',  batchAnalysisParams.analysisDirectory)
      openUnitFigs(fig2CheckList(11:22), 'imPSTH_*', [],  batchAnalysisParams.analysisDirectory)
      
      %20210624Sam001_Ch13U1 - Task transitions. 
      %20201204Mo003_Ch12U1 - Ramping to reward.
      %20201120Mo003_Ch9U1 20201120Mo003Ch2U1, 20201202MoRun003_Ch29U1 - Cognitive control? Ramping
      
      % Cognitive control units
      cogCon = {'20201120Mo003_Ch9U1','20201120Mo003_Ch2U1','20201202MoRun003_Ch29U1'};
      openUnitFigs(cogCon, 'subEventPSTH_*', '_rewardAnt',  batchAnalysisParams.analysisDirectory)
      openUnitFigs(cogCon, 'imPSTH_*', [],  batchAnalysisParams.analysisDirectory)
            
    end
  end
end

end
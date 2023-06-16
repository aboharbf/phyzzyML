function stimOffsetRewardAnalysis(spikePathBank, batchAnalysisParams)
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
      units2Find = 'rewardAnt';
      switch units2Find
        case 'stimOffset'
          selectedUnitTable = selTableParadigmUnit(selTableParadigmUnit.subSel_stimOffset_selInd, :);
          unitScore = selectedUnitTable.subSel_stimOffset_pVal;
          figTag2Find = 'subEventPSTH_*';
          filtTag = 'stimOffset';
        case 'rewardAnt'
          selectedUnitTable = selTableParadigmUnit(selTableParadigmUnit.subSel_rewardAnt_selInd, :);
          unitScore = selectedUnitTable.subSel_rewardAnt_pVal;
          figTag2Find = 'subEventPSTH_*';
          filtTag = 'rewardAnt';
        case 'reward'
          selectedUnitTable = selTableParadigmUnit(selTableParadigmUnit.subSel_reward_selInd, :);
          unitScore = selectedUnitTable.subSel_reward_pVal;
          figTag2Find = 'subEventPSTH_*';
          filtTag = 'reward';
      end
      
      unitLabels = strcat(selectedUnitTable.dateSubj, selectedUnitTable.runNum, '_', selectedUnitTable.channel, selectedUnitTable.unitType);
      [~, sortInd] = sort(unitScore);
      unitLabelsSorted = unitLabels(sortInd);
      
      %openUnitFigs(unitLabelsSorted(2), figTag2Find, filtTag, batchAnalysisParams.analysisDirectory)
      %openUnitFigs(unitLabelsSorted(2), 'imPSTH_*', [], batchAnalysisParams.analysisDirectory)

      % Make a Venn Diagram, of Stimulus Offset units (stimOffset w/ -
      % cohensD), Reward Anticipation (rewardAnt with +cohensD), Reward
      % (reward with +cohensD).
      stimOffsetSel = selTableParadigmUnit.subSel_stimOffset_selInd & (selTableParadigmUnit.subSel_stimOffset_cohensD < 0);
      rewardAntSel = selTableParadigmUnit.subSel_rewardAnt_selInd & (selTableParadigmUnit.subSel_rewardAnt_cohensD > 0);
      rewardSel = selTableParadigmUnit.subSel_rewardCombo_selInd;% & (selTableParadigmUnit.subSel_rewardCombo_cohensD > 0);
      %vennMat = [stimOffsetSel, rewardAntSel, rewardSel];
      %vennMat(900:end, :) = true;
      
      % Stim Offset + and Reward Ant +
%       A = selTableParadigmUnit.subSel_stimOffset_selInd & (selTableParadigmUnit.subSel_stimOffset_cohensD > 0);
%       B = selTableParadigmUnit.subSel_rewardAnt_selInd & (selTableParadigmUnit.subSel_rewardAnt_cohensD > 0);
%       vennXExpanded([A B], 'Stimulus Offset and Reward Activity', {'Stimulus Offset', 'Reward Anticipation'});
%       CheckInd = find([A & B]);
%       openUnitFigs(selTableParadigmUnit(CheckInd(1:10), :), figTag2Find, 'rewardAnt', batchAnalysisParams.analysisDirectory)
%       openUnitFigs(selTableParadigmUnit(CheckInd(1:10), :), figTag2Find, 'stimOffset', batchAnalysisParams.analysisDirectory)
      %openUnitFigs(selTableParadigmUnit(stimOffsetSel, :), figTag2Find, 'stimOffset', batchAnalysisParams.analysisDirectory)
      rewardAntInd = [4 6 16 17 18];
      
      index2Show = find(rewardAntSel);
      openUnitFigs(selTableParadigmUnit(index2Show(rewardAntInd), :), 'imPSTH_*', [], batchAnalysisParams.analysisDirectory)
      openUnitFigs(selTableParadigmUnit(index2Show(rewardAntInd), :), figTag2Find, 'rewardAnt', batchAnalysisParams.analysisDirectory)

      vennMat = [stimOffsetSel, rewardAntSel, rewardSel];
      vennXExpanded(vennMat, 'Stimulus Offset and Reward Activity', {'Stimulus Offset', 'Reward Anticipation', 'Reward'});
      
      vennMat = [stimOffsetSel, rewardAntSel];
      vennXExpanded(vennMat, 'Stimulus Offset and Reward Activity', {'Stimulus Offset', 'Reward Anticipation'});
      
      vennMat = [stimOffsetSel, rewardSel];
      vennXExpanded(vennMat, 'Stimulus Offset and Reward Activity', {'Stimulus Offset', 'Reward'});
      
      stimAny = any([selTableParadigmUnit.baseV_stimEarly_selInd, selTableParadigmUnit.baseV_stimLate_selInd], 2);
      vennMat = [rewardSel stimAny selTableParadigmUnit.baseV_Fix_selInd];
      vennXExpanded(vennMat, 'Reward Activity with Epoch selectivity', {'Reward Event', 'Stimulus Epoch', 'Fixation Epoch'});
      
      % Reward + Reward Ant?
      curiousPop = rewardAntSel & rewardSel;
      selectedUnitTable = selTableParadigmUnit(curiousPop,:);
      unitLabels = strcat(selectedUnitTable.dateSubj, selectedUnitTable.runNum, '_', selectedUnitTable.channel, selectedUnitTable.unitType);
      openUnitFigs("20201204Mo003_Ch12U1", figTag2Find, 'reward', batchAnalysisParams.analysisDirectory)
      % openUnitFigs(unitLabels, figTag2Find, filtTag, batchAnalysisParams.analysisDirectory)
      % Reward is Rewarded trials vs Non-rewarded, Reward Ant is Pre vs
      % Post Reward period, so a unit can have both, but its more likely
      % that the 'Reward Ant' is the stronger signal.
      % Anticipation + Reward response - "20201204Mo003_Ch12U1"
      
      % All Event vs All Epoch
      eventSelEvents = selTableParadigmUnit.Properties.VariableNames(contains(selTableParadigmUnit.Properties.VariableNames, 'subSel_') & contains(selTableParadigmUnit.Properties.VariableNames, '_selInd'));
      eventSelEvents = eventSelEvents(~contains(eventSelEvents, 'sacc'));
      eventSelEvents = eventSelEvents(~contains(eventSelEvents, 'rewardA'));
      eventSelEvents = eventSelEvents(~contains(eventSelEvents, 'reward_'));
      eventSelEvents = eventSelEvents(~contains(eventSelEvents, 'blinks'));
      
      epochSelEvents = selTableParadigmUnit.Properties.VariableNames(contains(selTableParadigmUnit.Properties.VariableNames, 'baseV_') & contains(selTableParadigmUnit.Properties.VariableNames, '_selInd'));

      eventSelInd = any(selTableParadigmUnit{:, eventSelEvents}, 2);
      epochSelInd = any(selTableParadigmUnit{:, epochSelEvents}, 2);
      
      vennMat = [eventSelInd, epochSelInd];
      vennXExpanded(vennMat, 'Event Selectivity and Epoch Modulation', {'Event', 'Epoch'});
      
      
    end
  end
end

end
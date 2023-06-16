function StimOnSetSaccadeAnalysis(spikePathBank, batchAnalysisParams)
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
    for unitType_i = 1:length(unitList)
      
      % Set output directory
      selParams.outputDir = fullfile(mainOutDir, paradigmList{par_i}, monkeyList{m_i}, unitList{unitType_i});
      
      unitIndex = contains(selTableParadigm.unitType, unitTags{unitType_i});
      selTableParadigmUnit = selTableParadigm(unitIndex, :);
      selParams.unitTag = unitList{unitType_i};
    
      %selTableSaccDir = selTableParadigmUnit(selTableParadigmUnit.subSel_pre_saccades_selInd, :);
      selTableSaccDir = selTableParadigmUnit(selTableParadigmUnit.subSel_stimOnset_selInd, :);
      
      runsWithUnits = strcat('S', unique(strcat(selTableSaccDir.dateSubj, selTableSaccDir.runNum)));
      
      fixDotspikePath = spikePathBank(runsWithUnits,:);
      
      % Load things
      [psthPerStimPerRun, eyeBehDataPerRun, psthParamsPerRun] = spikePathLoad(fixDotspikePath, {'psthByImage', 'eyeDataStruct', 'psthParams'}, batchAnalysisParams.spikePathLoadParams);
      psthPrePerRun = cellfun(@(x) x.psthPre, psthParamsPerRun);
      
      % Create storage structures for this
      meanStimOnsetPeak = [];
      meanFirstSaccadeTime = [];
      for run_i = 1:length(psthPerStimPerRun)
        runSelTable = fixDotspikePath.selTable{run_i};
        runSelTableStimOnset = runSelTable(runSelTable.subSel_stimOnset_selInd,:);
        unitIndex = contains(runSelTableStimOnset.unitType, unitTags{unitType_i});
        runSelTableStimOnsetUnit = runSelTableStimOnset(unitIndex,:);
        
        unitTag = strcat(runSelTableStimOnsetUnit.dateSubj, runSelTableStimOnsetUnit.runNum, '_', runSelTableStimOnsetUnit.channel, runSelTableStimOnsetUnit.unitType);
        
        % Create ref vector for channels
        chanNumVecAll = sort(str2double(extractAfter(unique(runSelTable.channel), 'Ch')));
        chanNumVecUnits = str2double(extractAfter(runSelTableStimOnsetUnit.channel, 'Ch'));
        unitNum = str2double(extractAfter(runSelTableStimOnsetUnit.unitType, 'U'))+1;
        psthPre = psthPrePerRun(run_i);
        
        for unit_i = 1:size(runSelTableStimOnsetUnit,1)
          chInd = chanNumVecAll == chanNumVecUnits(unit_i);
          psthPerStimUnit = psthPerStimPerRun{run_i}{chInd}{unitNum(unit_i)};
          eyeBehPerStim = eyeBehDataPerRun{run_i};
          stimEarlyActivity = psthPerStimUnit(:,psthPre:psthPre+500);
          
          peakTimePerStim = nan(size(psthPerStimUnit,1),1);
          for stim_i = 1:size(psthPerStimUnit,1)
            [~, peakTimePerStim] = max(stimEarlyActivity(stim_i,:));
            
            saccadeTimes = eyeBehPerStim.saccadeByStim{stim_i}(:,psthPre:psthPre+500) == 3;     
            saccStim = nan(size(saccadeTimes,1),1);
            for ii = 1:size(saccadeTimes,1)
              saccTime = find(saccadeTimes(ii,:), 1);
              if ~isempty(saccTime)
                saccStim(ii) = find(saccadeTimes(ii,:), 1);
              end
            end
            
            meanFirstSaccadeTime = [meanFirstSaccadeTime; mean(saccStim, 'omitnan')];
            meanStimOnsetPeak = [meanStimOnsetPeak; peakTimePerStim];

          end
          
          
        end
        
      end
      
      figure(); 
      hist1 = histogram(meanFirstSaccadeTime); 
      hold on; 
      hist2 = histogram(meanStimOnsetPeak);
      
      hist2.BinEdges = hist1.BinEdges;
      
      figure();
      hist3 = histogram(meanStimOnsetPeak - meanFirstSaccadeTime);

    end
  end
end

end
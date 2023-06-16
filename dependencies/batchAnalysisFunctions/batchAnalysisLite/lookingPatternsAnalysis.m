function lookingPatternsAnalysis(spikePathBank, batchAnalysisParams)
% A function which creates bar plots represeting how often specific objects
% in particular scenes are attended to.

paradigmList = unique(spikePathBank.paradigmName);

allObjPres = {'Face1', 'Face2', 'Body1', 'Body2', 'HandR1', 'HandR2', 'bkg'};
% monkeyList = {'Sam', 'Mo', 'Combo'};
monkeyList = {'Combo'};

statsPerMonkeyPar = cell(length(monkeyList), length(paradigmList));
parStim = cell(1, length(paradigmList));
for m_i = 1:length(monkeyList)
  for par_i = 1:length(paradigmList)
    
    pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
    if ~strcmp(monkeyList{m_i}, 'Combo')
      mInd = contains(spikePathBank.Row, monkeyList{m_i});
    else
      mInd = true(size(pInd));
    end
    
    spikePathBankParadigm = spikePathBank(pInd & mInd, :);
    
    % Compile a list of all the stimuli present
    uniqueStimList = unique([spikePathBankParadigm.stimuli{:}]);
    [parStim{par_i}, newSortInd] = createPlotNamesForStimuli(uniqueStimList);
    uniqueStimList = uniqueStimList(newSortInd);
    
    % Load the eyeDataStruct for all the runs.
    [eyeDataStructByRun] = spikePathLoad(spikePathBankParadigm, {'eyeDataStruct'}, batchAnalysisParams.spikePathLoadParams);
    
    objAttStats = cell(size(uniqueStimList));
    [objAttStats{:}] = deal(zeros(0, length(allObjPres)));
    for run_i = 1:length(eyeDataStructByRun)
      
      attendObjVecRun = eyeDataStructByRun{run_i}.attendedObjData.attendedObjVect;
      
      for stim_i = 1:length(attendObjVecRun)
        stimName = spikePathBankParadigm.stimuli{run_i}{stim_i};
        
        % Create a grid to store the desired stats
        trialCount = size(attendObjVecRun{stim_i}, 1);
        statGrid = nan(trialCount, length(allObjPres));
        
        for trial_i = 1:size(attendObjVecRun{stim_i},1)
          for obj_i = 1:length(allObjPres)
            statGrid(trial_i, obj_i) = sum(cellfun(@(x) strcmp(x, allObjPres(obj_i)), attendObjVecRun{stim_i}(trial_i,:)), 'all');
          end
        end
        
        stimInd = strcmp(uniqueStimList,stimName);
        if ~isempty(stimInd)
          objAttStats{stimInd} = cat(1, objAttStats{stimInd}, statGrid);
        end
        
      end
    end
    
    statsPerMonkeyPar{m_i, par_i} = objAttStats;
    
    % for each run, and each stimulus, convert the grid of objects/trials into a fraction per object present.
  end
end

% In the case of headTurnCon, merge the stimuli appropriately.
htcInd = strcmp(paradigmList, 'headTurnCon');
if any(htcInd)
  for m_i = 1:length(monkeyList)
    tmp = statsPerMonkeyPar{m_i, htcInd};
    [uniqueStimList, ~, uniB] = unique(parStim{htcInd});
    parStim{htcInd} = uniqueStimList;
    
    % Concatonate
    stimuliPlotListTmp = cell(size(uniqueStimList));
    for stim_i = 1:length(uniqueStimList)
      stim2Merge = uniB == stim_i;
      stimuliPlotListTmp{stim_i} = cat(1, tmp{stim2Merge});
    end
    
    statsPerMonkeyPar{m_i, htcInd} = stimuliPlotListTmp;
  end
end

% Convert those frame counts into means and SDs.
meanStatsPerStim = initNestedCellArray(statsPerMonkeyPar, 'zeros', [2 length(allObjPres)], 100);
for m_i = 1:length(monkeyList)
  for par_i = 1:length(paradigmList)
    for stim_i = 1:length(statsPerMonkeyPar{m_i, par_i})
      stimData = statsPerMonkeyPar{m_i, par_i}{stim_i};
      
      meanStatsPerStim{m_i, par_i}{stim_i}(1,:) = mean(stimData, 1);
      meanStatsPerStim{m_i, par_i}{stim_i}(2,:) = std(stimData, 0, 1)/size(stimData,1);
      
    end
  end
end


% Plot time
collapseObjs = true;      % Collapse across 1 and 2.
for par_i = 1:length(paradigmList)
  
  if contains(paradigmList{par_i}, 'Natural')
    keepInd = contains(parStim{par_i}, 'monkey');
  else
    keepInd = ~contains(parStim{par_i}, {'land', 'obj'});
  end
  
  for m_i = 1:length(monkeyList)
    
    data2Plot = meanStatsPerStim{m_i, par_i};
    data2Plot = cat(3, data2Plot{:});
    
    % Means
    meanPerObjStim = squeeze(data2Plot(1,:,:));
    meanPerObjStim = meanPerObjStim(1:end-1,:)./84; % Get rid of bkg.
        
    % Filter
    stimuliPlotList = strrep(parStim{par_i}(keepInd), '_', ' ');
    meanPerObjStim = meanPerObjStim(:, keepInd);
    
    if collapseObjs
      params.labels = {'Face', 'Body', 'Hand'};
      params.colors = [255 0 0; 0 0 255; 0 255 0]/255;
      clear tmp2
      tmp2(1,:) = sum(meanPerObjStim([1 2], :),1);
      tmp2(2,:) = sum(meanPerObjStim([3 4], :),1);
      tmp2(3,:) = sum(meanPerObjStim([5 6], :),1);
      meanPerObjStim = tmp2;
    else
      params.labels = {'Face1', 'Face2', 'Body1', 'Body2', 'HandR1', 'HandR2', 'bkg'};
      params.colors = [255 0 0; 130 0 0; 0 0 255; 0 0 130; 0 255 0; 0 130 0; 130 130 130]/255;
    end
    params.textSwitch = false;
    
    % Switch to percent as out of 100
    meanPerObjStim = meanPerObjStim * 100;
        
    %   subplot(2,1,1)
    segmentCount = 4;
    perRowCount = ceil(length(stimuliPlotList)/segmentCount);
    segmentIndStarts = [1:perRowCount:length(stimuliPlotList)];
    segmentEndInds = [segmentIndStarts(2:end)-1 length(stimuliPlotList)];
    
    for ii = 1:length(segmentIndStarts)
      figH = createBarPlotWithChanceLine(stimuliPlotList(segmentIndStarts(ii):segmentEndInds(ii)), meanPerObjStim(:, segmentIndStarts(ii):segmentEndInds(ii)), 0, 0, 'Looking time per Object', params.labels, params);
      delete(figH.Children(1))
      figH.Children(end).FontSize = 8;
      figH.Position = [0.1184    0.2035    0.4398    0.1986];
      if ii ~= 1
        figH.Children(end).Title.String = '';        
      end
      figH.Children.FontSize
      ylabel('% Looking time');
      ylim([0 80])
    end
    
    
  end
end

end
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
collapseObjs = true;      % Collapse across 1 and 2 (i.e. 'Face1' and 'Face2' -> 'Face')
mergeHands = true;        % For consistency, merge all 'Hand' counts for chasing/fighting into body
handBodyMergeIdx = contains(uniqueStimList, {'Chasing', 'Fighting'});

meanStatsPerStim = initNestedCellArray(statsPerMonkeyPar, 'zeros', [2 length(allObjPres)], 100);
for m_i = 1:length(monkeyList)
  for par_i = 1:length(paradigmList)
    for stim_i = 1:length(statsPerMonkeyPar{m_i, par_i})
      stimData = statsPerMonkeyPar{m_i, par_i}{stim_i};
      
      % Merge hands into bodies for chasing and fighting vids.
      if mergeHands
        if handBodyMergeIdx(stim_i)
          newB1 = stimData(:, 3) + stimData(:, 5);
          newB2 = stimData(:, 4) + stimData(:, 6);
          stimData = [stimData(:,1:2), newB1, newB2, nan(size(newB1)), nan(size(newB1)), stimData(:,7)];
        end
      end
      
      if collapseObjs
        stimData = [sum(stimData(:,1:2), 2), sum(stimData(:,3:4), 2), sum(stimData(:,5:6), 2), stimData(:,7)];
      end
      
      meanStatsPerStim{m_i, par_i}{stim_i}(1,:) = mean(stimData, 1);
      meanStatsPerStim{m_i, par_i}{stim_i}(2,:) = std(stimData, 0, 1)/size(stimData,1);
      
    end
  end
end


% Plot time
if collapseObjs
  params.labels = {'Face', 'Body', 'Hand'};
  params.colors = [255 0 0; 0 0 255; 0 255 0]/255;
else
  params.labels = {'Face1', 'Face2', 'Body1', 'Body2', 'HandR1', 'HandR2', 'bkg'};
  params.colors = [255 0 0; 130 0 0; 0 0 255; 0 0 130; 0 255 0; 0 130 0; 130 130 130]/255;
end
    
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
    meanPerObjStim = meanPerObjStim(1:end-1,:)./83; % Get rid of bkg, make percent of frames
        
    % std
    stdPerObjStim = squeeze(data2Plot(2,:,:));
    stdPerObjStim = stdPerObjStim(1:end-1,:)./83; % Get rid of bkg, make percent of frames
    
    % Filter
    stimuliPlotList = strrep(parStim{par_i}(keepInd), '_', ' ');
    meanPerObjStim = meanPerObjStim(:, keepInd);
    stdPerObjStim = stdPerObjStim(:, keepInd);
    
    params.textSwitch = false;
    
    % Switch to percent as out of 100
    meanPerObjStim = meanPerObjStim * 100;
    stdPerObjStim = stdPerObjStim * 100;
        
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
      ylabel('% Looking time');
      ylim([0 80])
    end
    
    
  end
end

end
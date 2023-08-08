function rewardPSTHallStack(spikePathBank, params)

outputDir = fullfile(params.subEventPSTHParams.outputDir, 'rewardPSTH');

normalizeSubEvent = 1;
plotPlots = 1;

% Extract data from meanPSTHStruct
paradigmList = unique(spikePathBank.paradigmName);

for par_i = 1:length(paradigmList)
  
  % Pull paradigm specific information
  paradigmIndex = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  spikePathBankParadigm = spikePathBank(paradigmIndex,:);
  selTablePerRun = spikePathBankParadigm.selTable;
  
  % Create a large stacked combination of psthes for the subEvents,
  % specifically the index to grab in the collected psthes per unit.
  subEventSigStructPerRun = spikePathLoad(spikePathBankParadigm, {'subEventSigStruct'}, params.spikePathLoadParams);
  if normalizeSubEvent
    [psthByImageFixAlignedPerRunUnProc, stimTimingByRun, spikeAlignParamsByRun] = spikePathLoad(spikePathBankParadigm, {'psthByCategoryFixAlign', 'stimTiming', 'spikeAlignParams'}, params.spikePathLoadParams);
    subEventSigStructPerRun = normalizeTraces(subEventSigStructPerRun, psthByImageFixAlignedPerRunUnProc, stimTimingByRun, spikeAlignParamsByRun);
    normTag = 'Normalized';
  else
    normTag = '';
  end
  
  % Stack
  subEventSigStructPerRun = [subEventSigStructPerRun{:}];
  
  % Find units per run
  unitsPerRun = cellfun(@(x) size(x,1), selTablePerRun);
  selTable = vertcat(selTablePerRun{:});
  traceStorage = cell(2, 2, 2);

  for reward_Type_i = 1
    
    if reward_Type_i == 1 % Grab the rewardAbsentTrace
      rewardAbsentIndex = cellfun(@(x) find(strcmp(x, 'rewardAbsent')), {subEventSigStructPerRun.events}, 'UniformOutput', false);
      rewardAbsentIndex(cellfun('isempty', rewardAbsentIndex)) = deal({NaN}); % Fill empty cells w/ NaN.
      rewardAbsentIndex = [rewardAbsentIndex{:}]';
      rewardAbsentIndex = arrayfun(@(x) repmat(rewardAbsentIndex(x), [unitsPerRun(x), 1]), 1:length(unitsPerRun), 'UniformOutput', false)';
      rewardAbsentIndex = vertcat(rewardAbsentIndex{:});
      index2Pull = rewardAbsentIndex;
      
      rewardTag = 'Reward vs Unrewarded trials';
      legendTitles = {'Rewarded Trials', 'Unrewarded Trials'};
      
    else
      rewardIndex = cellfun(@(x) find(strcmp(x, 'reward')), {subEventSigStructPerRun.events}, 'UniformOutput', false);
      rewardIndex = [rewardIndex{:}]';
      rewardIndex = arrayfun(@(x) repmat(rewardIndex(x), [unitsPerRun(x), 1]), 1:length(unitsPerRun), 'UniformOutput', false)';
      rewardIndex = vertcat(rewardIndex{:});
      index2Pull = rewardIndex;
      
      rewardTag = 'Reward vs Scramble';
      legendTitles = {'Rewarded Trials', 'Scramble'};

    end
    
    % Expand to logical
    logicalIndex2Pull = false(length(index2Pull), max(index2Pull));
    for jj = 1:length(index2Pull)
      if ~isnan(index2Pull(jj))
        logicalIndex2Pull(jj,index2Pull(jj)) = true;
      end
    end
    
    
    % Collect those traces
    subEventPSTH = {subEventSigStructPerRun.subEventPSTH}';
    subEventPSTH = vertcat(subEventPSTH{:});
    subEventDataPSTH = vertcat(subEventPSTH{:,1});
    subEventDataErrPSTH = vertcat(subEventPSTH{:,2});
    
    subEventPSTHNull = {subEventSigStructPerRun.subEventNullPSTH};
    subEventPSTHNull = vertcat(subEventPSTHNull{:});
    subEventNullPSTH = vertcat(subEventPSTHNull{:,1});
    subEventNullErrPSTH = vertcat(subEventPSTHNull{:,2});
    
    data2Process = {subEventDataPSTH, subEventDataErrPSTH, subEventNullPSTH, subEventNullErrPSTH};
    for ii = 1:length(data2Process)
      tmp = data2Process{ii};
      data2Process{ii} = arrayfun(@(x) tmp{x}(logicalIndex2Pull(x,:),:), 1:length(tmp), 'UniformOutput', false)';
      emptyInd = ~any(logicalIndex2Pull,2);
    end
    
    % Params
    psthParams.psthPre = -subEventSigStructPerRun(1).psthWindow(1);
    psthParams.psthImDur = subEventSigStructPerRun(1).psthWindow(2);
    psthParams.psthPost = 100;
    
    % For each run, identify which units are to be plotted
    for unitType_i = 1:2
      for sel_i = 1:2
        
        if sel_i == 1
          selIndex = selTable.subSel_rewardCombo_selInd;
          selTag = 'Selective';
        else
          selIndex = true(size(selTable.subSel_rewardCombo_selInd));
          selTag = 'All';
        end
        
        if unitType_i == 1
          unitIndex = contains(selTable.unitType, 'MUA') & selIndex & ~emptyInd;
          unitTag = 'MUA';
        else
          unitIndex = contains(selTable.unitType, digitsPattern) & selIndex & ~emptyInd;
          unitTag = 'U';
        end
        
        % Grab only relevant activity.
        
        data2ProcessUnit = cell(size(data2Process));
        for ii = 1:length(data2Process)
          data2ProcessUnit{ii} = vertcat(data2Process{ii}{unitIndex});
        end
        
        figTitle = sprintf('%s %s %s (%d) %s', rewardTag, selTag, unitTag, sum(unitIndex), normTag);
        dataTraces = [nanmean(data2ProcessUnit{1}, 1); nanmean(data2ProcessUnit{3}, 1)];
        %       errorTraces = [mean(data2Process{2}, 1); mean(data2Process{4}, 1)];
        errorTraces = [nanstd(data2ProcessUnit{1}, 1)/sqrt(sum(unitIndex)); nanstd(data2ProcessUnit{3}, 1)/sqrt(sum(unitIndex))];
        
        traceStorage{reward_Type_i, unitType_i, sel_i} = {dataTraces, errorTraces};
        
        if plotPlots
          figH = figure('Name', figTitle, 'Units', 'normalized', 'position', [0.35 0.4 0.35 0.5]);
          [axesH, ~, ~, legendH] = plotPSTH(dataTraces, errorTraces, [], psthParams, 'line', figTitle, legendTitles);
          legendH.Location = 'northeast';
          xlimSize = xlim();
          xlim([xlimSize(1), xlimSize(2)-300]);
          xlabel('Time from Reward Delivery');
          axesH.FontSize = 16;
          if normalizeSubEvent
            ylabel('Normalized activity');
          end
          
          saveFigure(outputDir, sprintf('1. %s - %s', paradigmList{par_i}, figTitle), [], params.figStruct, []);
        end
        
      end
      
    end
  end
  
  % Combine traces correctly
%   unitType_i = 1; unitTag = 'MUA';
%   reward_Type_i = 1;                  % rewardAbsent traces
%   
%   rewardSelTrace = traceStorage{reward_Type_i, unitType_i, 1}{1}(1,:);
%   rewardSelTraceErr = traceStorage{reward_Type_i, unitType_i, 1}{2}(1,:);
%   
%   rewardSelTraceNull = traceStorage{reward_Type_i, unitType_i, 1}{1}(2,:);
%   rewardSelTraceNullErr = traceStorage{reward_Type_i, unitType_i, 1}{2}(2,:);
%   
%   rewardNonSelTrace = traceStorage{reward_Type_i, unitType_i, 2}{1}(1,:);
%   rewardNonSelTraceErr = traceStorage{reward_Type_i, unitType_i, 2}{2}(1,:);
%   
%   dataTraces = [rewardSelTrace; rewardSelTraceNull; rewardNonSelTrace];
%   errorTraces = [rewardSelTraceErr; rewardSelTraceNullErr; rewardNonSelTraceErr];
%   
%   figTitle = sprintf('Reward vs Non-Reward Activity - %s %s', unitTag, normTag);
%   figH = figure('Name', figTitle, 'Units', 'normalized', 'position', [0.35 0.4 0.35 0.5]);
%   [~, ~, ~, legendH] = plotPSTH(dataTraces, errorTraces, [], psthParams, 'line', figTitle, legendTitles);
%   legendH.Location = 'northeast';
%   xlimSize = xlim();
%   xlim([xlimSize(1), xlimSize(2)-300]);
%   xlabel('Time from Reward Delivery');
%   figH.Children(2).FontSize = 16;
%   if normalizeSubEvent
%     ylabel('Normalized activity');
%   end
%   
%   saveFigure(outputDir, ['1. ' figTitle], [], params.figStruct, []);
  
end
function rewardPSTH(spikePathBank, params)

outputDir = fullfile(params.subEventPSTHParams.outputDir, 'rewardPSTH');

if ~exist(outputDir, 'dir')
  mkdir(outputDir);
end

% Extract data from meanPSTHStruct
paradigmList = unique(spikePathBank.paradigmName);

for par_i = 1:length(paradigmList)
  
  % Pull paradigm specific information
  paradigmIndex = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  spikePathBankParadigm = spikePathBank(paradigmIndex,:);
    
  % Create a large stacked combination of psthes for the subEvents,
  % specifically the index to grab in the collected psthes per unit.
  subEventSigStructPerRun = spikePathLoad(spikePathBankParadigm, {'subEventSigStruct'}, params.spikePathLoadParams);
  subEventSigStructPerRun = [subEventSigStructPerRun{:}];
  
  % Identify indices for rewards
  rewardAbsentIndex = cellfun(@(x) find(strcmp(x, 'rewardAbsent')), {subEventSigStructPerRun.events}, 'UniformOutput', false);
  rewardAbsentIndex(cellfun('isempty', rewardAbsentIndex)) = deal({NaN}); % Fill empty cells w/ NaN.
  rewardIndex = cellfun(@(x) find(strcmp(x, 'reward')), {subEventSigStructPerRun.events}, 'UniformOutput', false);
  
  % Collect those traces
  subEventPSTH = {subEventSigStructPerRun.subEventPSTH};
  subEventNullPSTH = {subEventSigStructPerRun.subEventNullPSTH};
  selTablePerRun = spikePathBankParadigm.selTable;
  
  % Params
  psthParams.psthPre = -subEventSigStructPerRun(1).psthWindow(1);
  psthParams.psthImDur = subEventSigStructPerRun(1).psthWindow(2);
  psthParams.psthPost = 100;
  
  % Cycle through unitTypes
    
  % For each run, identify which units are to be plotted
  for run_i = 1:length(subEventSigStructPerRun)
    % Grab relevant data
    subEventParams = subEventSigStructPerRun(run_i);
    subEventRunPSTH = vertcat(subEventPSTH{run_i}{:,1});
    subEventRunPSTHErr = vertcat(subEventPSTH{run_i}{:,2});
    subEventNullRunPSTH = vertcat(subEventNullPSTH{run_i}{:,1});
    subEventNullRunPSTHErr = vertcat(subEventNullPSTH{run_i}{:,2});
    
    selTableRun = selTablePerRun{run_i};
    assert(length(subEventRunPSTH) == size(selTableRun,1), 'mismatch in run data')
    
    % check which reward it has
    rewardAbsentParadigm = ~all(selTableRun.subSel_rewardAbsent_pVal == 1);
    
    for unitType_i = 2%1:2
      
      if unitType_i
        selTableRunUnit = selTableRun(contains(selTableRun.unitType, 'MUA'), :);
      else
        selTableRunUnit = selTableRun(~contains(selTableRun.unitType, 'MUA'), :);
      end
      
      % Extract the reward index to be used
      if rewardAbsentParadigm
        rewardSelIndex = selTableRunUnit.subSel_rewardAbsent_selInd;
        rewardIndex = strcmp(subEventParams.events, 'rewardAbsent');
      else
        rewardSelIndex = selTableRunUnit.subSel_reward_selInd;
        rewardIndex = strcmp(subEventParams.events, 'reward');
      end
      
      % Use the selectivity index to plot selective units
      selSites = rewardSelIndex;
      siteInfo = strcat(selTableRunUnit.dateSubj, selTableRunUnit.runNum, '-', selTableRunUnit.channel, selTableRunUnit.unitType);
      siteInfo = siteInfo(selSites);
      
      extractTraces = @(x) x(rewardIndex,:);
      rewardTracePerUnit = cellfun(extractTraces, subEventRunPSTH, 'UniformOutput', false);
      rewardTracePerUnit = vertcat(rewardTracePerUnit{rewardSelIndex});
      
      rewardTraceErrPerUnit = cellfun(extractTraces, subEventRunPSTHErr, 'UniformOutput', false);
      rewardTraceErrPerUnit = vertcat(rewardTraceErrPerUnit{rewardSelIndex});
      
      rewardTraceNullPerUnit = cellfun(extractTraces, subEventNullRunPSTH, 'UniformOutput', false);
      rewardTraceNullPerUnit = vertcat(rewardTraceNullPerUnit{rewardSelIndex});
      
      rewardTraceNullErrPerUnit = cellfun(extractTraces, subEventNullRunPSTHErr, 'UniformOutput', false);
      rewardTraceNullErrPerUnit = vertcat(rewardTraceNullErrPerUnit{rewardSelIndex});
      
      for ii = 1:size(rewardTracePerUnit, 1)
        figTitle = sprintf('Reward selective unit - %s', siteInfo(ii));
        dataTraces = [rewardTracePerUnit(ii,:); rewardTraceNullPerUnit(ii,:)];
        errorTraces = [rewardTraceErrPerUnit(ii,:); rewardTraceNullErrPerUnit(ii,:)];
        
        figH = figure('Name', figTitle, 'Units', 'normalized', 'position', [0.35 0.5 0.5 0.4]);
        plotPSTH(dataTraces, errorTraces, [], psthParams, 'line', figTitle, {'Event', 'Null'});
        
      end
      
    end
  end
  
end

end
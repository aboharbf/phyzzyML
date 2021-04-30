function spikePathBank = processAppendSelTable(spikePathBank, params)
% A function which loads the selTable variable, processes it for use, and
% appends it back to the spikePathBank.

% Collect unit selectivity
[selTablePerRun, anovaTablePerRun] = spikePathLoad(spikePathBank, {'selTable', 'anovaTable'}, params.spikePathLoadParams);
alpha = params.selParam.alpha;
strThres = params.selParam.stretchThreshold;
objectStretches = params.selParam.objectStretches;                             

% Keep the tables seperate, since they have to be returned to
% spikePathBank.
processedSelTableArray = cell(size(selTablePerRun));

% Variables to exclude
varExcludeanova = {'slidingWin_socIntTest_Eye'};
varExcludesel = {'_cohensD'};

for run_i = 1:length(selTablePerRun)
  
  selTableRun = selTablePerRun{run_i};
  
  % Remove undesired variables
  var2Remove = contains(selTableRun.Properties.VariableNames, varExcludesel);
  selTableRun(:, var2Remove) = [];
  
  % Check for pVal rows.
  pValIndex = contains(selTableRun.Properties.VariableNames, '_pVal');
  pValVars = selTableRun.Properties.VariableNames(pValIndex)';
  selIndexVars = strcat(extractBefore(pValVars, '_pVal'), '_selInd');
  
  % Collect pVals and convert them to selectivity indicies, based on
  % previous defined alpha.
  pValValues = selTableRun{:, pValVars};
  selIndexValues = pValValues < alpha;
  
  % Add them back to the table.
  selTableRun{:, selIndexVars} = selIndexValues;
   
  % Process the ANOVA table
  anovaTableRun = anovaTablePerRun{run_i};
  var2Remove = contains(anovaTableRun.Properties.VariableNames, varExcludeanova);
  anovaTableRun(:, var2Remove) = []; 
  
  % Check for pVal rows in the ANOVA table
  pValIndex = contains(anovaTableRun.Properties.VariableNames, '_pVal');
  pValVars = anovaTableRun.Properties.VariableNames(pValIndex)';
  selIndexVars = strcat(extractBefore(pValVars, '_pVal'), '_selInd');
  selIndexValues = zeros(size(anovaTableRun,1), length(pValVars));
  
  for ii = 1:length(pValVars)
    pValValues = anovaTableRun{:, pValVars{ii}};
    dataStack = vertcat(pValValues{:});
    pValInd2Keep = dataStack < alpha;
    
    % Find stretches of significant runs, threshold, and store.
    selIndexValues(:, ii) = findStretchRuns(pValInd2Keep) >= strThres;
    
  end
  
  % Return to results to selTable instead.
  selTableRun{:, selIndexVars} = selIndexValues;
    
  % For the preferred objects, things are more complicated. collect the
  % data, stack it, and expand to a unit*bin*uniqueObj logical array.
  prefStimInd = find(contains(anovaTableRun.Properties.VariableNames, '_prefStim'));
  prefStimName = anovaTableRun.Properties.VariableNames(prefStimInd);
  selIndexVars = strcat(extractBefore(prefStimName, '_prefStim'), '_selInd');
  preStimObjVec = strcat(extractBefore(prefStimName, '_prefStim'), '_prefStimObj');

  for p_val_ind = 1:length(prefStimInd)
    prefStimArray = anovaTableRun{:, prefStimInd(p_val_ind)};
    dataStack = vertcat(prefStimArray{:});
    
    uniqueStimLabels = unique(dataStack(:));
    uniqueStimLabels = uniqueStimLabels(~contains(uniqueStimLabels, 'None'));

    dataStackFinalTmp = zeros([size(dataStack), length(uniqueStimLabels)]);
    for ii = 1:length(uniqueStimLabels)
      dataStackFinalTmp(:,:,ii) = strcmp(dataStack, uniqueStimLabels{ii}) * ii;
    end
    
    % This section may be expanded in the future to keep stock of which
    % objects selectivity was shown for. Not necessary here.
    sigStretchPerObj = squeeze(findStretchRuns(dataStackFinalTmp) >= strThres);
    
    % store the objects
    objVec = cell(size(sigStretchPerObj,1),1);
    uniqueStimLabelCombos = uniqueStimLabels;
    uniqueStimLabelCombos(2:end) = strcat(',', uniqueStimLabelCombos(2:end));
    
    for obj_i = 1:length(uniqueStimLabelCombos)
      objVec(sigStretchPerObj(:, obj_i)) = strcat(objVec(sigStretchPerObj(:, obj_i)), uniqueStimLabelCombos{obj_i});
    end

    % If you only care about stretches for the same objects, filter here.
    if objectStretches
      sigObjStretchInd = any(sigStretchPerObj,2);
    else
      sigObjStretchInd = findStretchRuns(any(dataStackFinalTmp,3)) >= strThres;
    end
    
    % Add the index for units with significant stretches + the objects.
    selTableRun.(selIndexVars{p_val_ind}) = sigObjStretchInd;
    selTableRun.(preStimObjVec{p_val_ind}) = objVec;
  end
  
  % Replace the channel names for 20201123Mo
  if strcmp(selTableRun{1,1}, '20201123Mo')
    selTableRun = replaceChanNum_1123(selTableRun);
  end
  
  % Expand the table with combined events.
  selTableRun = expandSelTableComboEvents(selTableRun, params.selParam);
  
  % For rewards, if rewardAbsent is present, use it, otherwise use the
  % other comparison.
  if all(selTableRun.subSel_rewardAbsent_pVal == 1)
    selTableRun.subSel_rewardCombo_selInd = selTableRun.subSel_reward_selInd;
  else
    selTableRun.subSel_rewardCombo_selInd = selTableRun.subSel_rewardAbsent_selInd;
  end
  % Add processed table to the array
  processedSelTableArray{run_i} = selTableRun;
  
end

spikePathBank.selTable = processedSelTableArray;

end
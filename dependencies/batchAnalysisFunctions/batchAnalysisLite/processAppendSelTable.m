function spikePathBank = processAppendSelTable(spikePathBank, params)
% A function which loads the selTable variable, processes it for use, and
% appends it back to the spikePathBank.

% Collect unit selectivity
[selTablePerRun, anovaTablePerRun, anovaBinParams, epochCatsParamsPerRun, epochStatsParamsPerRun, psthParamsPerRun] = spikePathLoad(spikePathBank, {'selTable', 'anovaTable', 'anovaBinParams', 'epochSWparams', 'epochTargParams', 'psthParams'}, params.spikePathLoadParams);
alpha = params.selParam.alpha;
mctMethod = params.selParam.mctmethod;
strThres = params.selParam.stretchThreshold;

% Variables to exclude
varExcludeanova = {'slidingWin_socIntTest_Eye'};
% varExcludesel = {'_cohensD'};
switch2StimWhole = true;            % Convert the sliding window test from parameters defined in the initial run to those defined below.

tables2Process = find(~cellfun('isempty', selTablePerRun))';

mergeAndReport = False;
if mergeAndReport
  % Combine all the tables
  % Generate variables for unifying the columns across all the runs.
  [columnNames, columnTypes] = deal([]);
  for run_i = tables2Process
    columnNames = [columnNames, selTablePerRun{run_i}.Properties.VariableNames];
    columnTypes = [columnTypes, varfun(@class, selTablePerRun{run_i},'OutputFormat','cell')];
  end

  [~, col_idx] = unique(columnNames);
  col_idx = sort(col_idx);
  uniqueCol = columnNames(col_idx);
  columnTypes = columnTypes(col_idx);

  selTableAll = [];
  for run_i = tables2Process

    % Pull out the table for the run

    % Add in missing columns with nans.
    [newCol, col_idx] = setdiff(uniqueCol, selTableRun.Properties.VariableNames);
    colDataTypes = columnTypes(col_idx);
    for col_i = 1:length(newCol)
      switch colDataTypes{col_i}
        case 'double'
          selTableRun{:, newCol(col_i)} = nan;
        case 'cell'
          selTableRun{:, newCol(col_i)} = {''};
        case 'string'
          selTableRun{:, newCol(col_i)} = '';
      end
    end

    selTableAll = [selTableAll; selTableRun];

  end

  % Once the selTable is combined run through generating the desired 'selInd'
  % columns.
  selTableUnits = selTableAll(contains(selTableAll.unitType, digitsPattern), :);

  % Variables - remove those not being considered for use. 
  pValVars = selTableUnits.Properties.VariableNames(contains(selTableUnits.Properties.VariableNames', '_pVal'))';
  % keepIdx = ~contains(pValVars, {'subSel_reward', 'epochSel_categories', 'saccDir', 'saccSel_categories', 'subSel_stim', 'subSel_fix', 'blinks', 'saccades', 'epochSel_socVNonSoc_Fix_pVal'});
  keepIdx = ~contains(pValVars, {'saccDir', 'subSel_stim', 'subSel_fix', 'blinks', 'saccades', 'epochSel_socVNonSoc_Fix_pVal', '_vBase_pVal', 'epochSel_categories'});
  pValVars = pValVars(keepIdx);

  selIndexValues = generateSelectivityIndex(selTableRun{:, pValVars}, alpha, mctMethod);

  % Summary stats
  TM_Idx = contains(pValVars, 'baseV_');
  TM_Idx = contains(pValVars, 'epochSel_socVNonSoc_') & contains(pValVars, '_vBase_pVal');
  tmIdx = any(selIndexValues(:, TM_Idx'),2);
  compT = horzcat(pValVars, num2cell(100*sum(selIndexValues(tmIdx,:))/sum(tmIdx))')
  
end

for run_i = tables2Process
  
  selTableRun = selTablePerRun{run_i};
  anovaTableRun = anovaTablePerRun{run_i};
  epochCParamsRun = epochCatsParamsPerRun{run_i};
  epochSParamsRun = epochStatsParamsPerRun{run_i};
  psthParamsRun = psthParamsPerRun{run_i};
  anovaParamsRun = anovaBinParams{run_i};
  
  unitCountRun = size(selTableRun,1);
  
  % Remove undesired variables
%   var2Remove = contains(selTableRun.Properties.VariableNames, varExcludesel);
%   selTableRun(:, var2Remove) = [];
    
  % Collect pVals and convert them to selectivity indicies
  pValIndex = contains(selTableRun.Properties.VariableNames, '_pVal') & ~contains(selTableRun.Properties.VariableNames, 'baseV_');
  pValVars = selTableRun.Properties.VariableNames(pValIndex)';
  
  % We don't care about all the tests run. Filter the pValVars variable for
  % the tests we actually want to report.
  keepIdx = ~contains(pValVars, {'saccDir', 'subSel_stim', 'subSel_fix', 'blinks', 'saccades', 'epochSel_socVNonSoc_Fix_pVal', '_vBase_pVal', 'epochSel_categories'});
  pValVars = pValVars(keepIdx);
  
  % Create and store the selectivity indicies.
  selIndexVars = strcat(extractBefore(pValVars, '_pVal'), '_selInd');
  selIndexValues = generateSelectivityIndex(selTableRun{:, pValVars}, alpha, mctMethod);
  selTableRun{:, selIndexVars} = selIndexValues;   % Add them back to the table.

  % Process the ANOVA table
  anovaVars = anovaTableRun.Properties.VariableNames';
  var2Remove = contains(anovaVars, varExcludeanova);
  anovaTableRun(:, var2Remove) = []; 
  
  % Check for pVal rows in the ANOVA table
  pValIndex = contains(anovaVars, '_pVal');
  pValVars = anovaVars(pValIndex)';
  selIndexVars = strcat(extractBefore(pValVars, '_pVal'), '_selInd');
  selIndexValues = zeros(unitCountRun, length(pValVars));
  
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
  selTableRun = convertSigRuns2Inds(anovaTableRun, selTableRun, params);
  
  % Create a task modulated index
  selTabVars = contains(selTableRun.Properties.VariableNames, 'baseV_') & contains(selTableRun.Properties.VariableNames, '_pVal');
  pValTaskMod = selTableRun{:, selTabVars};
  taskModInd = any(pValTaskMod <= alpha, 2);
  selTableRun.taskModulated_selInd = taskModInd;

  % Create a task modulated stim index
  selTabVars = contains(selTableRun.Properties.VariableNames, 'baseV_') & contains(selTableRun.Properties.VariableNames, '_pVal') & contains(selTableRun.Properties.VariableNames, 'stim');
  pValTaskMod = selTableRun{:, selTabVars};
  taskModInd = any(pValTaskMod <= alpha, 2);
  selTableRun.taskModStim_selInd = taskModInd;

  % Sliding Window to Epoch conversion.
  if switch2StimWhole
%     epochSParamsRun.times = [-500 0; -800 0; 0 2800; 2800 3150];
%     epochSParamsRun.labels = {'preFix', 'Fix', 'Stimulus', 'Reward'};
    epochSParamsRun.times = [0 2800];
    epochSParamsRun.labels = {'Stimulus'};
  end
  
  selTableRun = slidingWindowToEpochs(anovaTableRun, selTableRun, anovaParamsRun, psthParamsRun, epochCParamsRun, epochSParamsRun, params);
  
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
    selTableRun.subSel_rewardCombo_cohensD = selTableRun.subSel_reward_cohensD;
  else
    selTableRun.subSel_rewardCombo_selInd = selTableRun.subSel_rewardAbsent_selInd;
    selTableRun.subSel_rewardCombo_cohensD = selTableRun.subSel_rewardAbsent_cohensD;
  end
  
  % Add a label to distinguish Mo neurons from Sam Neurons easily
  if contains(spikePathBank.Row{run_i}, 'Mo')
    selTableRun.M1unit_selInd = true(size(selTableRun,1),1);
  else
    selTableRun.M1unit_selInd = false(size(selTableRun,1),1);
  end
  
  % Add processed table to the output array
  processedSelTableArray{run_i} = selTableRun;
  
end

% Combine tables 
selTableAll = [];
for run_i = 1:length(processedSelTableArray)
  selTableAll = [selTableAll; processedSelTableArray{run_i}];
end

% for each of the required comparisons, collect the correct values for
% comparisons, and the method for MCT.
fprintf('spikePathBank reduces to only runs with');

pVal_idx = contains(selTableAll.Properties.VariableNames, 'pVal');
pVal_idx(1:7) = 1;
selTableUnits = selTableAll(contains(selTableAll.unitType, digitsPattern),pVal_idx);

% Create lists of groups to be compared. the total number of comparisons
% here (columns * units) will be used for corrections.

end

function selIndMat = generateSelectivityIndex(pValMat, alpha, mctMethod)
% Function which takes in a vector or matrix of p values (shape n_neurons * n_tests), which need to be
% corrected for multiple comparisons. Returns a binary array of the same size and
% shape, indicating whether a null hypothesis is rejected.

% Testing done per neuron.

selIndMat = zeros(size(pValMat));

for unit_i = 1:size(pValMat,1)
  unit_tests = pValMat(unit_i,:);
  sorted_tests = sort(unit_tests);
  m = sum(~isnan(sorted_tests));
  
  switch mctMethod
    case 'fdr'
        % Benjamini Hochberg technique (Intro to Stat Learning Ed2. Pg 573)
        j = 1:m;
        thresholdVal = (alpha.* j)./m;

        compVec = sorted_tests(1:m) < thresholdVal;
        pValIdx = find(compVec, 1, 'last');
        pValThreshold = sorted_tests(pValIdx);

    case 'fwe'
      % Holm's Step-Down Procedure (ISL Ed2, Pg 566)
      % Calculate adjusted p-values (q-values)
      qvalues = zeros(1, m);
      for j = 1:m
          qvalues(j) = (alpha / (m + 1 - j));
      end

      compVec = sorted_tests(1:m) < qvalues;
      pValIdx = find(compVec,1, 'last');
      pValThreshold = sorted_tests(pValIdx);
      
     case 'bon'
      pValThreshold = alpha/m;
  end
  
  % if some value is below the threshold, compare the rest of vector
  if ~isempty(pValThreshold)
    selIndMat(unit_i, :) = unit_tests <= pValThreshold;
  end
  
end

selIndMat = logical(selIndMat);

end

function selTableRun = slidingWindowToEpochs(anovaTableRun, selTableRun, anovaParamsRun, psthParamsRun, epochCParamsRun, epochSParamsRun, params)
% Function which takes in the grid of neurons*bins and turns it into
% epochs*bin.
% Input
% - neurons*bins logical array (2D or 3D, output will match)
% - bin start times.
% - epoch windows (epoch * 2, start and end).

BaseArrayTests = {{'chasing', 'fighting', 'goalDirected', 'grooming', 'idle', 'mounting', 'objects', 'scene'}, {'goalDirected', 'idle', 'objects', 'scene', 'socialInteraction'}};
BaseArrayTestNames = {'categoriesTest', 'broadCatTest'};
% BaseArrayResort = {{'chasing', 'fighting', 'mounting', 'grooming', 'goalDirected', 'idle', 'objects', 'scene'}, {'socialInteraction', 'goalDirected', 'idle', 'objects'}};
% Resorted in another function

strThres = params.selParam.stretchThreshold;
objectStretches = params.selParam.objectStretches;   

% Generate the start times and end times for the bins. 
binSize = epochCParamsRun.binSize;
binStarts = anovaParamsRun.binStarts;
binStarts(1) = binStarts(1) + 1;
binEnds = anovaParamsRun.binEnds;

% Covert to epochWindows
epochNames = epochSParamsRun.labels;
epochTimes = epochSParamsRun.times;

% Remove pre stimulus epoch data.
epochNames = epochNames(epochTimes(:,1) >= -800);
epochTimes = epochTimes(epochTimes(:,1) >= -800, :);

% Shift the end times for epochs so bins which extend more than 75 ms into
% the next epoch are classified as such
epochTimes(2:end, 1) = epochTimes(2:end, 1) - binSize/2;
epochTimes(1:end-1, 2) = epochTimes(1:end-1, 2) - binSize/2;
epochTimes(end,end) = epochTimes(end,end) + 100;

% for each epoch, see which bins fall within.
bin2Epoch = zeros(size(binStarts));
for epoch_i = 1:length(epochNames)
  epochBins = epochTimes(epoch_i, 1) < binStarts & epochTimes(epoch_i, 2) >= binStarts;
  bin2Epoch(epochBins) = deal(epoch_i);
end

% For each test performed, cycle through, creating the necessary variable names 
% and check if there is a significant stretch within, and store its index
prefStimInd = find(contains(anovaTableRun.Properties.VariableNames, '_prefStim'));
prefStimName = anovaTableRun.Properties.VariableNames(prefStimInd);
testNames = extractBetween(prefStimName, '_', '_');

for test_i = 1:length(testNames)
  
  testInd = strcmp(testNames{test_i}, BaseArrayTestNames);
  BaseArray = BaseArrayTests{testInd};
  testData = anovaTableRun{:, prefStimInd(test_i)};  
  testData = vertcat(testData{:});
  [testDataIndicies, A] = cellArray2Indicies(testData);
  testEpochNames = strcat('SW_', testNames{test_i}, '_', epochNames);
    
  for epoch_i = 1:length(epochNames)
    
    testEpochIndicies = testDataIndicies(:, bin2Epoch == epoch_i, :);
    
    sigStretchPerObj = squeeze(findStretchRuns(testEpochIndicies) >= strThres);

    % Find run lengths for each object
    if objectStretches
      testEpochSigObjAny = any(sigStretchPerObj,2);
    else
      testEpochSigObjAny = findStretchRuns(any(testEpochIndicies,3)) >= strThres;
    end
    
    testEpochRunLen = squeeze(findStretchRuns(testEpochIndicies)) >= strThres;
    
    if size(testEpochRunLen,2) ~= length(BaseArray)
      % To make sure run sizes always match, expand cells in cases where
      % not all selectivities are represented.
      testEpochRunLenTmp = zeros(size(testEpochRunLen,1), length(BaseArray));
      [~, B] = intersect(BaseArray, A);
      testEpochRunLenTmp(:, B) = testEpochRunLen;
      testEpochRunLen = testEpochRunLenTmp;
      
    end
    
    % Condence into a structure which can be stored per unit.
    objIndPerUnit = mat2cell(testEpochRunLen, ones(size(testEpochRunLen,1), 1), size(testEpochRunLen,2), ones(size(testEpochRunLen,3), 1));

    % Add the index for units with significant stretches + the objects.
    selTableRun.([testEpochNames{epoch_i} '_sigStretch']) = testEpochSigObjAny;
    selTableRun.([testEpochNames{epoch_i} '_sigObj']) = objIndPerUnit;
    
  end
  
end

end

function selTableRun = convertSigRuns2Inds(anovaTableRun, selTableRun, params)
% Finds significant stretches of object selectivity, as defined in
% anovaTable, and adds new columns to selTable based on those results.

strThres = params.selParam.stretchThreshold;
objectStretches = params.selParam.objectStretches;   

% Find the column indicies and their names
prefStimInd = find(contains(anovaTableRun.Properties.VariableNames, '_prefStim'));
prefStimName = anovaTableRun.Properties.VariableNames(prefStimInd);

% Create the new column names
selIndexVars = strcat(extractBefore(prefStimName, '_prefStim'), '_selInd');
preStimObjVec = strcat(extractBefore(prefStimName, '_prefStim'), '_prefStimObj');

% Iterate across each instances of preferred stimuli column and collapse it
% into the selInd and the prefStimObj columns.
for p_val_ind = 1:length(prefStimInd)
  prefStimArray = anovaTableRun{:, prefStimInd(p_val_ind)};
  dataStack = vertcat(prefStimArray{:});
  
  [dataStackFinalTmp, uniqueStimLabels] = cellArray2Indicies(dataStack);
  
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

end

function selTableRun = replaceChanNum_1123(selTableRun)
% Run 20201123Mo001 - 004 has 32 channels, labeled Ch33 - Ch64. This causes
% problems. 

runList = selTableRun.dateSubj;

if any(strcmp(runList, "20201123Mo"))
  runInd2Change = strcmp(runList, "20201123Mo");
  table2Change = selTableRun(runInd2Change, :);
  % Ch33 - Ch64 --> Ch1 - Ch32
  table2Change.channel = strcat("Ch", string(double(extractAfter(table2Change.channel, "Ch"))-32));
  selTableRun(runInd2Change, :) = table2Change;
end

end

function [indicies, labels] = cellArray2Indicies(cellArray)
% Turns a 2D cell array of entries into a 3D Logical matrix of entries, one
% layer per unique entry in the cell array.

labels = unique(cellArray(:));
labels = labels(~contains(labels, 'None'));

indicies = zeros([size(cellArray), length(labels)]);
for ii = 1:length(labels)
  indicies(:,:,ii) = strcmp(cellArray, labels{ii}) * ii;
end

end
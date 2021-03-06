function spikePathBank = processAppendSelTable(spikePathBank, params)
% A function which loads the selTable variable, processes it for use, and
% appends it back to the spikePathBank.

% Collect unit selectivity
[selTablePerRun, anovaTablePerRun, anovaBinParams, epochCatsParamsPerRun, epochStatsParamsPerRun, psthParamsPerRun] = spikePathLoad(spikePathBank, {'selTable', 'anovaTable', 'anovaBinParams', 'epochSWparams', 'epochTargParams', 'psthParams'}, params.spikePathLoadParams);
alpha = params.selParam.alpha;
strThres = params.selParam.stretchThreshold;

% Keep the tables seperate, since they have to be returned to
% spikePathBank.
processedSelTableArray = cell(size(selTablePerRun));

% Variables to exclude
varExcludeanova = {'slidingWin_socIntTest_Eye'};
% varExcludesel = {'_cohensD'};

for run_i = 1:length(selTablePerRun)
  
  % Pull out the table for the run
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
  
  % Check for pVal rows.
  pValIndex = contains(selTableRun.Properties.VariableNames, '_pVal');
  pValVars = selTableRun.Properties.VariableNames(pValIndex)';
  selIndexVars = strcat(extractBefore(pValVars, '_pVal'), '_selInd');
  
  % Collect pVals and convert them to selectivity indicies, based on
  % previously defined alpha.
  pValValues = selTableRun{:, pValVars};
  selIndexValues = pValValues < alpha;
  
  % Add them back to the table.
  selTableRun{:, selIndexVars} = selIndexValues;
   
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
  
  % Sliding Window to Epoch conversion.
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
  else
    selTableRun.subSel_rewardCombo_selInd = selTableRun.subSel_rewardAbsent_selInd;
  end
  
  % Add processed table to the output array
  processedSelTableArray{run_i} = selTableRun;
  
end

spikePathBank.selTable = processedSelTableArray;

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
  
  BaseArray = BaseArrayTests{strcmp(testNames{test_i}, BaseArrayTestNames)};
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

function selTableRun = convertSigRuns2Inds(anovaTableRun, selTableRun, params);
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
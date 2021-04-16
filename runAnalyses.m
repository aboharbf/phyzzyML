function [analysisOutFilename] = runAnalyses(inputs)
% runAnalyses should be the main site for customization
% - inputs is a struct, which will be unpacked
%  this version of runAnalyses does the following:
%   - makes psth's and evoked potentials for (possibly overlapping)categories, 
%     specified at stimParamsFilename
%   - RF analyses: spike rate, evoked power, mean spike latency
%   - performs time-frequency and spike-field coherence analyses for each
%     channel, using the multitaper method implemented in Chronux
%   - performs across channel coupling analyses: spike-spike, field-field,
%     and spike-field, using Chronux
%   - makes tuning curves for parameterized images or categories, as
%     specified in the stim param file
%% unpack inputs
%List of inputs to be unpacked
inputFields = fields(inputs);

%Cycle through and unpack inputs
for input_ind = 1:length(inputFields)
  eval(sprintf('%s = inputs.%s;', inputFields{input_ind}, inputFields{input_ind}));
end

clear inputs
clear spikesByCategoryForTF
clear spikesByEventForTF
clear taskDataAll
clear analogInByCategory

%%% Load parameters from analysis and stim param files.
load(analysisParamFilename);
eventLabelsTmp = eventLabels; %note: hack to avoid overwriting list of not presented stimuli
refStructTmp = refStruct;
load(stimParamsFilename);
eventLabels = eventLabelsTmp; % conclusion of hack
refStruct = refStructTmp;

if ~plotSwitch.neuroGLM
  clear lfpByEvent
end

figStruct.channelUnitNames = channelUnitNames;
figStruct.figDir = outDir;
figStruct.figTag = sprintf('%s,Run%s', dateSubject, runNum);

channelNames = ephysParams.channelNames;
spikeChannels = ephysParams.spikeChannels;
lfpChannels = ephysParams.lfpChannels;

analogInChannels = analogInParams.analogInChannels;
analogInChannelNames = analogInParams.channelNames;
analogInChannelUnits = analogInParams.channelUnits;

psthPre = psthParams.psthPre;
ITI = psthParams.ITI;
psthPost = psthParams.psthPost;
smoothingWidth = psthParams.smoothingWidth;
psthErrorType = psthParams.errorType;
psthErrorRangeZ = psthParams.errorRangeZ;
psthBootstrapSamples = psthParams.bootstrapSamples;
psthColormapFilename = psthParams.psthColormapFilename;

lfpPreAlign = lfpAlignParams.msPreAlign; 
lfpPostAlign = lfpAlignParams.msPostAlign;
lfpPaddedBy = tfParams.movingWin(1)/2;

movingWin = tfParams.movingWin;
specgramRowAve = tfParams.specgramRowAve;
samPerMS = ephysParams.samPerMS;

frEpochs = zeros(length(frEpochsCell),2);
for epoch_i = 1:length(frEpochsCell)
  for range_i = 1:2
    tmp = frEpochsCell{epoch_i}{range_i};
    if isa(tmp,'function_handle')
      frEpochs(epoch_i,range_i) = tmp(psthImDur);
    else
      frEpochs(epoch_i,range_i) = tmp;
    end
  end
end

gridHole = taskData.gridHole;
recDepth = taskData.recDepth;
analysisOutFilename = strcat(outDir,'analyzedData.mat');
save(analysisOutFilename, 'dateSubject', 'runNum', 'analysisParamFilename', 'gridHole', 'recDepth'); % Task data is already in preprocessed, not sure why I'd want to save it here.

colors = {[0.55 0.13 0.16];[0.93 .2 0.15];[.98 0.65 0.13];[0 0.55 0.25];[0.15, 0.20, 0.5]};
chColors = [{'b'}, {[0 .6 0]} , {'m'}];

%% check calcSwitch definitions, set defaults, and apply field name updates as needed
calcSwitchFields = fields(calcSwitch);
% update deprecated fieldnames: syntax is:
% {{{'oldName1_1';'oldName1_2';...};'newName1'};{{'oldName2_1';'oldName2_2';...};'newName2'};...}
calcSwitchFieldnameUpdateMap = {{{'faceSelectIndex'};'selectivityIndices'};{{'faceSelectIndexEarly'};'selectivityIndicesEarly'};...
  {{'faceSelectIndexLate'};'selectivityIndicesLate'}};
for currentField_i = 1:length(calcSwitchFieldnameUpdateMap)
  for oldField_i = 1:length(calcSwitchFieldnameUpdateMap{currentField_i}{1})
    if isfield(calcSwitch,calcSwitchFieldnameUpdateMap{currentField_i}{1}{oldField_i})
      calcSwitch.(calcSwitchFieldnameUpdateMap{currentField_i}{2}) = calcSwitch.(calcSwitchFieldnameUpdateMap{currentField_i}{1}{oldField_i});
    end
  end
end

for field_i = 1:length(calcSwitchFields)
  if ~isfield(calcSwitch,calcSwitchFields{field_i}) || ~ismember(calcSwitch.(calcSwitchFields{field_i}),[0,1])
    calcSwitch.(calcSwitchFields{field_i}) = 0;
  end
end

% check analysisGroups definitions and set defaults
analysisGroupFields = fields(analysisGroups);
for field_i = 1:length(analysisGroupFields)
  field = analysisGroupFields{field_i};
  if ~isfield(analysisGroups,field) || ~isstruct(analysisGroups.(field)) || ~isfield(analysisGroups.(field),'groups')
    analysisGroups.(field).groups = {};
  end
end

% remove images and categories not presented from analysis groups
% assumes that either, (1) the groups field will contain a depth=1 cell array of
% stimulus/category names, or (2) the field groupDepth will exist and give
% the correct cell array depth of groups. Future version can check depth
% automatically. Also, asumes all fields other than groupDepth are cell
% arrays.

% Determine whether groups need to be removed.
analysisGroupFields = fieldnames(analysisGroups);
for field_i = 1:length(analysisGroupFields)
  field = analysisGroupFields{field_i};
  
  % If the field doesn't conform to this checker's syntax, or is more
  % complicated, then do not process it in this way.
  if isempty(analysisGroups.(field).groups)
    continue
  end
  if isfield(analysisGroups.(field),'groupDepth') && analysisGroups.(field).groupDepth > 1
    continue
  end
  
  subfields = fieldnames(analysisGroups.(field));
  itemDelimitedSubfields = cell(length(subfields),1); %we'll use this to remove all subfield entries corresponding to items not presented
  itemDelimitedSubfield_i = 0;
  for subfield_i = 1:length(itemDelimitedSubfields)
    subfield = subfields{subfield_i};
    if length(analysisGroups.(field).(subfield)) == length(analysisGroups.(field).groups) ...
        && length(analysisGroups.(field).(subfield){1}) == length(analysisGroups.(field).groups{1})... 
        && (iscell(analysisGroups.(field).(subfield){1}) || isnumeric(analysisGroups.(field).(subfield){1}))
      itemDelimitedSubfield_i = itemDelimitedSubfield_i + 1;
      itemDelimitedSubfields{itemDelimitedSubfield_i} = subfield;
    end
  end
  
  itemDelimitedSubfields = itemDelimitedSubfields(1:itemDelimitedSubfield_i);
  for group_i = 1:length(analysisGroups.(field).groups)
    newStruct = struct();
    for subfield_i = 1:length(itemDelimitedSubfields)
      if iscell(analysisGroups.(field).(itemDelimitedSubfields{subfield_i}){group_i})
        newStruct.(itemDelimitedSubfields{subfield_i}) = cell(size(analysisGroups.(field).groups{group_i}));
      else
        newStruct.(itemDelimitedSubfields{subfield_i}) = zeros(size(analysisGroups.(field).groups{group_i}));
      end
    end
    newItem_i = 0;
    for item_i = 1:length(analysisGroups.(field).groups{group_i})
      if ismember(analysisGroups.(field).groups{group_i}{item_i}, categoryList) || ismember(analysisGroups.(field).groups{group_i}{item_i}, eventLabels)
        newItem_i = newItem_i + 1;
        for subfield_i = 1:length(itemDelimitedSubfields)
          if iscell(analysisGroups.(field).(itemDelimitedSubfields{subfield_i}){group_i})
            newStruct.(itemDelimitedSubfields{subfield_i}){newItem_i} = analysisGroups.(field).(itemDelimitedSubfields{subfield_i}){group_i}{item_i};
          else
            newStruct.(itemDelimitedSubfields{subfield_i})(newItem_i) = analysisGroups.(field).(itemDelimitedSubfields{subfield_i}){group_i}(item_i);
          end
        end
      end
    end
    for subfield_i = 1:length(itemDelimitedSubfields)
      tmp = newStruct.(itemDelimitedSubfields{subfield_i});
      analysisGroups.(field).(itemDelimitedSubfields{subfield_i}){group_i} = tmp(1:newItem_i);
    end
  end
  
  % now, remove any groups that are empty after empty-item exclusion
  newStruct = analysisGroups.(field);
  validGroups = 0;
  for group_i = 1:length(analysisGroups.(field).groups)
    if ~isempty(analysisGroups.(field).groups{group_i})
      validGroups = validGroups + 1;
      for subfield_i = 1:length(subfields)
        newStruct.(subfields{subfield_i})(validGroups) = analysisGroups.(field).(subfields{subfield_i})(group_i);
      end
    else
      fprintf('The analysisGroups.%s group with index %d was removed because no item had a valid trial\n',field,group_i);
      disp(analysisGroups.(field).groups{group_i})
%       disp('If you think there should be valid trials, check that the item naming convention in the task log file matches that in stimulusParamFile.');
    end
  end
  for subfield_i = 1:length(subfields)
    newStruct.(subfields{subfield_i}) = newStruct.(subfields{subfield_i})(1:validGroups);
  end
  analysisGroups.(field) = newStruct;
end

for field_i = 1:length(analysisGroupFields)
  field = analysisGroupFields{field_i};
  if isempty(analysisGroups.(field).groups)
    continue
  end
  if (~isfield(analysisGroups.(field),'groupDepth')) || analysisGroups.(field).groupDepth == 1
    continue
  end
  if analysisGroups.(field).groupDepth > 2
    analysisGroups.(field).groups = {};
    continue
  end
  subfields = fieldnames(analysisGroups.(field));
  itemDelimitedSubfields = cell(size(subfields));
  itemDelimitedSubfield_i = 0;
  for subfield_i = 1:length(itemDelimitedSubfields)
    subfield = subfields{subfield_i};
    if strcmp(subfield, 'groupDepth')
      continue
    end
    itemDelimited = 1;
    for supergroup_i = 1:length(analysisGroups.(field).groups)
      if ~(length(analysisGroups.(field).(subfield){supergroup_i}) == length(analysisGroups.(field).groups{supergroup_i}) ...
          && iscell(analysisGroups.(field).(subfield){supergroup_i}) ...
          && length(analysisGroups.(field).(subfield){supergroup_i}{1}) == length(analysisGroups.(field).groups{supergroup_i}{1})...
          && (iscell(analysisGroups.(field).(subfield){supergroup_i}{1}) || isnumeric(analysisGroups.(field).(subfield){supergroup_i}{1})))
        itemDelimited = 0;
      end
    end
    if itemDelimited
      itemDelimitedSubfield_i = itemDelimitedSubfield_i + 1;
      itemDelimitedSubfields{itemDelimitedSubfield_i} = subfield;
    end
  end
  itemDelimitedSubfields = itemDelimitedSubfields(1:itemDelimitedSubfield_i);
  for supergroup_i = 1:length(analysisGroups.(field).groups)
    for group_i = 1:length(analysisGroups.(field).groups{supergroup_i})
      newStruct = struct();
      for subfield_i = 1:length(itemDelimitedSubfields)
        newStruct.(itemDelimitedSubfields{subfield_i}) = cell(size(analysisGroups.(field).groups{supergroup_i}{group_i}));
      end
      newItem_i = 0;
      for item_i = 1:length(analysisGroups.(field).groups{supergroup_i}{group_i})
        if ismember(analysisGroups.(field).groups{supergroup_i}{group_i}{item_i},categoryList) || ismember(analysisGroups.(field).groups{supergroup_i}{group_i}{item_i},eventLabels)
          newItem_i = newItem_i + 1;
          for subfield_i = 1:length(itemDelimitedSubfields)
            newStruct.(itemDelimitedSubfields{subfield_i}){newItem_i} = analysisGroups.(field).(itemDelimitedSubfields{subfield_i}){supergroup_i}{group_i}{item_i};
          end
        end
      end
      for subfield_i = 1:length(itemDelimitedSubfields)
        tmp = newStruct.(itemDelimitedSubfields{subfield_i});
        analysisGroups.(field).(itemDelimitedSubfields{subfield_i}){supergroup_i}{group_i} = tmp(1:newItem_i);
      end
    end
  end
end

% build category index map
catInds = struct();
for cat_i = 1:length(categoryList)
  catInds.(categoryList{cat_i}) = cat_i;
end

imInds = struct();
for image_i = 1:length(eventLabels)
  imInds.(eventLabels{image_i}) = image_i;
end
save(analysisOutFilename,'catInds','imInds', '-append');

catIndStruct = struct();
catIndMat = zeros(length(eventIDs), length(categoryList));
for cat_i = 1:length(categoryList)
  catIndMat(:, cat_i) = cellfun(@(x) any(strcmp(categoryList(cat_i), x)), eventCategories);
end
catIndStruct.catIndMat = logical(catIndMat);
catIndStruct.eventIDs = eventIDs;
catIndStruct.categoryList = categoryList;

% Generate groupImagesByLabel
groupLabelsByImage = zeros(length(eventLabels), length(analysisGroups.stimulusLabelGroups.groups));
groupLabelColorsByImage = ones(length(eventLabels),3,length(analysisGroups.stimulusLabelGroups));
for group_i = 1:1:length(analysisGroups.stimulusLabelGroups.groups)
  group = analysisGroups.stimulusLabelGroups.groups{group_i};
  groupColors = analysisGroups.stimulusLabelGroups.colors{group_i};
  if iscell(groupColors)
    colorArray = zeros(length(groupColors),3);
    for item_i = 1:length(groupColors)
      % this silly line converts from colorspec letters to the corresponding rgb values
      % Silly line now checks to see if conversion is needed.
      if ischar(groupColors{item_i})
        colorArray(item_i,:) = rem(floor((strfind('kbgcrmyw', groupColors{item_i}) - 1) * [0.25 0.5 1]), 2);
      else
        colorArray(item_i,:) = groupColors{item_i};
      end
    end
    groupColors = colorArray;
    analysisGroups.stimulusLabelGroups.colors{group_i} = groupColors;
  end
  
  if group_i == 1
    %Maybe just make the above colors into cells?
    %Swap colors defined above to numbers
    colors_tmp = zeros(length(colors),3);
    for item_i = 1:length(colors)
      if ischar(colors{item_i})
        colors_tmp(item_i,:) = rem(floor((strfind('kbgcrmyw', colors{item_i}) - 1) * [0.25 0.5 1]), 2);
      else
        colors_tmp(item_i,:) = colors{item_i};
      end
    end
    colors = colors_tmp;
  end
  
  for image_i = 1:length(eventLabels)
    for item_i = 1:length(group)
      if any(strcmp(eventCategories{image_i},group{item_i})) || strcmp(eventLabels{image_i},group{item_i})
        groupLabelsByImage(image_i,group_i) = item_i;
        groupLabelColorsByImage(image_i,:,group_i) = groupColors(item_i,:);
        break
      end
    end
    if groupLabelsByImage(image_i,group_i) == 0
      Output.VERBOSE(sprintf('no stim category match found for %s\n',eventLabels{image_i}));
    end
  end
end

save(analysisOutFilename,'groupLabelsByImage','groupLabelColorsByImage', '-append');

%% Eye Signal Processing
eyeStatsParams.eventIDs = eventIDs;
eyeStatsParams.eventLabels = eventLabels;
eyeDataStruct = struct();

if plotSwitch.pupilDilation
  eyeDataStruct = pupilDilation(analogInByEvent, eyeStatsParams, eyeDataStruct, catIndStruct, figStruct);
end

if plotSwitch.eyeStatsAnalysis
  [eyeDataStruct, eyeBehStatsByStim, eyeInByEvent] = eyeStatsAnalysis(analogInByEvent, eyeStatsParams, taskData, eyeDataStruct);
  save(analysisOutFilename, 'eyeBehStatsByStim', 'eyeInByEvent', '-append');
else
  % Reshapes analogInByEvent into eye signal. not smoothed
  eyeInByEvent = cellfun(@(n) squeeze(n(:,1:2,:,lfpPaddedBy+1:length(n)-(lfpPaddedBy+1))), analogInByEvent, 'UniformOutput', false);
  eyeBehStatsByStim = [];
  save(analysisOutFilename, 'eyeInByEvent', '-append');
end

if plotSwitch.attendedObject
  % See whats being attended
  eyeDataStruct = calcEyeObjectTrace(eyeInByEvent, channelUnitNames, psthParams, eventIDs, taskData, eyeDataStruct);

  % Assign each Saccade labeled a target.
  if plotSwitch.eyeStatsAnalysis
    eyeBehStatsByStim = saccadeTarget(eyeDataStruct, eyeBehStatsByStim);
  end
  
end

save(analysisOutFilename, 'eyeBehStatsByStim', 'eyeDataStruct','-append');

if plotSwitch.eyeCorrelogram 
  eyeCorrelogram(eyeInByEvent, psthParams, eventLabels, figStruct)
end

%% Spike Processing

if ~calcSwitch.spikeTimes %use 1 ms bins for spikes
  spikesByEventBinned = calcSpikeTimes(spikesByEvent, psthParams);
  spikesByCategoryBinned = calcSpikeTimes(spikesByCategory, psthParams);
  save(analysisOutFilename,'spikesByEventBinned', '-append');
  analysisOutFilenameBig = fullfile(fileparts(analysisOutFilename), 'analyzedDataBig.mat');
  save(analysisOutFilenameBig,'spikesByCategoryBinned', '-v7.3');
end

if calcSwitch.imagePSTH && calcSwitch.spikeTimes
  [psthByImage, psthErrByImage] = calcStimPSTH(spikesByEvent, psthEmptyByEvent, calcSwitch.spikeTimes, psthParams, spikeAlignParams);
  save(analysisOutFilename,'psthByImage','psthErrByImage', '-append');
elseif calcSwitch.imagePSTH
  [psthByImage, psthErrByImage] = calcStimPSTH(spikesByEventBinned, psthEmptyByEvent, calcSwitch.spikeTimes, psthParams, spikeAlignParams);
  save(analysisOutFilename,'psthByImage','psthErrByImage', '-append');
end

if calcSwitch.categoryPSTH && calcSwitch.spikeTimes
  [psthByCategory, psthErrByCategory] = calcStimPSTH(spikesByCategory, psthEmptyByCategory, calcSwitch.spikeTimes, psthParams, spikeAlignParams);
  save(analysisOutFilename,'psthByCategory','psthErrByCategory', '-append');
elseif calcSwitch.categoryPSTH
  [psthByCategory, psthErrByCategory] = calcStimPSTH(spikesByCategoryBinned, psthEmptyByCategory, calcSwitch.spikeTimes, psthParams, spikeAlignParams);
  save(analysisOutFilename,'psthByCategory','psthErrByCategory', '-append');
  clear spikesByCategoryBinned
end

%% Spike, Task data, Eye Data Analysis

if plotSwitch.eyeStimOverlay
  [eyeInByEventDS] = eyeStimOverlay(eyeInByEvent, eyeDataStruct, spikesByEventBinned, eventIDs, taskData, psthParams, eyeStimOverlayParams);
  save(analysisOutFilename,'eyeInByEventDS', '-append');
end

% Initalize a table which the 3 sensitivity functions will use.
selTable = initializeSelTable(channelUnitNames, channelNames, dateSubject, runNum, gridHole, recDepth);

% Determine selectivity for events labeled in eventData, + blinks & rewards.
if plotSwitch.subEventAnalysis
  ephysParams.spikeTimes = calcSwitch.spikeTimes;
  ephysParams.channelUnitNames = channelUnitNames;
  subEventAnalysisParams.onsetsByEvent = onsetsByEvent;
  subEventAnalysisParams.eventIDs = eventIDs;
    
  [subEventSigStruct, specSubEventStruct, selTable] = subEventAnalysis(eyeBehStatsByStim, spikesByChannel, taskData, ephysParams, subEventAnalysisParams, selTable, figStruct);
  save(analysisOutFilename,'subEventSigStruct', 'specSubEventStruct','-append');
%   clear subEventSigStruct spikesByChannel
end

epochStatsParams.groupLabelsByImage = groupLabelsByImage;
selTable = saccadePerUnit(spikesByEventBinned, eyeBehStatsByStim, psthParams, analysisGroups.stimulusLabelGroups, eventIDs, refStruct, epochStatsParams, selTable);

% Determine which units are selective for saccades (direction).
selTable = saccadeSel(spikesByEventBinned, eyeBehStatsByStim, psthParams.psthPre, selTable);

% Tests between epochs of the conditions
selTable = epochStats(spikesByEvent, selTable, eventIDs, taskData.paradigm, epochStatsParams);

if ~strcmp(taskData.paradigm, 'familiarFace')
  selTable = epochCats(spikesByEventBinned, eyeDataStruct.saccadeByStim, selTable, eventIDs, taskData.paradigm, psthParams, epochCatsParams);
end

save(analysisOutFilename, 'selTable', '-append');
% error('Done with selTable')

% NeuroGLM time
% if plotSwitch.neuroGLM
%   neuroGLMStruct = runNeuroGLM(spikesByEvent, lfpByEvent, taskData, trialIDsByEvent, catIndStruct, eyeDataStruct, eyeBehStatsByStim, eyeInByEvent, neuroGLMParams);
% end

%% Plotting and further analyses

if plotSwitch.imagePsth
  for channel_i = 1:length(channelNames)
    for unit_i = 1:length(channelUnitNames{channel_i})
      
      % if no isolated unit defined, plot just MUA, (unsorted is identical)
      if length(channelUnitNames{channel_i}) == 2 && unit_i == 1
        continue;
      end
      
      % figure out stimOrder
      if isfield(psthParams, 'sortStim') && psthParams.sortStim
        if ~exist('NewStimOrder')
          sortOrder = psthParams.(taskData.paradigm).sortOrder;
          groupLabelsByImagePSTH = zeros(length(eventIDs),1);
          for image_i = 1:length(eventLabels)
            for item_i = 1:length(sortOrder)
              if any(strcmp(eventCategories{image_i},sortOrder{item_i})) || strcmp(eventLabels{image_i}, sortOrder{item_i})
                groupLabelsByImagePSTH(image_i) = item_i;
                break
              end
            end
            if groupLabelsByImagePSTH(image_i) == 0
              Output.VERBOSE(sprintf('no stim category match found for %s\n',eventLabels{image_i}));
            end
          end
          %Sorts based on group membership.
          [~, NewStimOrder] = sort(groupLabelsByImagePSTH);
        end
      else
        NewStimOrder = 1:length(eventLabels);
      end
      
      % If you want to Z-score the data before it runs, do so here.
      psthTitle = sprintf('Per Image PSTH %s, %s',channelNames{channel_i}, channelUnitNames{channel_i}{unit_i});
      if isfield(psthParams,'Zscore') && psthParams.Zscore
        %Calculate Baseline in PSTH by taking the average of all bins in the Pre-stimulus period
        psthTitle = [psthTitle ' - Z scored'];
        for ii = 1:size(psthByImage{channel_i}{unit_i},1)
          if psthParams.Zscore == 1
            preStimBaseline = mean(psthByImage{channel_i}{unit_i}(ii,1:ITI)); %(Baseline) Mean Subtraction, during pre-fixation period (pre-successful fixation).
            psthByImage{channel_i}{unit_i}(ii,:) = (psthByImage{channel_i}{unit_i}(ii,:) - preStimBaseline)/std(psthByImage{channel_i}{unit_i}(ii,:)); %Z score.
          else
            stimBaseline = mean(psthByImage{channel_i}{unit_i}(ii,:)); % Mean Subtraction.
            psthByImage{channel_i}{unit_i}(ii,:) = (psthByImage{channel_i}{unit_i}(ii,:) - stimBaseline)/std(psthByImage{channel_i}{unit_i}(ii,:)); %Z score.
          end
        end
      end
      
      % Plot 
      figure('Name', psthTitle, 'NumberTitle','off','units','normalized','outerposition',figStruct.figPos);
      eventLabelsPlot = strrep(eventLabels, '_', ' ');
      plotPSTH(psthByImage{channel_i}{unit_i}(NewStimOrder,:), [], [], psthParams, 'color', psthTitle, eventLabelsPlot(NewStimOrder));
      clear figData
      title(psthTitle);
      figData.z = psthByImage{channel_i}{unit_i};
      figData.x = -psthPre:psthImDur+psthPost;
      saveFigure(outDir, sprintf('imPSTH_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}, runNum), figData, figStruct, figStruct.figTag);
      
    end
  end
end

% If you want to Z-score Catagory data before it runs, do so here.
if psthParams.Zscore
  %Calculate Baseline in PSTH by taking the average of all bins in the Pre-stimulus period
  zTag = ' Z scored';
  for channel_i = 1:length(channelNames)
    for unit_i = 1:length(channelUnitNames{channel_i})
      for ii = 1:size(psthByCategory{channel_i}{unit_i},1)
        if psthParams.Zscore == 1
          preStimBaseline = mean(psthByCategory{channel_i}{unit_i}(ii,1:ITI)); %(Baseline) Mean Subtraction, during pre-fixation period (pre-successful fixation).
          psthByCategory{channel_i}{unit_i}(ii,:) = (psthByCategory{channel_i}{unit_i}(ii,:) - preStimBaseline)/std(psthByCategory{channel_i}{unit_i}(ii,:)); %Z score.
        else
          stimBaseline = mean(psthByCategory{channel_i}{unit_i}(ii,:)); % Mean Subtraction.
          psthByCategory{channel_i}{unit_i}(ii,:) = (psthByCategory{channel_i}{unit_i}(ii,:) - stimBaseline)/std(psthByCategory{channel_i}{unit_i}(ii,:)); %Z score.
        end
      end
    end
  end
else
  zTag = '';
end

if plotSwitch.categoryPsth
  categoryListPlot = strrep(categoryList, '_', ' ');
  for channel_i = 1:length(channelNames)
    for unit_i = 1:length(channelUnitNames{channel_i})
      
      %if no isolated unit defined, plot just MUA, not also unsorted (since it's identical)
      if length(channelUnitNames{channel_i}) == 2 && unit_i == 1 
        continue;
      end
           
      psthTitle = sprintf('Per Catagory PSTH %s, %s, %s',channelNames{channel_i}, channelUnitNames{channel_i}{unit_i}, zTag);
      figure('Name',psthTitle,'NumberTitle','off','units','normalized', 'outerposition', figStruct.figPos);
      plotPSTH(psthByCategory{channel_i}{unit_i}, [],  [], psthParams, 'color', psthTitle, categoryListPlot);
      clear figData
      title(psthTitle);
      figData.z = psthByCategory{channel_i}{unit_i};
      figData.x = -psthPre:psthImDur+psthPost;
      saveFigure(outDir, sprintf('catPSTH_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, figStruct, figStruct.figTag);
    end
  end
end

if plotSwitch.analysisGroupsPsth && isfield(analysisGroups.analysisGroupPSTH, taskData.paradigm)
  % Identify which PSTHes to stack.
%   analysisPSTH = fields(analysisGroups.analysisGroupPSTH);
  analysisPSTHStruct = analysisGroups.analysisGroupPSTH.(taskData.paradigm);
  analysisPSTH = fields(analysisPSTHStruct);
  categoryListPlot = strrep(categoryList, '_', ' ');
  
  for ii = 1:length(analysisPSTH)
    stimOrder = [];
    categoryListSet = analysisPSTHStruct.(analysisPSTH{ii});
    
    % Pick the members which belong in this group
    for jj = 1:length(categoryListSet)
      stimOrder = [stimOrder; find(strcmp(categoryList, categoryListSet(jj))')];
    end
    
    if ~isempty(stimOrder)
      
      if sum(stimOrder) > 10
        psthType = 'color';
      else
        psthType = 'line';
      end
      
      
      for channel_i = 1:length(channelNames)
        for unit_i = 1:length(channelUnitNames{channel_i})
          if length(channelUnitNames{channel_i}) == 2 && unit_i == 1 %if no isolated unit defined, plot just MUA, not also unsorted (since it's identical)
            continue;
          end
          
          psthTitle = sprintf('Analysis Group %s PSTH %s, %s, %s', analysisPSTH{ii}, channelNames{channel_i}, channelUnitNames{channel_i}{unit_i}, zTag);
          figure('Name',psthTitle,'NumberTitle','off','units','normalized','outerposition',figStruct.figPos);

          plotPSTH(psthByCategory{channel_i}{unit_i}(stimOrder,:), psthErrByCategory{channel_i}{unit_i}(stimOrder,:),  [], psthParams, psthType, psthTitle, categoryListPlot(stimOrder));

          clear figData
          title(psthTitle);
          figData.z = psthByCategory{channel_i}{unit_i};
          figData.x = -psthPre:psthImDur+psthPost;
          saveFigure(outDir, sprintf('analysisGroupPSTH_%s_%s_%s_Run%s', analysisPSTH{ii}, channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, figStruct, figStruct.figTag);
        end
      end
    end
  end
  
end

[spikeCountsByImageByEpoch, firingRatesByImageByEpoch, firingRateErrsByImageByEpoch] = deal(cell(size(frEpochs,1),1));
for epoch_i = 1:size(frEpochs,1)
  [spikeCounts, fr, frErr] = spikeCounter(spikesByEvent, frEpochs(epoch_i,1), frEpochs(epoch_i,2));
  spikeCountsByImageByEpoch{epoch_i} = spikeCounts;
  firingRatesByImageByEpoch{epoch_i} = fr;
  firingRateErrsByImageByEpoch{epoch_i} = frErr;
end

if ~isempty(spikesByCategory)
  [firingRatesByCategoryByEpoch, firingRateErrsByCategoryByEpoch, spikeCountsByCategoryByEpoch] = deal(cell(size(frEpochs,1),1));
  for epoch_i = 1:size(frEpochs,1)
    [spikeCounts, fr, frErr] = spikeCounter(spikesByCategory, frEpochs(epoch_i,1), frEpochs(epoch_i,2));
    firingRatesByCategoryByEpoch{epoch_i} = fr;
    firingRateErrsByCategoryByEpoch{epoch_i} = frErr;
    spikeCountsByCategoryByEpoch{epoch_i} = spikeCounts;
  end

%   catFr = firingRatesByCategoryByEpoch{1};
%   catFrErr = firingRateErrsByCategoryByEpoch{1};
%   catSpikeCounts = spikeCountsByCategoryByEpoch{1};  
  trialCountsByCategory = cellfun(@(x) length(x{1}{1}), spikesByCategory);
  
  save(analysisOutFilename,'frEpochs','firingRatesByCategoryByEpoch','firingRateErrsByCategoryByEpoch','spikeCountsByCategoryByEpoch','trialCountsByCategory', '-append');
   
end

trialCountsByImage = cellfun(@(x) length(x{1}{1}), spikesByEvent);

save(analysisOutFilename, 'firingRatesByImageByEpoch', 'firingRateErrsByImageByEpoch', 'spikeCountsByImageByEpoch', '-append');

%% Unit Selectivity Analyses

[sigStruct, imageSortOrder, nullModelPvalues, nullTraceMeans, nullTraceSD] = nullModelStats(psthByImage, spikeCountsByImageByEpoch, firingRatesByImageByEpoch, firingRateErrsByImageByEpoch, ...
  trialCountsByImage, analysisGroups, epochLabels, eventIDs, ephysParams);

save(analysisOutFilename,'sigStruct', 'selTable', '-append');

for epoch_i = 1:length(firingRatesByImageByEpoch)
  epochTag = sprintf('%d-%d ms',frEpochs(epoch_i,1), frEpochs(epoch_i,2));
  % preferred images
  for channel_i = 1:length(spikeChannels)
    for unit_i = 1:length(channelUnitNames{channel_i})
      chanUnitTag = sprintf('%s, %s', channelNames{channel_i}, channelUnitNames{channel_i}{unit_i});
      if length(channelUnitNames{channel_i}) == 2 && unit_i == 1 %if no isolated unit defined, plot just MUA, not also unsorted (since it's identical)
        continue;
      end
      
      imageSortOrderInd = imageSortOrder{epoch_i}{channel_i}{unit_i};
      imageSortedRates  = firingRatesByImageByEpoch{epoch_i}{channel_i}(unit_i,imageSortOrderInd);  %todo: write the firing rates to file
      if sum(imageSortedRates) == 0
        imageSortedRates(1) = 1; %Defense against this causing problems later.
      end
      imFrErrSorted = firingRateErrsByImageByEpoch{epoch_i}{channel_i}(unit_i,imageSortOrderInd);
      sortedImageLabels = eventLabels(imageSortOrderInd);
      sortedEventIDs = eventIDs(imageSortOrderInd);
      trialCountsByImageSorted = trialCountsByImage(imageSortOrderInd);
      if exist('groupLabelColorsByImage','var')
        sortedGroupLabelColors = groupLabelColorsByImage(imageSortOrderInd,:,:);
      else
        sortedGroupLabelColors = ones(length(eventLabels),3);
      end
      if epoch_i == 1
%         fprintf('\n\n\nPreferred Images: %s, %s\n\n',channelNames{channel_i},channelUnitNames{channel_i}{unit_i});
%         for i = 1:min(10,length(eventLabels))
%           fprintf('%d) %s: %.2f +/- %.2f Hz\n',i,sortedImageLabels{i},imageSortedRates(i),imFrErrSorted(i));
%         end
%         fprintf('\nLeast Preferred Images: %s, %s\n\n',channelNames{channel_i},channelUnitNames{channel_i}{unit_i});
%         for i = 0:min(4,length(eventLabels)-1)
%           fprintf('%d) %s: %.2f +/- %.2f Hz\n',i,sortedImageLabels{end-i},imageSortedRates(end-i), imFrErrSorted(end-i));
%         end
        % Preferred images raster plot
        if isfield(plotSwitch,'prefImRaster') && plotSwitch.prefImRaster
          if isfield(plotSwitch,'topStimToPlot')
            topStimToPlot = plotSwitch.topStimToPlot;
          else
            topStimToPlot = 5;
          end
          prefImRasterTitle = sprintf('Preferred Image Raster - %s, %s',chanUnitTag, epochTag);
          fh = figure('Name',prefImRasterTitle,'NumberTitle','off','units','normalized','outerposition',figStruct.figPos);
          clear figData
          figData.z = spikesByEvent(imageSortOrderInd(1:topStimToPlot));
          figData.x = -psthPre:psthImDur+psthPost;
          raster(spikesByEvent(imageSortOrderInd(1:topStimToPlot)), sortedImageLabels(1:topStimToPlot), psthParams, stimTiming.ISI, channel_i, unit_i, colors);
          title(sprintf('Preferred Images, %s %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
          saveFigure(outDir, sprintf('prefImRaster_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, figStruct, figStruct.figTag );
        end
        % Preferred Image Raster plots, Color coded
        if isfield(plotSwitch,'prefImRasterColorCoded') && plotSwitch.prefImRasterColorCoded
          colorTagList = {'ColorCodedSpikes', 'AttendedObjShaded', 'SaccadeShaded', 'PupilDilationShaded'};
          colorTag = sprintf(colorTagList{plotSwitch.prefImRasterColorCoded});
          if isfield(plotSwitch,'topStimToPlot') && plotSwitch.topStimToPlot == 0
            topStimToPlot = length(eventLabels);
          elseif isfield(plotSwitch,'topStimToPlot')
            topStimToPlot = plotSwitch.topStimToPlot;
          else
            topStimToPlot = 5;
          end
          prefImRasterTitle = sprintf('Preferred Image Raster, Color coded - %s, %s',chanUnitTag, epochTag);
          fh = figure('Name',prefImRasterTitle,'NumberTitle','off','units','normalized','outerposition',figStruct.figPos);
          clear figData
          figData.z = spikesByEvent(imageSortOrderInd(1:topStimToPlot));
          figData.x = -psthPre:psthImDur+psthPost;
          rasterColorCoded(fh, spikesByEvent(imageSortOrderInd(1:topStimToPlot)), sortedEventIDs(1:topStimToPlot), psthParams, stimTiming.ISI, channel_i, unit_i, eyeDataStruct,  plotSwitch.prefImRasterColorCoded);
          title(sprintf('Preferred Images, %s %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
          saveFigure(outDir, sprintf('prefImRaster_%s_%s_%s_Run%s',colorTag, chanUnitTag,epochTag,runNum), figData, figStruct, figStruct.figTag );
        end
        % preferred images raster-evoked overlay
        if isfield(plotSwitch,'prefImRasterEvokedOverlay') && plotSwitch.prefImRasterEvokedOverlay
          prefImRasterEvokedOverlayTitle = sprintf('Preferred Image Raster, Evoked Potential Overlay - %s, %s',chanUnitTag, epochTag);
          fh = figure('Name',prefImRasterEvokedOverlayTitle,'NumberTitle','off');
          rasterEvoked(spikesByEvent(imageSortOrderInd), lfpByEvent(imageSortOrderInd), sortedImageLabels, psthPre, psthPost, psthImDur, stimTiming.ISI, lfpPaddedBy, channel_i, colors, 1)
          title(sprintf('Preferred Images, from top, %s %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
          saveFigure(outDir, sprintf('prefImRaster-LFP_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, figStruct, figStruct.figTag );
        end
        % preferred images average evoked
        if isfield(plotSwitch,'prefImAverageEvoked') && plotSwitch.prefImRasterAverageEvokedOverlay
          prefImRasterAverageEvokedOverlayTitle = sprintf('Preferred Image Raster, Average Evoked Potential Overlay - %s, %s',chanUnitTag, epochTag);
          fh = figure('Name',prefImRasterAverageEvokedOverlayTitle,'NumberTitle','off');
          averageEvoked(spikesByEvent(imageSortOrderInd), lfpByEvent(imageSortOrderInd), sortedImageLabels, psthPre, psthPost, psthImDur, stimTiming.ISI, lfpPaddedBy, channel_i, colors)
          title(sprintf('Preferred Images - Average Evoked, from top, %s %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
          saveFigure(outDir, sprintf('prefImAverage-LFP_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, figStruct, figStruct.figTag );
        end
        % preferred images raster-evoked overlay, with other channels
        if isfield(plotSwitch,'prefImMultiChRasterEvokedOverlay') && plotSwitch.prefImMultiChRasterEvokedOverlay
          prefImMultiChRasterEvokedOverlayTitle = sprintf('Preferred Image Raster, Multichannel, Evoked Potential Overlay - %s, %s',chanUnitTag, epochTag);
          fh = figure('Name',prefImMultiChRasterEvokedOverlayTitle,'NumberTitle','off');
          rasterEvokedMultiCh(spikesByEvent(imageSortOrderInd), lfpByEvent(imageSortOrderInd), sortedImageLabels, psthPre, psthPost, psthImDur, stimTiming.ISI, lfpPaddedBy, 1:length(lfpChannels), channelNames, colors)
          title(sprintf('Preferred Images, from top, %s %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
          saveFigure(outDir, sprintf('prefImRaster-LFP-MultiChannel_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum), figData, figStruct, figStruct.figTag );
        end
      end
      % Image preference barplot
      if isfield(plotSwitch,'imageTuningSorted') && plotSwitch.imageTuningSorted
        for group_i = 1:length(analysisGroups.stimulusLabelGroups.groups)
          imageTuningSortedTitle = sprintf('Image Tuning, Sorted - %s, %s',chanUnitTag, epochTag);
          fh = figure('Name',imageTuningSortedTitle,'NumberTitle','off');
          groupName = analysisGroups.stimulusLabelGroups.names{group_i};
          
          % Plot Bars
          [HB, HE, ~, ~, ~] = superbar(imageSortedRates,'E',imFrErrSorted,'P',nullModelPvalues{epoch_i}{channel_i}{unit_i}, 'PStarShowNS', 0,'BarFaceColor',sortedGroupLabelColors(:,:,group_i), 'ErrorbarLineWidth', 1.5, 'ErrorbarColor', [0 0 .1]);
          set(gca,'XTickLabel',sortedImageLabels,'XTickLabelRotation',45,'XTick',1:length(eventLabels),'TickDir','out');
          
          % Add Null model Line
          lineProps.col = {'k'};
          lineProps.width = 1;
          h = mseb(1:length(imageSortedRates),nullTraceMeans{epoch_i}{channel_i}{unit_i}, nullTraceSD{epoch_i}{channel_i}{unit_i},lineProps);
          h.patch.FaceAlpha = '0.5';
          
          % Reshape plot
          ylim(max(ylim(),0));
          ylabel('Firing rate [Hz]');
          xlabel('Stimulus Video Title');
          ylim([ylim()*1.2]) %provides additional room for legend.
          title(imageTuningSortedTitle);
          % Legends are made, not born
          [~, firstgroupMember] = unique(groupLabelsByImage(imageSortOrderInd)); %Get the right handles to put in the legend
          legendHandleArray = [HB(firstgroupMember); h.mainLine; HE(firstgroupMember(end))];
          legendHandleLabel = [analysisGroups.stimulusLabelGroups.groups{group_i}(:)' {'Null model'} {sprintf('SE Bars, n = %d', mode(trialCountsByImageSorted))}]';
          legend(legendHandleArray,legendHandleLabel,'FontSize',7,'NumColumns', 2);
          legend('boxoff');
          
          % Save the Figure
          clear figData
          figData.y = imageSortedRates;
          figData.e = imFrErrSorted;
          saveFigure(outDir, sprintf('imageTuningSorted, %s %s %s R%s',groupName, chanUnitTag, epochTag,runNum), figData, figStruct, figStruct.figTag );
        end
      end
    end
  end
  
  % multi-channel MUA image preference
  %(todo: all these, but with mean-sigma as the sort criterion)
  if length(channelNames) > 1
    multiChSpikesMin = zeros(length(eventLabels));
    multiChSpikesMinNorm = zeros(length(eventLabels));
    multiChSpikesMeanNorm = zeros(length(eventLabels));
    multiChMua = zeros(length(spikeChannels),length(eventLabels));
    multiChMuaNorm = zeros(length(spikeChannels),length(eventLabels));
    for channel_i = 1:length(spikeChannels) %todo: deal with MUA vs. units
      imFr = firingRatesByImageByEpoch{epoch_i};
      multiChMua(channel_i,:) = imFr{channel_i}(end,:);
      multiChMuaNorm(channel_i,:) = imFr{channel_i}(end,:)/max(imFr{channel_i}(end,:));
    end
    multiChSpikesMin = min(multiChMua);
    multiChSpikesMinNorm = min(multiChMuaNorm);
    multiChSpikesMeanNorm = mean(multiChMuaNorm);
    % multi-channel preferred images
    [imageSortedRatesMulti, imageSortOrderMulti] = sort(multiChSpikesMin,2,'descend');
    
    sortedImageLabels = eventLabels(imageSortOrderMulti);
    multiChMinFrSort.imageSortedRates = imageSortedRatesMulti';
    multiChMinFrSort.imageSortOrder = imageSortOrder';
    multiChMinFrSort.sortedImageLabels = sortedImageLabels';
    save(analysisOutFilename,'multiChMinFrSort','-append');
    
    fprintf('\n\n\nMulti-channel Preferred Images, Maximin\n\n');
    for i = 1:min(10,length(eventLabels))
      fprintf('%d) %s: %.2f Hz\n',i,sortedImageLabels{i},imageSortedRatesMulti(i));
    end
    fprintf('\n\nMulti-channel Least preferred Images, Maximin\n\n');
    for i = 0:min(4,length(eventLabels)-1)
      fprintf('%d) %s: %.2f Hz \n',i,sortedImageLabels{end-i},imageSortedRatesMulti(end-i));
    end
    
    [imageSortedRatesMulti, imageSortOrderMulti] = sort(multiChSpikesMinNorm,2,'descend');
    sortedImageLabels = eventLabels(imageSortOrderMulti);
    multiChMinNormFrSort.imageSortedRates = imageSortedRatesMulti';
    multiChMinNormFrSort.imageSortOrder = imageSortOrderMulti';
    multiChMinNormFrSort.sortedImageLabels = sortedImageLabels';
    save(analysisOutFilename,'multiChMinNormFrSort','-append');
    fprintf('\n\n\nMulti-channel Preferred Images, Channel-Normalized Maximin\n\n');
    for i = 1:min(10, length(sortedImageLabels))
      fprintf('%d) %s: %.2f Hz\n',i,sortedImageLabels{i},imageSortedRatesMulti(i));
    end
    fprintf('\n\nMulti-channel Least preferred Images, Channel-Normalized Maximin\n\n');
    for i = 0:min(4, length(sortedImageLabels)-1)
      fprintf('%d) %s: %.2f Hz \n',i,sortedImageLabels{end-i},imageSortedRatesMulti(end-i));
    end
    
    [imageSortedRatesMulti, imageSortOrderMulti] = sort(multiChSpikesMeanNorm,2,'descend');
    sortedImageLabels = eventLabels(imageSortOrderMulti);
    multiChMeanNormFrSort.imageSortedRates = imageSortedRatesMulti';
    multiChMeanNormFrSort.imageSortOrder = imageSortOrderMulti';
    multiChMeanNormFrSort.sortedImageLabels = sortedImageLabels';
    save(analysisOutFilename,'multiChMeanNormFrSort','-append');
    
    fprintf('\n\n\nMulti-channel Preferred Images, Channel-Normalized Mean\n\n');
    for i = 1:min(10,length(eventLabels))
      fprintf('%d) %s: %.2f Hz\n',i,sortedImageLabels{i},imageSortedRatesMulti(i));
    end
    fprintf('\n\nMulti-channel Least preferred Images, Channel-Normalized Mean\n\n');
    for i = 0:min(4,length(eventLabels)-1)
      fprintf('%d) %s: %.2f Hz \n',i,sortedImageLabels{end-i},imageSortedRatesMulti(end-i));
    end
  end
  % selectivity indices
  calcSwitches = [calcSwitch.selectivityIndices, calcSwitch.selectivityIndicesEarly, calcSwitch.selectivityIndicesLate];
  selectivityIndices = struct();
  %todo: add parametric and nonparametric z scores
  for calc_i = 1:length(calcSwitches)
    if ~calcSwitches(calc_i)
      continue
    end
    catFr = firingRatesByCategoryByEpoch{calc_i};
    imFr = firingRatesByImageByEpoch{calc_i};
    frCalcOn = frEpochs(calc_i,1);
    frCalcOff = frEpochs(calc_i,2);
    if calcSwitch.faceSelectIndex
      for channel_i = 1:length(spikeChannels)
        for unit_i = 1:length(channelUnitNames{channel_i})
          fprintf('\n\n%s %s Preference Indices, %dms - %d ms\n',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},frCalcOn,frCalcOff)
          for group_i = 1:length(analysisGroups.selectivityIndex.groups)
            group = analysisGroups.selectivityIndex.groups{group_i};
            if length(group) ~= 2
              Output.INFO(sprintf('Selectivity index only implemented for item pairs. Group %s contains %d items. Skipping...',group_i,length(group)));
              continue
            end
            if isfield(catInds,group{1})
              response1 = catFr{channel_i}(unit_i,catInds.(group{1}));
            else
              response1 = imFr{channel_i}(unit_i,imInds.(group{1}));
            end
            if isfield(catInds,group{2})
              response2 = catFr{channel_i}(unit_i,catInds.(group{2}));
            else
              response2 = imFr{channel_i}(unit_i,imInds.(group{2}));
            end
            selectIndex = (response1 - response2) / (response1 + response2);
            fprintf('%s > %s: %.3f\n',group{1},group{2},selectIndex);
            selectivityIndices.(analysisGroups.selectivityIndex.names{group_i}).(sprintf('from%dms_to_%dms',frEpochs(calc_i,1),frEpochs(calc_i,2))) = selectIndex;
          end
        end
      end
    end
  end
  save(analysisOutFilename,'selectivityIndices','-append');
  
  % category preference bar plot
  calcSwitches = [plotSwitch.stimPrefBarPlot, plotSwitch.stimPrefBarPlotEarly, plotSwitch.stimPrefBarPlotLate];
  for calc_i = 1:length(calcSwitches)
    if ~calcSwitches(calc_i)
      continue
    end
    catFr = firingRatesByCategoryByEpoch{calc_i};
    imFr = firingRatesByImageByEpoch{calc_i};
    catFrErr = firingRateErrsByCategoryByEpoch{calc_i};
    imFrErr = firingRateErrsByImageByEpoch{calc_i};
    frCalcOn = frEpochs(calc_i,1);
    frCalcOff = frEpochs(calc_i,2);
    for group_i = 1:length(analysisGroups.stimPrefBarPlot.groups)
      group = analysisGroups.stimPrefBarPlot.groups{group_i};
      groupName = analysisGroups.stimPrefBarPlot.names{group_i};
      for channel_i = 1:length(spikeChannels)
        fh = figure('Name','Catagory Preference Bar plot','NumberTitle','off');
        for subgroup_i = 1:length(group)
          itemInds = zeros(length(group{subgroup_i}),1);
          itemNames = cell(length(group{subgroup_i}),1);
          itemTypes = cell(length(group{subgroup_i}),1);
          responses = zeros(length(group{subgroup_i}),1);
          responseErrs = zeros(length(group{subgroup_i}),1);
          for item_i = 1:length(group{subgroup_i})
            if isfield(catInds,group{subgroup_i}{item_i})
              itemInds(item_i) = catInds.(group{subgroup_i}{item_i});
              itemNames{item_i} = categoryList{itemInds(item_i)};
              itemTypes{item_i} = 'cat';
            else
              itemInds(item_i) = imInds.(group{subgroup_i}{item_i});
              itemNames{item_i} = eventLabels{itemInds(item_i)};
              itemTypes{item_i} = 'im';
            end
          end
          for unit_i = 1:length(channelUnitNames{channel_i})
            if length(channelUnitNames{channel_i}) == 2
              if unit_i == 1
                continue
              else
                ax = subplot(length(group),1,subgroup_i);
              end
            else
              ax = subplot(length(group), length(channelUnitNames{channel_i}), unit_i+length(channelUnitNames{channel_i})*(subgroup_i-1));
            end
            for item_i = 1:length(itemInds)
              if strcmp(itemTypes{item_i},'cat')
                responses(item_i) = catFr{channel_i}(unit_i,itemInds(item_i));
                responseErrs(item_i) = catFrErr{channel_i}(unit_i,itemInds(item_i));
              else
                responses(item_i) = imFr{channel_i}(unit_i,itemInds(item_i));
                responseErrs(item_i) = imFrErr{channel_i}(unit_i,itemInds(item_i));
              end
            end
            superbar(ax, responses,'E',responseErrs,'BarFaceColor',analysisGroups.stimPrefBarPlot.colors{group_i}{subgroup_i});
            set(gca,'XTickLabel',itemNames,'XTickLabelRotation',45,'XTick',1:length(itemNames),'TickDir','out');
            tmp = ylim();
            ylim([0, tmp(2)]);
            ylabel('Firing rate, Hz');
            title(channelUnitNames{channel_i}{unit_i});
          end
        end
        suptitle(sprintf('%s spikes category tuning, %dms-%dms post-onset',channelNames{channel_i},frCalcOn,frCalcOff));
        saveFigure(outDir, sprintf('catPrefBars_%s_%s_%dms-%dms_Run%s',groupName,channelNames{channel_i},frCalcOn,frCalcOff,runNum), [], figStruct, figStruct.figTag );
      end
    end
  end
  
  % tuning curves
  calcSwitches = [isfield(plotSwitch,'tuningCurves') && plotSwitch.tuningCurves, isfield(plotSwitch,'tuningCurvesEarly') && plotSwitch.tuningCurvesEarly, ...
    isfield(plotSwitch,'tuningCurvesLate') && plotSwitch.tuningCurvesLate];
  %todo: compute and save measures like pref view, tuning depth, fisher info, fwhm, (mirror) symmetry
  for calc_i = 1:length(calcSwitches)
    if ~calcSwitches(calc_i)
      continue
    end
    catFr = firingRatesByCategoryByEpoch{calc_i};
    imFr = firingRatesByImageByEpoch{calc_i};
    catFrErr = firingRateErrsByImageByEpoch{calc_i};
    imFrErr = firingRateErrsByCategoryByEpoch{calc_i};
    frCalcOn = frEpochs(calc_i,1);
    frCalcOff = frEpochs(calc_i,2);
    for group_i = 1:length(analysisGroups.tuningCurves.groups)
      group = analysisGroups.tuningCurves.groups{group_i};
      groupName = analysisGroups.tuningCurves.names{group_i};
      muaTcFig = figure();
      for channel_i = 1:length(channelNames)
        channelTcFig = figure();
        for unit_i = 1:length(channelUnitNames{channel_i})
          tcFrs = zeros(length(group),1);
          tcFrErrs = zeros(length(group),1);
          for item_i = 1:length(group)
            if isfield(catInds,group{item_i})
              tcFrs(item_i) = catFr{channel_i}(unit_i,catInds.(group{item_i}));
              tcFrErrs(item_i) = catFrErr{channel_i}(unit_i,catInds.(group{item_i}));
            else
              tcFrs(item_i) = imFr{channel_i}(unit_i,imInds.(group{item_i}));
              tcFrErrs(item_i) = imFrErr{channel_i}(unit_i,imInds.(group{item_i}));
            end
          end
          if length(channelUnitNames{channel_i}) == 2  % don't make the single-channel plot if no unit defined
            close(channelTcFig);
            break
          end
          subplot(1,length(channelUnitNames{channel_i}),unit_i);
          errorbar(analysisGroups.tuningCurves.paramValues{group_i},tcFrs,tcFrErrs,'linestyle','-','linewidth',4);
          xlabel(analysisGroups.tuningCurves.paramLabels{group_i});
          ylabel('firing rate (Hz)');
          title(channelUnitNames{channel_i}{unit_i});
        end
        if length(channelUnitNames{channel_i}) > 2
          suptitle(sprintf('%s, %s, %dms - %d ms post-onset',groupName, channelNames{channel_i}, frCalcOn, frCalcOff));
          clear figData
          figData.x = analysisGroups.tuningCurves.paramValues{group_i};
          figData.y = tcFrs;
          figData.e = tcFrErrs;
          drawnow;
          saveFigure(outDir, sprintf('%s_%s_%dms-%dms_Run%s',regexprep(groupName,' ',''),channelNames{channel_i},frCalcOn,frCalcOff,runNum), figData, figStruct, figStruct.figTag );
        end
        figure(muaTcFig);
        subplot(1,length(channelNames),channel_i);
        errorbar(analysisGroups.tuningCurves.paramValues{group_i},tcFrs,tcFrErrs,'linestyle','-','linewidth',4);
        xlabel(analysisGroups.tuningCurves.paramLabels{group_i});
        ylabel('firing rate (Hz)');
        title(sprintf('%s MUA',channelNames{channel_i}));
      end
      suptitle(sprintf('%s, %dms - %d ms post-onset',groupName, frCalcOn, frCalcOff));
      drawnow;
      saveFigure(outDir, sprintf('%s_%dms-%dms_Run%s',regexprep(groupName,' ',''),frCalcOn,frCalcOff,runNum), figData, figStruct, figStruct.figTag );
    end
  end
end

if isfield(plotSwitch, 'stimPSTHoverlay') && plotSwitch.stimPSTHoverlay
  % (psthByImage, sortMask, inclusionMask, epochLabels, stimDir, psthParams, ephysParams, lfpPaddedBy, taskEventList, outDir)
  stimPSTHoverlay(psthByImage, imageSortOrder, nullModelPvalues, stimDir, psthParams, lfpPaddedBy, eventIDs, outDir);
end

if isfield(plotSwitch, 'spikeCorr') && plotSwitch.spikeCorr
  spikeCorrTraces = spikeCorr(spikesByEventBinned);
  save(analysisOutFilename,'spikeCorrTraces','-append');
end

% firing rate RF color plot
if taskData.RFmap
  load AkinoriColormap.mat;
  rfGrid  = unique(taskData.stimJumps,'rows');
  gridX = unique(rfGrid(:,1));
  gridsize = 2*mean(gridX(2:end,1)-gridX(1:end-1,1));
  % these are the x and y values at which to interpolate RF values
  xi = linspace(min(rfGrid(:,1)),max(rfGrid(:,1)),200);
  yi = linspace(min(rfGrid(:,2)),max(rfGrid(:,2)),200);
  calcSwitches = [isfield(plotSwitch,'RF') && plotSwitch.RF, isfield(plotSwitch,'rfEarly') && plotSwitch.rfEarly,...
    isfield(plotSwitch,'rfLate') && plotSwitch.rfLate];  
  for calc_i = 1:length(calcSwitches)
    if ~calcSwitches(calc_i)
      continue
    end
    frCalcOn = frEpochs(calc_i,1);
    frCalcOff = frEpochs(calc_i,2);
    for channel_i = 1:length(channelNames)
      for unit_i = 1:length(channelUnitNames{channel_i})
        if length(channelUnitNames{channel_i}) == 2  %skip to MUA if no isolated unit
          if unit_i == 1
            continue
          end
        end
        meanRF = zeros(length(rfGrid),1);
        allTrialSpikeCounts = cell(size(meanRF));
        for image_i = 1:length(eventLabels)
          imageRF = zeros(length(rfGrid),1);
          imageRFerrs = zeros(length(rfGrid),1);
          spikeLatencyRF = zeros(length(rfGrid),1);
          %also calc evoked peak time? power-weighted latency? peak of mean evoked correlogram?
          evokedPowerRF = zeros(length(rfGrid),1);
          for grid_i = 1:length(rfGrid)
            gridPointTrials = ismember(jumpsByImage{image_i},rfGrid(grid_i,:),'rows');
            if sum(gridPointTrials) == 0 
              imageRF(grid_i) = 0;
              disp('no data for this image and location');
              continue;
            end
            trialSpikes = spikesByEvent{image_i}{channel_i}{unit_i}(gridPointTrials);
            trialSpikeCounts = zeros(size(trialSpikes,1),1);
            totalSpikes = 0;
            spikeLatency = 0;
            evokedPotential = zeros(size(lfpByEvent{image_i}(1,channel_i,1,samPerMS*(frCalcOn+psthPre):samPerMS*(frCalcOff+psthPre)))); 
            gridPointTrialInds = 1:length(gridPointTrials);
            gridPointTrialInds = gridPointTrialInds(gridPointTrials);
            for trial_i = 1:length(trialSpikes)
              numSpikes = sum(trialSpikes(trial_i).times > frCalcOn & trialSpikes(trial_i).times < frCalcOff);
              trialSpikeCounts(trial_i) = numSpikes;
              totalSpikes = totalSpikes + numSpikes;
              % bad latency measure; try weighting spike times by 1/ISI
              spikeLatency = spikeLatency + mean(trialSpikes(trial_i).times(trialSpikes(trial_i).times > frCalcOn & trialSpikes(trial_i).times < frCalcOff));
              evokedPotential = evokedPotential + lfpByEvent{image_i}(1,channel_i,gridPointTrialInds(trial_i),samPerMS*(frCalcOn+psthPre):samPerMS*(frCalcOff+psthPre));
            end
            imageRF(grid_i) = (1000/(frCalcOff-frCalcOn))*totalSpikes/length(trialSpikes);
            imageRFerrs(grid_i) = (1000/(frCalcOff-frCalcOn))*std(trialSpikeCounts)/sqrt(length(trialSpikeCounts));
            spikeLatencyRF(grid_i) = 1/sum(gridPointTrials)*spikeLatency;
            if unit_i == length(channelUnitNames{channel_i})
              evokedPotential = 1/sum(gridPointTrials)*(evokedPotential - mean(evokedPotential));
              evokedPowerRF(grid_i) = sum(evokedPotential.^2);
            end
            if isempty(allTrialSpikeCounts{grid_i})
              allTrialSpikeCounts{grid_i} = trialSpikeCounts;
            else
              allTrialSpikeCounts{grid_i} = vertcat(allTrialSpikeCounts{grid_i},trialSpikeCounts);
            end
          end
          meanRF = meanRF + imageRF;
          display_map(rfGrid(:,1),rfGrid(:,2),imageRF,xi,yi,2.2857*gridsize,0,saveFig,sprintf('%s %s, %s RF',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},eventLabels{image_i}),...
            [outDir sprintf('RF_%s_%s_%s_Run%s.png',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},eventLabels{image_i},runNum)]);
          
          fh = plotRF(rfGrid, imageRF, imageRFerrs, 'heat', AkinoriColormap);
          drawnow;
%           saveFigure(outDir, sprintf('psthEvokedOverlay_%s_Run%s',groupName,runNum), figData, figStruct, figStruct.figTag );    
          
          if isfield(plotSwitch,'latencyRF') && plotSwitch.latencyRF
            display_map(rfGrid(:,1),rfGrid(:,2),spikeLatencyRF,xi,yi,2.2857*gridsize,0,saveFig,sprintf('%s %s, %s Latency RF',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},eventLabels{image_i}),...
              [outDir sprintf('LatencyRF_%s_%s_%s_Run%s.png',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},eventLabels{image_i},runNum)]);
          end
          if isfield(plotSwitch,'evokedPowerRF') && plotSwitch.evokedPowerRF
            display_map(rfGrid(:,1),rfGrid(:,2),evokedPowerRF,xi,yi,2.2857*gridsize,0,saveFig,sprintf('%s %s, %s Evoked Power RF',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},eventLabels{image_i}),...
              [outDir sprintf('EvokedPowerRF_%s_%s_%s_Run%s.png',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},eventLabels{image_i},runNum)]);
          end
          % todo: add background subtracted version
          % todo: add subplots version
          % todo: coherency RFs (cc, cpt, ptpt)
          % todo: granger RF
          % todo: bandpassed power RFs
          % todo: stimulus decodability RF
        end
        meanRF = meanRF/length(eventLabels);
        display_map(rfGrid(:,1),rfGrid(:,2),meanRF,xi,yi,2.2857*gridsize,0,saveFig,sprintf('%s %s, Mean RF',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}),...
          [outDir sprintf('MeanRF_%s_%s_Run%s.png',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},runNum)]);
        
        meanRFerrs = zeros(length(rfGrid),1);
        for grid_i = 1:length(rfGrid)
          meanRFerrs(grid_i) = std(allTrialSpikeCounts{grid_i})/sqrt(length(allTrialSpikeCounts{grid_i}));
        end
        fh = plotRF(rfGrid, meanRF, meanRFerrs, 'heat', AkinoriColormap);
        drawnow;
        saveFigure(outDir, sprintf('psthEvokedOverlay_%s_Run%s',groupName,runNum), figData, figStruct, figStruct.figTag ); 
      end
    end
  end
  disp('Receptive field mapping run; returning before coupling analyses');
  return
end

% evokedTF says to do TF analysis on the trial average psth and evoked potentials
if calcSwitch.meanEvokedTF && calcSwitch.spikeTimes
  % note: because I'm currently working with spike times rather than
  % binned spikes, I need to do this loop to make a struct whose 'times'
  % field has all the spikes from a category
  allSpikesByCategoryForTF = cell(length(spikesByCategory),1);
  for cat_i = 1:length(spikesByCategory)
    catChannelSpikes = cell(length(channelNames),1);
    for channel_i = 1:length(channelNames)
      channelUnitSpikes = cell(length(channelUnitNames{channel_i}),1);
      for unit_i = 1:length(channelUnitNames{channel_i})
        s.times = [];
        for trial_i = 1:length(spikesByCategoryForTF{cat_i}{channel_i}{unit_i})
          s.times = vertcat(s.times,spikesByCategoryForTF{cat_i}{channel_i}{unit_i}(trial_i).times);
        end
        channelUnitSpikes{unit_i} = s;
      end
      catChannelSpikes{channel_i} = channelUnitSpikes;
    end
    allSpikesByCategoryForTF{cat_i} = catChannelSpikes;
  end
end

% if we're going to need induced psth's and lfp's, build them now

%by event
if (calcSwitch.inducedSpectra || calcSwitch.inducedImageTF) && ~calcSwitch.spikeTimes && ~calcSwitch.inducedTrialMagnitudeCorrection
  % subtract the event mean psth from every trial to obtain induced psth
  spikesByEventBinnedInduced = spikesByEventBinned;
  for event_i = 1:length(spikesByEventBinned)
    for channel_i = 1:length(channelNames)
      for unit_i = 1:length(channelUnitNames{channel_i})
        spikesByEventBinnedInduced{event_i}{channel_i}{unit_i} = spikesByEventBinned{event_i}{channel_i}{unit_i} - mean(spikesByEventBinned{event_i}{channel_i}{unit_i},1); %note: matlab automatically applies subtraction to every row
      end
    end
  end
  lfpByEventInduced = lfpByEvent;
  for event_i = 1:length(lfpByEvent)
    lfpByEventInduced{event_i} = lfpByEvent{event_i} - mean(lfpByEvent{event_i},3); %note: matlab automatically applies subtraction to every trial/channel appropriately
  end
end

% by category - spikes
if (calcSwitch.inducedSpectra || calcSwitch.inducedCatSpikeTF || calcSwitch.inducedCoupling) && ~calcSwitch.spikeTimes && ~calcSwitch.inducedTrialMagnitudeCorrection
  % subtract the category mean psth from every trial to obtain induced psth
  spikesByCategoryBinnedInduced = spikesByCategoryBinned;
  for cat_i = 1:length(spikesByCategoryBinned)
    for channel_i = 1:length(channelNames)
      for unit_i = 1:length(channelUnitNames{channel_i})
        spikesByCategoryBinnedInduced{cat_i}{channel_i}{unit_i} = spikesByCategoryBinned{cat_i}{channel_i}{unit_i} - mean(spikesByCategoryBinned{cat_i}{channel_i}{unit_i},1); %note: matlab automatically applies subtraction to every row
      end
    end
  end
end

% by category - LFP (removed spikeTimes conditional since the LFP component
% seems divorced from that).
if (calcSwitch.inducedSpectra || calcSwitch.inducedCatLFPTF || calcSwitch.inducedCoupling) && ~calcSwitch.inducedTrialMagnitudeCorrection
  lfpByCategoryInduced = lfpByCategory;
  for cat_i = 1:length(lfpByCategory)
    lfpByCategoryInduced{cat_i} = lfpByCategory{cat_i} - mean(lfpByCategory{cat_i},3); %note: matlab automatically applies subtraction to every trial/channel appropriately
  end
end

% if we're going to need induced, response-magnitude-corrected psth's and lfp's, build them now
%by event
if (calcSwitch.inducedSpectra || calcSwitch.inducedImageTF) && ~calcSwitch.spikeTimes && calcSwitch.inducedTrialMagnitudeCorrection
  % subtract the category mean psth from every trial to obtain induced psth
  spikesByEventBinnedInduced = spikesByEventBinned;
  for event_i = 1:length(spikesByEventBinned)
    for channel_i = 1:length(channelNames)
      for unit_i = 1:length(channelUnitNames{channel_i})
        unitEventMean = mean(spikesByEventBinned{event_i}{channel_i}{unit_i},1);
        unitEventMeanSpikeBkgnd = sum(unitEventMean(1:50)); %todo: set this according to psth pre etc. to get best background estimate
        unitEventMeanBkgndSub = unitEventMean - unitEventMeanSpikeBkgnd;
        unitEventMeanTotalEvokedSpikes = sum(unitEventMeanBkgndSub);
        for trial_i = 1:size(spikesByEventBinned{event_i}{channel_i}{unit_i},1)
          %magnitudeCorrector = mean(spikesByEventBinned{event_i}{channel_i}{unit_i},1)*sum(spikesByEventBinned{event_i}{channel_i}{unit_i});
          trialSpikes = spikesByEventBinned{event_i}{channel_i}{unit_i}(trial_i,:);
          trialSpikesBkgndSub = trialSpikes - mean(trialSpikes(1:50));
          evokedMagnitudeCorrector =  sum(trialSpikesBkgndSub)/unitEventMeanTotalEvokedSpikes;
          spikesByEventBinnedInduced{event_i}{channel_i}{unit_i} = trialSpikesBkgndSub - evokedMagnitudeCorrector*unitEventMeanBkgndSub;
        end
      end
    end
  end
  lfpByEventInduced = lfpByEvent;
  for event_i = 1:length(lfpByEvent)
    for channel_i = 1:length(channelNames)
      channelEventMean = squeeze(mean(lfpByEvent{event_i}(1,channel_i,:,:),3));
      channelEventMeanMagnitude = mean(abs(channelEventMean));
      for trial_i = 1:size(lfpByEvent{event_i},3)
        magnitudeCorrector = mean(squeeze(abs(lfpByEvent{event_i}(1,channel_i,trial_i,:))))/channelEventMeanMagnitude;
        lfpByEventInduced{event_i}(1,channel_i,trial_i,:) = lfpByEvent{event_i}(1,channel_i,trial_i,:) - magnitudeCorrector*reshape(channelEventMean,1,1,1,[]);
      end
    end
  end
end

% by category - Spikes
if (calcSwitch.inducedSpectra || calcSwitch.inducedCatSpikeTF || calcSwitch.inducedCoupling) && ~calcSwitch.spikeTimes && calcSwitch.inducedTrialMagnitudeCorrection
  % subtract the category mean psth from every trial to obtain induced psth
  spikesByCategoryBinnedInduced = spikesByCategoryBinned;
  for cat_i = 1:length(spikesByCategoryBinned)
    for channel_i = 1:length(channelNames)
      for unit_i = 1:length(channelUnitNames{channel_i})
        unitCatMean = mean(spikesByCategoryBinned{cat_i}{channel_i}{unit_i},1);
        unitCatMeanSpikeBkgnd = sum(unitCatMean(1:50)); %todo: set this according to psth pre etc. to get best background estimate
        unitCatMeanBkgndSub = unitCatMean - unitCatMeanSpikeBkgnd;
        unitCatMeanTotalEvokedSpikes = sum(unitCatMeanBkgndSub);
        for trial_i = 1:size(spikesByCategoryBinned{cat_i}{channel_i}{unit_i},1)
          %magnitudeCorrector = mean(spikesByCategoryBinned{cat_i}{channel_i}{unit_i},1)*sum(spikesByCategoryBinned{cat_i}{channel_i}{unit_i});
          trialSpikes = spikesByCategoryBinned{cat_i}{channel_i}{unit_i}(trial_i,:);
          trialSpikesBkgndSub = trialSpikes - mean(trialSpikes(1:50));
          evokedMagnitudeCorrector =  sum(trialSpikesBkgndSub)/unitCatMeanTotalEvokedSpikes;
          spikesByCategoryBinnedInduced{cat_i}{channel_i}{unit_i} = trialSpikesBkgndSub - evokedMagnitudeCorrector*unitCatMeanBkgndSub;
        end
      end
    end
  end
end

% by category - LFP (Removed spikesTimes conditional).
if (calcSwitch.inducedSpectra || calcSwitch.inducedCatLFPTF || calcSwitch.inducedCoupling) && calcSwitch.inducedTrialMagnitudeCorrection
  lfpByCategoryInduced = lfpByCategory;
  for cat_i = 1:length(lfpByCategory)
    for channel_i = 1:length(channelNames)
      channelCatMean = squeeze(mean(lfpByCategory{cat_i}(1,channel_i,:,:),3));
      channelCatMeanMagnitude = mean(abs(channelCatMean));
      for trial_i = 1:size(lfpByCategory{cat_i},3)
        magnitudeCorrector = mean(squeeze(abs(lfpByCategory{cat_i}(1,channel_i,trial_i,:))))/channelCatMeanMagnitude;
        lfpByCategoryInduced{cat_i}(1,channel_i,trial_i,:) = lfpByCategory{cat_i}(1,channel_i,trial_i,:) - magnitudeCorrector*reshape(channelCatMean,1,1,1,[]);
      end
    end
  end
end

%%%% evoked potential plots
times = -psthPre:psthImDur+psthPost;
if isfield(plotSwitch,'evokedPsthMuaMultiCh') && plotSwitch.evokedPsthMuaMultiCh
  for group_i = 1:length(analysisGroups.evokedPsthOnePane.groups)
    group = analysisGroups.evokedPsthOnePane.groups{group_i};
    groupName = analysisGroups.evokedPsthOnePane.names{group_i};
    fh = figure();
    handlesForLegend = gobjects(2*length(lfpChannels),1);
    forLegend = cell(2*length(lfpChannels),1);
    for item_i = 1:length(group)
      subplot(length(group),1,item_i);
      title(group{item_i});
      xlabel('time (ms)');
      yyaxis right
      ylabel('lfp (normalized)'); 
      yyaxis left
      ylabel('firing rate (normalized)');
      hold on
      if isfield(catInds,group{item_i})
        lfpByItem = lfpByCategory;
        psthByItem = psthByCategory;
        itemNum = catInds.(group{item_i});
      else
        lfpByItem = lfpByEvent;
        psthByItem = psthByImage;
        itemNum = imInds.(group{item_i});
      end
      for channel_i = 1:length(lfpChannels)
        yyaxis right
        itemLfp = squeeze(mean(lfpByItem{itemNum}(1,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3));
        lfpHandle = plot(times,itemLfp/max(itemLfp),'color',chColors{channel_i}, 'linestyle','-','linewidth',2);  
        yyaxis left
        itemPSTH = psthByItem{channel_i}{end}(itemNum,:);
        spikeHandle = plot(times,itemPSTH/max(itemPSTH),'color',chColors{channel_i},'linestyle','--','linewidth',2); 
        handlesForLegend(2*channel_i-1) = lfpHandle;
        handlesForLegend(2*channel_i) = spikeHandle;
        forLegend{2*channel_i-1} = strcat(channelNames{channel_i},' LFP');
        forLegend{2*channel_i} = strcat(channelNames{channel_i},' MUA');
        if channel_i == length(channelNames)
          if item_i == 1
            legend(handlesForLegend,forLegend,'Location','northeastoutside');
            % align x axes with colorbar; see http://stackoverflow.com/questions/5259853/matlab-how-to-align-the-axes-of-subplots-when-one-of-them-contains-a-colorbar
            % note: drawnow line synchronizes the rendering thread with the main
            drawnow;
            pos1 = get(gca,'Position');
          else
            pos2 = get(gca,'Position');
            pos2(1) = pos1(1);
            pos2(3) = pos1(3);
            set(gca,'Position',pos2);
          end
        end
        h = get(gca,'ylim');
      end
      h = plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3,'linestyle','-');
      %this turns off legend for this line. Source: https://www.mathworks.com/matlabcentral/answers/406-how-do-i-skip-items-in-a-legend
      set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
      xlim([min(times), max(times)]);
    end
    drawnow;
    saveFigure(outDir, sprintf('psthEvokedOverlay_%s_Run%s',groupName,runNum), figData, figStruct, figStruct.figTag );
  end
end

%evoked potentials
for channel_i = 1:length(lfpChannels)
  if isfield(plotSwitch,'evokedByCategory') && plotSwitch.evokedByCategory
    for group_i = 1:length(analysisGroups.evokedPotentials.groups)
      group = analysisGroups.evokedPotentials.groups{group_i};
      groupName = analysisGroups.evokedPotentials.names{group_i};
      fh = figure();
      responses = zeros(length(group),length(times));
      responseErrs = zeros(length(group),length(times));
      for item_i = 1:length(group)
        if isfield(catInds,group{item_i})
          responses(item_i,:) = squeeze(mean(lfpByCategory{catInds.(group{item_i})}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3))';
          responseErrs(item_i,:) = squeeze(std(lfpByCategory{catInds.(group{item_i})}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),[],3)/sqrt(size(lfpByCategory{catInds.(group{item_i})},3)))';
        else
          responses(item_i,:) = squeeze(mean(lfpByEvent{imInds.(group{item_i})}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3))';
          responseErrs(item_i,:) = squeeze(std(lfpByEvent{imInds.(group{item_i})}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),[],3)/sqrt(size(lfpByEvent{imInds.(group{item_i})},3)))';
        end
      end
      lineProps.width = 3;
      lineProps.col = analysisGroups.evokedPotentials.colors{group_i};
      hold on
      mseb(repmat(times,length(group),1),responses,responseErrs,lineProps);
      legend(group);
      h = get(gca,'ylim');
      h = plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
      %this turns off legend for this line. Source: https://www.mathworks.com/matlabcentral/answers/406-how-do-i-skip-items-in-a-legend
      set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
      xlim([min(times), max(times)]);
      hold off
      title(sprintf('%s evoked potentials',channelNames{channel_i}), 'FontSize',18);
      xlabel('time after stimulus (ms)', 'FontSize',18);
      ylabel('lfp (uV)', 'FontSize',18);
      set(gca,'fontsize',18);
      clear figData
      figData.y = responses;
      figData.e = responseErrs;
      figData.x = times;
      drawnow;
      saveFigure(outDir,sprintf('Evoked_byCat_%s_%s_Run%s',channelNames{channel_i},groupName,runNum), figData, figStruct, figStruct.figTag );
    end
  end
end

%evoked analog in potentials
if isfield(plotSwitch,'analogInByItem') && plotSwitch.analogInByItem
  for group_i = 1:length(analysisGroups.analogInPotentials.groups)
    group = analysisGroups.analogInPotentials.groups{group_i};
    groupChannels = analysisGroups.analogInPotentials.channels{group_i};
    groupName = analysisGroups.analogInPotentials.names{group_i};
    groupUnit = analysisGroups.analogInPotentials.units{group_i};
    fh = figure();
    responses = zeros(length(group),length(times));
    responseErrs = zeros(length(group),length(times));
    for item_i = 1:length(group)
      if isfield(catInds,group{item_i})
        responses(item_i,:) = squeeze(mean(sqrt(sum(analogInByCategory{catInds.(group{item_i})}(:,groupChannels,:,lfpPaddedBy+1:end-lfpPaddedBy).^2,2)),3))';
        responseErrs(item_i,:) = squeeze(std(sqrt(sum(analogInByCategory{catInds.(group{item_i})}(:,groupChannels,:,lfpPaddedBy+1:end-lfpPaddedBy).^2,2)),[],3)/sqrt(size(analogInByCategory{catInds.(group{item_i})},3)))';
      else
        for channel_i = 1:length(groupChannels)
        end
        responses(item_i,:) = squeeze(mean(sqrt(sum(analogInByEvent{imInds.(group{item_i})}(:,groupChannels,:,lfpPaddedBy+1:end-lfpPaddedBy).^2,2)),3))';
        responseErrs(item_i,:) = squeeze(std(sqrt(sum(analogInByEvent{imInds.(group{item_i})}(:,groupChannels,:,lfpPaddedBy+1:end-lfpPaddedBy).^2,2)),[],3)/sqrt(size(analogInByEvent{imInds.(group{item_i})},3)))';
      end
    end
    lineProps.width = 3;
    lineProps.col = analysisGroups.analogInPotentials.colors{group_i};
    hold on
    mseb(repmat(times,length(group),1),responses,responseErrs,lineProps);
    legend(group);
    h = get(gca,'ylim');
    h = plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
    %this turns off legend for this line. Source: https://www.mathworks.com/matlabcentral/answers/406-how-do-i-skip-items-in-a-legend
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    xlim([min(times), max(times)]);
    hold off
    title(groupName, 'FontSize',18);
    xlabel('time after stimulus (ms)', 'FontSize',18);
    ylabel(groupUnit, 'FontSize',18);
    set(gca,'fontsize',18);
    clear figData
    figData.y = responses;
    figData.e = responseErrs;
    figData.x = times;
    drawnow;
    saveFigure(outDir,sprintf('%s_Run%s',groupName,runNum), figData, figStruct, figStruct.figTag );
  end
end

%evoked analog in derivatives
if isfield(plotSwitch,'analogInDerivativesByItem') && plotSwitch.analogInDerivativesByItem
  for group_i = 1:length(analysisGroups.analogInPotentials.groups)
    group = analysisGroups.analogInDerivatives.groups{group_i};
    groupChannels = analysisGroups.analogInDerivatives.channels{group_i};
    groupName = analysisGroups.analogInDerivatives.names{group_i};
    groupUnit = analysisGroups.analogInDerivatives.units{group_i};
    fh = figure();
    responses = zeros(length(group),length(times)-1);
    responseErrs = zeros(length(group),length(times)-1);
    for item_i = 1:length(group)
      if isfield(catInds,group{item_i})
        responses(item_i,:) = squeeze(mean(sqrt(sum(diff(analogInByCategory{catInds.(group{item_i})}(:,groupChannels,:,lfpPaddedBy+1:end-lfpPaddedBy),1,4).^2,2)),3))';
        responseErrs(item_i,:) = squeeze(std(sqrt(sum(diff(analogInByCategory{catInds.(group{item_i})}(:,groupChannels,:,lfpPaddedBy+1:end-lfpPaddedBy),1,4).^2,2)),[],3)/sqrt(size(analogInByCategory{catInds.(group{item_i})},3)))';
      else
        for channel_i = 1:length(groupChannels)
        end
        responses(item_i,:) = squeeze(mean(sqrt(sum(diff(analogInByEvent{imInds.(group{item_i})}(:,groupChannels,:,lfpPaddedBy+1:end-lfpPaddedBy),1,4).^2,2)),3))';
        responseErrs(item_i,:) = squeeze(std(sqrt(sum(diff(analogInByEvent{imInds.(group{item_i})}(:,groupChannels,:,lfpPaddedBy+1:end-lfpPaddedBy),1,4).^2,2)),[],3)/sqrt(size(analogInByEvent{imInds.(group{item_i})},3)))';
      end
    end
    lineProps.width = 3;
    lineProps.col = analysisGroups.analogInDerivatives.colors{group_i};
    hold on
    mseb(repmat(times(1:end-1),length(group),1),responses,responseErrs,lineProps);
    legend(group);
    h = get(gca,'ylim');
    h = plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
    %this turns off legend for this line. Source: https://www.mathworks.com/matlabcentral/answers/406-how-do-i-skip-items-in-a-legend
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    xlim([min(times), max(times)]);
    hold off
    title(groupName, 'FontSize',18);
    xlabel('time after stimulus (ms)', 'FontSize',18);
    ylabel(groupUnit, 'FontSize',18);
    set(gca,'fontsize',18);
    clear figData
    figData.y = responses;
    figData.e = responseErrs;
    figData.x = times;
    drawnow;
    saveFigure(outDir,sprintf('%s_Run%s',groupName,runNum), figData, figStruct, figStruct.figTag );
  end
end

% lfp - (color) psth subplot
for channel_i = 1:length(lfpChannels)  
  if isfield(plotSwitch,'colorPsthEvoked') && plotSwitch.colorPsthEvoked
    for group_i = 1:length(analysisGroups.colorPsthEvoked.groups)
      group = analysisGroups.colorPsthEvoked.groups{group_i};
      groupName = analysisGroups.colorPsthEvoked.names{group_i};
      lfpResponses = zeros(length(group),length(times));
      lfpResponseErrs = zeros(length(group),length(times));
      spkResponses = zeros(length(group),length(times));
      for item_i = 1:length(group)
        if isfield(catInds,group{item_i})
          lfpResponses(item_i,:) = squeeze(mean(lfpByCategory{catInds.(group{item_i})}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3))';
          lfpResponseErrs(item_i,:) = squeeze(std(lfpByCategory{catInds.(group{item_i})}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),[],3)/sqrt(size(lfpByCategory{catInds.(group{item_i})},3)))';
          spkResponses(item_i,:) = psthByCategory{channel_i}{end}(catInds.(group{item_i}),:);
        else
          lfpResponses(item_i,:) = squeeze(mean(lfpByEvent{imInds.(group{item_i})}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3))';
          lfpResponseErrs(item_i,:) = squeeze(std(lfpByEvent{imInds.(group{item_i})}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),[],3)/sqrt(size(lfpByEvent{imInds.(group{item_i})},3)))';
          spkResponses(item_i,:) = psthByImage{channel_i}{end}(imInds.(group{item_i}),:);
        end
      end
      fh = figure();
      ah1 = subplot(2,1,1);
      psthTitle = sprintf('%s, %s',channelNames{channel_i}, channelUnitNames{channel_i}{end});
      plotPSTH(spkResponses, [], [], psthParams, 'color', psthTitle, group);
%       yyaxis right
%       plot(times,mean(spkResponses,1),'Color',[0.8,0.8,0.9],'LineWidth',2); %todo: this mean should probably be weighted by number of images per category
%       ylabel('Mean PSTH (Hz)');
      xlim([min(times), max(times)]);
      hold off
      
      ah2 = subplot(2,1,2);
      lineProps.width = 3;
      lineProps.col = analysisGroups.colorPsthEvoked.colors{group_i};
      hold on
      mseb(repmat(times,length(group),1),lfpResponses,lfpResponseErrs,lineProps);
      legend(group,'Location','northeastoutside','FontSize',10);
      h = get(gca,'ylim');
      h = plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
      %this turns off legend for this line. Source: https://www.mathworks.com/matlabcentral/answers/406-how-do-i-skip-items-in-a-legend
      set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
      xlim([min(times), max(times)]);
      xlabel('time after stimulus (ms)', 'FontSize',14);
      ylabel('lfp (uV)', 'FontSize',14);
      set(gca,'fontsize',14);
      % align x axes with colorbar; see http://stackoverflow.com/questions/5259853/matlab-how-to-align-the-axes-of-subplots-when-one-of-them-contains-a-colorbar
      % note: drawnow line synchronizes the rendering thread with the main
      drawnow;
      pos1 = get(ah1,'Position');
      pos2 = get(ah2,'Position');
      pos2(1) = max(pos1(1),pos2(1)); % right limit
      pos1(1) = max(pos1(1),pos2(1));
      pos2(3) = min(pos1(3),pos2(3)); % left limit
      pos1(3) = min(pos1(3),pos2(3));
      set(ah2,'Position',pos2);
      set(ah1,'Position',pos1);
      linkaxes([ah1,ah2],'x');
      title(ah1,sprintf('%s PSTH and Evoked Potentials',channelNames{channel_i}));
      clear figData
      figData.ax11.z = spkResponses;
      figData.ax11.x = times;
      figData.ax21.y = lfpResponses;
      figData.ax21.e = lfpResponseErrs;
      figData.ax21.x = times;
      drawnow;
      saveFigure(outDir,sprintf('PSTH_(color)_Evoked_%s_%s_Run%s',channelNames{channel_i},groupName,runNum), figData, figStruct, figStruct.figTag );
    end
  end
end

% lfp - (line) psth subplot
for channel_i = 1:length(lfpChannels)  
  if length(channelUnitNames{channel_i}) == 2
    numSubplots = 2;
  else
    numSubplots = length(channelUnitNames{channel_i}) + 1;
  end
  if isfield(plotSwitch,'linePsthEvoked') && plotSwitch.linePsthEvoked
    for group_i = 1:length(analysisGroups.linePsthEvoked.groups)
      group = analysisGroups.linePsthEvoked.groups{group_i};
      groupName = analysisGroups.linePsthEvoked.names{group_i};
      fh = figure();
      subplotNum = 0;
      for unit_i = 1:length(channelUnitNames{channel_i})
        if length(channelUnitNames{channel_i}) == 2 && unit_i == 1
          continue
        end
        subplotNum = subplotNum + 1;
        lfpResponses = zeros(length(group),length(times));
        lfpResponseErrs = zeros(length(group),length(times));
        spkResponses = zeros(length(group),length(times));
        spkResponseErrs = zeros(length(group),length(times));
        for item_i = 1:length(group)
          if isfield(catInds,group{item_i})
            spkResponses(item_i,:) = psthByCategory{channel_i}{unit_i}(catInds.(group{item_i}),:);
            spkResponseErrs(item_i,:) = psthErrByCategory{channel_i}{unit_i}(catInds.(group{item_i}),:);
          else
            spkResponses(item_i,:) = psthByImage{channel_i}{unit_i}(imInds.(group{item_i}),:);
            spkResponseErrs(item_i,:) = psthErrByImage{channel_i}{unit_i}(imInds.(group{item_i}),:);
          end
        end
        subplot(numSubplots,1,subplotNum);
        lineProps.width = 3;
        lineProps.col = analysisGroups.linePsthEvoked.colors{group_i};
        hold on
        mseb(repmat(times,length(group),1),spkResponses,cat(3,spkResponseErrs,min(spkResponseErrs,spkResponses)),lineProps);
        xlim([min(times), max(times)]);
        ylabel('Mean PSTH (Hz)');
        title(channelUnitNames{channel_i}{unit_i},'Fontsize',9);
        if subplotNum == 1
          legend(group,'Location','northeastoutside','FontSize',10);
          drawnow;
          pos1 = get(gca,'Position');
        else
          pos2 = get(gca,'Position');
          pos2(1) = pos1(1);
          pos2(3) = pos1(3);
          set(gca,'Position',pos2);
        end
      end
      for item_i = 1:length(group)
        if isfield(catInds,group{item_i})
          lfpResponses(item_i,:) = squeeze(mean(lfpByCategory{catInds.(group{item_i})}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3))';
          lfpResponseErrs(item_i,:) = squeeze(std(lfpByCategory{catInds.(group{item_i})}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),[],3)/sqrt(size(lfpByCategory{catInds.(group{item_i})},3)))';
        else
          lfpResponses(item_i,:) = squeeze(mean(lfpByEvent{imInds.(group{item_i})}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),3))';
          lfpResponseErrs(item_i,:) = squeeze(std(lfpByEvent{imInds.(group{item_i})}(:,channel_i,:,lfpPaddedBy+1:end-lfpPaddedBy),[],3)/sqrt(size(lfpByEvent{imInds.(group{item_i})},3)))';
        end
      end
      subplot(numSubplots,1,numSubplots);
      lineProps.width = 3;
      lineProps.col = analysisGroups.linePsthEvoked.colors{group_i};
      hold on
      mseb(repmat(times,length(group),1),lfpResponses,lfpResponseErrs,lineProps);
      h = get(gca,'ylim');
      h = plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
      %this turns off legend for this line. Source: https://www.mathworks.com/matlabcentral/answers/406-how-do-i-skip-items-in-a-legend
      set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
      xlim([min(times), max(times)]);
      xlabel('time after stimulus (ms)');
      ylabel('lfp (uV)');
      title('LFP','Fontsize',9);
      pos2 = get(gca,'Position');
      pos2(1) = pos1(1);
      pos2(3) = pos1(3);
      set(gca,'Position',pos2);
      suptitle(sprintf('%s PSTH and Evoked Potentials',channelNames{channel_i}));
      clear figData
      figData.ax11.z = spkResponses;
      figData.ax11.x = times;
      figData.ax21.y = lfpResponses;
      figData.ax21.e = lfpResponseErrs;
      figData.ax21.x = times;
      drawnow;
      saveFigure(outDir,sprintf('PSTH_(line)_Evoked_%s_%s_Run%s',channelNames{channel_i},groupName,runNum), figData, figStruct, figStruct.figTag );
    end
  end
end

if isfield(plotSwitch,'colorPsthItemMarginals') && plotSwitch.colorPsthItemMarginals
  for channel_i = 1:length(channelUnitNames) 
    for unit_i = 1:length(channelUnitNames{channel_i})
      for group_i = 1:length(analysisGroups.colorPsthItemMarginals.groups)
        group = analysisGroups.colorPsthItemMarginals.groups{group_i};
        groupName = analysisGroups.colorPsthItemMarginals.names{group_i};
        spkResponses = zeros(length(group),length(times));
        for item_i = 1:length(group)
          if isfield(catInds,group{item_i})
            spkResponses(item_i,:) = psthByCategory{channel_i}{unit_i}(catInds.(group{item_i}),:);
          else
            spkResponses(item_i,:) = psthByImage{channel_i}{unit_i}(imInds.(group{item_i}),:);
          end
        end
        fh = figure();
        psthTitle = '';%suptitle(sprintf('%s PSTH and Evoked Potentials',channelNames{channel_i}));
        plotPsthMarginal(spkResponses, psthPre, psthPost, psthImDur,[frCalcOn,frCalcOff], psthTitle, group);
        drawnow;
        figData = spkResponses;
        saveFigure(outDir,sprintf('PSTH_itemMarginals_%s_%s_%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},groupName,runNum), figData, figStruct);
      end
    end
  end
end

if isfield(plotSwitch,'runSummary') && plotSwitch.runSummary
  for channel_i = 1:length(channelNames)
    fh = figure();
    numSubplots = length(imSpikeCounts{channel_i})-1+4;
    axisHandles = gobjects(numSubplots,1);
    axisNum = 1;
    %juice
    ah = subplot(numSubplots,1,numSubplots);
    ylabel('juice');
    xlabel('time (ms)');
    axisHandles(axisNum) = ah;
    axisNum = axisNum + 1;
    hold on
    assert(length(taskDataAll.juiceOnTimes) == length(taskDataAll.juiceOffTimes),'different number of juice on and off triggers; do not have defense');
    assert(min(taskDataAll.juiceOnTimes) < min(taskDataAll.juiceOffTimes),'juice on at beginning of trial; do not have defense');
    for juice_i = 1:length(taskDataAll.juiceOnTimes)
      plot([taskDataAll.juiceOnTimes(juice_i) taskDataAll.juiceOffTimes(juice_i)],[0 0]);
    end
    % fix spot flash
    ah = subplot(numSubplots,1,numSubplots-1);
    ylabel('flash');
    axisHandles(axisNum) = ah;
    axisNum = axisNum + 1;
    hold on
    assert(min(taskDataAll.fixSpotFlashStartTimes) < min(taskDataAll.fixSpotFlashEndTimes),'fix spot flash off trigger before on; do not have defense');
    %handle case where flash on at end of run with indexing on end times
    for flash_i = 1:length(taskDataAll.fixSpotFlashEndTimes)
      plot([taskDataAll.fixSpotFlashStartTimes(flash_i) taskDataAll.fixSpotFlashEndTimes(flash_i)],[0 0]);
    end
    %fixation
    ah = subplot(numSubplots,1,numSubplots-2);
    ylabel('fix');
    axisHandles(axisNum) = ah;
    axisNum = axisNum + 1;
    hold on
    assert(min(taskDataAll.fixationInTimes) < min(taskDataAll.fixationOutTimes),'fixation out trigger before in; do not have defense');
    %handle case where fix in at end of run
    if length(taskDataAll.fixationInTimes) == length(taskDataAll.fixationOutTimes)+1
      fixationOutTimesTmp = vertcat(taskDataAll.fixationOutTimes,taskDataAll.stimEndTimes(end));
    else
      fixationOutTimesTmp = taskDataAll.fixationOutTimes;
    end
    for fix_i = 1:length(taskDataAll.fixationInTimes)
      plot([taskDataAll.fixationInTimes(fix_i) fixationOutTimesTmp(fix_i)],[0 0]);
    end
    %spikes and lfp
    for image_i = 1:length(onsetsByEvent)
      onsets = onsetsByEvent{image_i};
      lfpByTrial = squeeze(lfpByEvent{image_i}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
      lfpMeanSubByTrial = lfpByTrial - mean(lfpByTrial,2);
      lfpPowerByTrial = sqrt(sum(lfpMeanSubByTrial.^2,2));
      
      for trial_i = 1:length(onsets)
        %trialSubplot = ceil((onsets(trial_i)-firstOnset)/(trialsPerRow*(psthImDur+stimTiming.ISI)));
        %subplot(trialSubplot,1,1);
        
        %lfp
        ah = subplot(numSubplots,1,numSubplots-3);
        if image_i == trial_i && trial_i == 1
          axisHandles(axisNum) = ah;
          axisNum = axisNum + 1;
          ylabel('lfp pwr');
        end
        hold on
        plot([onsets(trial_i) onsets(trial_i)+psthImDur], [lfpPowerByTrial(trial_i) lfpPowerByTrial(trial_i)], 'color','blue');
        %spikes
        for unit_i = 1:length(imSpikeCounts{channel_i})-1 %for now, don't show MUA in addition to units and hash
          ah = subplot(numSubplots,1,unit_i);
          if image_i == trial_i && trial_i == 1
            axisHandles(axisNum) = ah;
            axisNum = axisNum + 1;
            ylabel(sprintf('fr, u%d',unit_i-1));
          end
          hold on
          plot([onsets(trial_i) onsets(trial_i)+psthImDur], [imSpikeCounts{channel_i}{unit_i}{image_i}.counts(trial_i) imSpikeCounts{channel_i}{unit_i}{image_i}.counts(trial_i)], 'color','blue')
        end
      end
    end
    linkaxes(axisHandles,'x');
    suptitle(channelNames{channel_i});
    figData = 'none';
    drawnow;
    saveFigure(outDir,sprintf('runSummary_%s_Run%s',channelNames{channel_i},runNum), figData, figStruct, figStruct.figTag );
  end
end

if isfield(plotSwitch,'runSummaryImMeanSub') && plotSwitch.runSummaryImMeanSub
  for channel_i = 1:length(channelNames)
    fh = figure();
    numSubplots = length(imSpikeCounts{channel_i})-1+4;
    axisHandles = gobjects(numSubplots,1);
    axisNum = 1;
    %juice
    ah = subplot(numSubplots,1,numSubplots);
    ylabel('juice');
    xlabel('time (ms)');
    axisHandles(axisNum) = ah;
    axisNum = axisNum + 1;
    hold on
    assert(length(taskDataAll.juiceOnTimes) == length(taskDataAll.juiceOffTimes),'different number of juice on and off triggers; do not have defense');
    assert(min(taskDataAll.juiceOnTimes) < min(taskDataAll.juiceOffTimes),'juice on at beginning of trial; do not have defense');
    for juice_i = 1:length(taskDataAll.juiceOnTimes)
      plot([taskDataAll.juiceOnTimes(juice_i) taskDataAll.juiceOffTimes(juice_i)],[0 0]);
    end
    % fix spot flash
    ah = subplot(numSubplots,1,numSubplots-1);
    ylabel('flash');
    axisHandles(axisNum) = ah;
    axisNum = axisNum + 1;
    hold on
    assert(min(taskDataAll.fixSpotFlashStartTimes) < min(taskDataAll.fixSpotFlashEndTimes),'fix spot flash off trigger before on; do not have defense');
    %handle case where flash on at end of run with indexing on end times
    for flash_i = 1:length(taskDataAll.fixSpotFlashEndTimes)
      plot([taskDataAll.fixSpotFlashStartTimes(flash_i) taskDataAll.fixSpotFlashEndTimes(flash_i)],[0 0]);
    end
    %fixation
    ah = subplot(numSubplots,1,numSubplots-2);
    ylabel('fix');
    axisHandles(axisNum) = ah;
    axisNum = axisNum + 1;
    hold on
    assert(min(taskDataAll.fixationInTimes) < min(taskDataAll.fixationOutTimes),'fixation out trigger before in; do not have defense');
    %handle case where fix in at end of run
    if length(taskDataAll.fixationInTimes) == length(taskDataAll.fixationOutTimes)+1
      fixationOutTimesTmp = vertcat(taskDataAll.fixationOutTimes,taskDataAll.stimEndTimes(end));
    end
    for fix_i = 1:length(taskDataAll.fixationInTimes)
      plot([taskDataAll.fixationInTimes(fix_i) fixationOutTimesTmp(fix_i)],[0 0]);
    end
    %spikes and lfp
    for image_i = 1:length(onsetsByEvent)
      onsets = onsetsByEvent{image_i};
      lfpByTrial = squeeze(lfpByEvent{image_i}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
      lfpMeanSubByTrial = lfpByTrial - mean(lfpByTrial,2);
      lfpPowerByTrial = sqrt(sum(lfpMeanSubByTrial.^2,2));
      lfpPowerByTrialImMeanSub = lfpPowerByTrial - mean(lfpPowerByTrial);
      
      for trial_i = 1:length(onsets)
        %trialSubplot = ceil((onsets(trial_i)-firstOnset)/(trialsPerRow*(psthImDur+stimTiming.ISI)));
        %subplot(trialSubplot,1,1);
        
        %lfp
        ah = subplot(numSubplots,1,numSubplots-3);
        if image_i == trial_i && trial_i == 1
          axisHandles(axisNum) = ah;
          axisNum = axisNum + 1;
          ylabel('lfp pwr fluct.');
        end
        hold on
        plot([onsets(trial_i) onsets(trial_i)+psthImDur], [lfpPowerByTrialImMeanSub(trial_i) lfpPowerByTrialImMeanSub(trial_i)], 'color','blue');
        %spikes
        for unit_i = 1:length(imSpikeCounts{channel_i})-1 %for now, don't show MUA in addition to units and hash
          ah = subplot(numSubplots,1,unit_i);
          if image_i == trial_i && trial_i == 1
            axisHandles(axisNum) = ah;
            axisNum = axisNum + 1;
            ylabel(sprintf('fr fluct., u%d',unit_i-1));
          end
          hold on
          imageUnitMeanFr = mean(imSpikeCounts{channel_i}{unit_i}{image_i}.counts);
          plot([onsets(trial_i) onsets(trial_i)+psthImDur], [imSpikeCounts{channel_i}{unit_i}{image_i}.counts(trial_i) imSpikeCounts{channel_i}{unit_i}{image_i}.counts(trial_i)] - imageUnitMeanFr, 'color','blue')
        end
      end
    end
    linkaxes(axisHandles,'x');
    suptitle(sprintf('LFP and FR trial fluctuations, %s',channelNames{channel_i}));
    figData = 'none';
    drawnow;
    saveFigure(outDir,sprintf('runSummary_fluct_%s_Run%s',channelNames{channel_i},runNum), figData, figStruct, figStruct.figTag );
  end
end

if isfield(plotSwitch,'runSummaryImMeanSubDiv') && plotSwitch.runSummaryImMeanSubDiv
  for channel_i = 1:length(channelNames)
    fh = figure();
    numSubplots = length(imSpikeCounts{channel_i})-1+4;
    axisHandles = gobjects(numSubplots,1);
    axisNum = 1;
    %juice
    ah = subplot(numSubplots,1,numSubplots);
    ylabel('juice');
    xlabel('time (ms)');
    axisHandles(axisNum) = ah;
    axisNum = axisNum + 1;
    hold on
    assert(length(taskDataAll.juiceOnTimes) == length(taskDataAll.juiceOffTimes),'different number of juice on and off triggers; do not have defense');
    assert(min(taskDataAll.juiceOnTimes) < min(taskDataAll.juiceOffTimes),'juice on at beginning of trial; do not have defense');
    for juice_i = 1:length(taskDataAll.juiceOnTimes)
      plot([taskDataAll.juiceOnTimes(juice_i) taskDataAll.juiceOffTimes(juice_i)],[0 0]);
    end
    % fix spot flash
    ah = subplot(numSubplots,1,numSubplots-1);
    ylabel('flash');
    axisHandles(axisNum) = ah;
    axisNum = axisNum + 1;
    hold on
    assert(min(taskDataAll.fixSpotFlashStartTimes) < min(taskDataAll.fixSpotFlashEndTimes),'fix spot flash off trigger before on; do not have defense');
    %handle case where flash on at end of run with indexing on end times
    for flash_i = 1:length(taskDataAll.fixSpotFlashEndTimes)
      plot([taskDataAll.fixSpotFlashStartTimes(flash_i) taskDataAll.fixSpotFlashEndTimes(flash_i)],[0 0]);
    end
    %fixation
    ah = subplot(numSubplots,1,numSubplots-2);
    ylabel('fix');
    axisHandles(axisNum) = ah;
    axisNum = axisNum + 1;
    hold on
    assert(min(taskDataAll.fixationInTimes) < min(taskDataAll.fixationOutTimes),'fixation out trigger before in; do not have defense');
    %handle case where fix in at end of run
    if length(taskDataAll.fixationInTimes) == length(taskDataAll.fixationOutTimes)+1
      fixationOutTimesTmp = vertcat(taskDataAll.fixationOutTimes,taskDataAll.stimEndTimes(end));
    end
    for fix_i = 1:length(taskDataAll.fixationInTimes)
      plot([taskDataAll.fixationInTimes(fix_i) fixationOutTimesTmp(fix_i)],[0 0]);
    end
    %spikes and lfp
    frCalcOn = frEpochs(1,1);
    frCalcOff = frEpochs(1,2);
    for image_i = 1:length(onsetsByEvent)
      onsets = onsetsByEvent{image_i};
      lfpByTrial = squeeze(lfpByEvent{image_i}(1,channel_i,:,lfpPaddedBy+1+psthPre+frCalcOn:lfpPaddedBy+1+psthPre+frCalcOff));
      lfpMeanSubByTrial = lfpByTrial - mean(lfpByTrial,2);
      lfpPowerByTrial = sqrt(sum(lfpMeanSubByTrial.^2,2));
      lfpPowerByTrialImMeanSubDiv = (lfpPowerByTrial - mean(lfpPowerByTrial))/mean(lfpPowerByTrial);
      
      for trial_i = 1:length(onsets)
        %trialSubplot = ceil((onsets(trial_i)-firstOnset)/(trialsPerRow*(psthImDur+stimTiming.ISI)));
        %subplot(trialSubplot,1,1);
        
        %lfp
        ah = subplot(numSubplots,1,numSubplots-3);
        if image_i == trial_i && trial_i == 1
          axisHandles(axisNum) = ah;
          axisNum = axisNum + 1;
          ylabel('lfp pwr fluct.');
        end
        hold on
        plot([onsets(trial_i) onsets(trial_i)+psthImDur], [lfpPowerByTrialImMeanSubDiv(trial_i) lfpPowerByTrialImMeanSubDiv(trial_i)], 'color','blue');
        %spikes
        for unit_i = 1:length(imSpikeCounts{channel_i})-1 %for now, don't show MUA in addition to units and hash
          ah = subplot(numSubplots,1,unit_i);
          if image_i == trial_i && trial_i == 1
            axisHandles(axisNum) = ah;
            axisNum = axisNum + 1;
            ylabel(sprintf('fr fluct., u%d',unit_i-1));
          end
          hold on
          imageUnitMeanFr = mean(imSpikeCounts{channel_i}{unit_i}{image_i}.counts);
          plot([onsets(trial_i) onsets(trial_i)+psthImDur], ([imSpikeCounts{channel_i}{unit_i}{image_i}.counts(trial_i) imSpikeCounts{channel_i}{unit_i}{image_i}.counts(trial_i)] - imageUnitMeanFr)/imageUnitMeanFr, 'color','blue');
        end
      end
    end
    linkaxes(axisHandles,'x');
    suptitle(sprintf('LFP and FR trial fractional fluctuations, %s',channelNames{channel_i}));
    figData = 'none';
    drawnow;
    saveFigure(outDir,sprintf('runSummary_fracfluct_%s_Run%s',channelNames{channel_i},runNum), figData, figStruct, figStruct.figTag );
  end
end

%%% scatter plots %%%
if isfield(plotSwitch,'lfpPowerMuaScatter') && plotSwitch.lfpPowerMuaScatter
  for epoch_i = 1:size(frEpochs,1)
    for channel_i = 1:length(channelNames)
      for unit_i = 1:length(channelUnitNames{channel_i})
        for group_i = 1:length(analysisGroups.stimulusLabelGroups.groups)
          fh = figure();
          group = analysisGroups.stimulusLabelGroups.groups{group_i};
          groupName = analysisGroups.stimulusLabelGroups.names{group_i};
          groupColors = analysisGroups.stimulusLabelGroups.colors{group_i};
          hold on;
          handlesForLegend = gobjects(length(group),1);
          for item_i = 1:length(group)
            if isfield(catInds,group{item_i})
              lfpByTrial = squeeze(lfpByCategory{catInds.(group{item_i})}(1,channel_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
              spikeCountsByTrial = spikeCountsByCategoryByEpoch{epoch_i}{channel_i}{unit_i}{catInds.(group{item_i})}.counts;
            else
              lfpByTrial = squeeze(lfpByEvent{imInds.(group{item_i})}(1,channel_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
              spikeCountsByTrial = spikeCountsByImageByEpoch{epoch_i}{channel_i}{unit_i}{imInds.(group{item_i})}.rates;
            end
            lfpMeanSubByTrial = lfpByTrial - mean(lfpByTrial,2);
            lfpPowerByTrial = sqrt(sum(lfpMeanSubByTrial.^2,2));
            h = scatter(lfpPowerByTrial,spikeCountsByTrial,36,repmat(groupColors(item_i,:),size(lfpByTrial,1),1),'filled');
            handlesForLegend(item_i) = h;
          end
          xlabel('total lfp power (uV)');
          ylabel('firing rate (Hz)');
          title(sprintf('LFP power vs Firing rate, %s %s, %d ms - %d ms',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}, frEpochs(epoch_i,1), frEpochs(epoch_i,2)))
          legend(handlesForLegend,group,'location','northeastoutside');
          drawnow;
          clear figData
          figData.x = [];
          figData.y = [];
          saveFigure(outDir,sprintf('scatter_lfpPwrVsFr_%s_%s_%s_%dms-%dms_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},groupName,frEpochs(epoch_i,1),frEpochs(epoch_i,2),runNum), figData, figStruct, figStruct.figTag );
        end
      end
    end
  end
end
  
if isfield(plotSwitch,'lfpPeakToPeakMuaScatter') && plotSwitch.lfpPeakToPeakMuaScatter
  for epoch_i = 1:size(frEpochs,1)
    for channel_i = 1:length(channelNames)
      for unit_i = 1:length(channelUnitNames{channel_i})
        for group_i = 1:length(analysisGroups.stimulusLabelGroups.groups)
          group = analysisGroups.stimulusLabelGroups.groups{group_i};
          groupName = analysisGroups.stimulusLabelGroups.names{group_i};
          groupColors = analysisGroups.stimulusLabelGroups.colors{group_i};
          fh = figure();
          hold on;
          handlesForLegend = gobjects(length(group),1);
          for item_i = 1:length(group)
            if isfield(catInds,group{item_i})
              lfpByTrial = squeeze(lfpByCategory{catInds.(group{item_i})}(1,channel_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
              spikeCountsByTrial = spikeCountsByCategoryByEpoch{epoch_i}{channel_i}{unit_i}{catInds.(group{item_i})}.counts;
            else
              lfpByTrial = squeeze(lfpByEvent{imInds.(group{item_i})}(1,channel_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
              spikeCountsByTrial = spikeCountsByImageByEpoch{epoch_i}{channel_i}{unit_i}{imInds.(group{item_i})}.rates;
            end
            h = scatter(max(lfpByTrial,[],2)-min(lfpByTrial,[],2),spikeCountsByTrial,36,groupColors(item_i,:),'filled');
            handlesForLegend(item_i) = h;
          end
          xlabel('Peak-to-Trough LFP Amplitude (uV)');
          ylabel('firing rate (Hz)');
          title(sprintf('LFP Amplitude vs Firing rate, %s %s, %d ms - %d ms',channelNames{channel_i},channelUnitNames{channel_i}{unit_i}, frEpochs(epoch_i,1),frEpochs(epoch_i,2)))
          legend(handlesForLegend,group,'location','northeastoutside');
          drawnow;
          clear figData
          figData.x = [];
          figData.y = [];
          saveFigure(outDir,sprintf('scatter_lfpP2PVsFr_%s_%s_%s_%dms-%dms_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},groupName,frEpochs(epoch_i,1),frEpochs(epoch_i,2),runNum), figData, figStruct, figStruct.figTag );
        end
      end
    end
  end
end 
  
if isfield(plotSwitch,'lfpLatencyMuaLatency') && plotSwitch.lfpLatencyMuaLatency
  for channel_i = 1:length(channelNames)
  end
end
  
if isfield(plotSwitch,'lfpPowerAcrossChannels') && plotSwitch.lfpPowerAcrossChannels && channel_i < length(lfpChannels)
  for epoch_i = 1:size(frEpochs,1)
    for channel_i = 1:length(channelNames)
      for channel2_i = channel_i+1:length(lfpChannels)
        for group_i = 1:length(analysisGroups.stimulusLabelGroups.groups)
          fh = figure();
          group = analysisGroups.stimulusLabelGroups.groups{group_i};
          groupName = analysisGroups.stimulusLabelGroups.names{group_i};
          groupColors = analysisGroups.stimulusLabelGroups.colors{group_i};
          hold on;
          handlesForLegend = gobjects(length(group),1);
          for item_i = 1:length(group)
            if isfield(catInds,group{item_i})
              lfpByTrialCh1 = squeeze(lfpByCategory{catInds.(group{item_i})}(1,channel_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
              lfpByTrialCh2 = squeeze(lfpByCategory{catInds.(group{item_i})}(1,channel2_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
            else
              lfpByTrialCh1 = squeeze(lfpByEvent{imInds.(group{item_i})}(1,channel_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
              lfpByTrialCh2 = squeeze(lfpByEvent{imInds.(group{item_i})}(1,channel2_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
            end
            lfpMeanSubByTrialCh1 = lfpByTrialCh1 - mean(lfpByTrialCh1,2);
            lfpPowerByTrialCh1 = sqrt(sum(lfpMeanSubByTrialCh1.^2,2));
            lfpMeanSubByTrialCh2 = lfpByTrialCh2 - mean(lfpByTrialCh2,2);
            lfpPowerByTrialCh2 = sqrt(sum(lfpMeanSubByTrialCh2.^2,2));
            h = scatter(lfpPowerByTrialCh1,lfpPowerByTrialCh2,36,groupColors(item_i,:),'filled');
            handlesForLegend(item_i) = h;
          end
          xlabel(sprintf('total lfp power, %s (uV)',channelNames{channel_i}));
          ylabel(sprintf('total lfp power, %s (uV)',channelNames{channel2_i}));
          title(sprintf('LFP Power, %s vs. %s, %d ms - %d ms',channelNames{channel_i},channelNames{channel2_i}, frEpochs(epoch_i,1), frEpochs(epoch_i,2)))
          legend(handlesForLegend,group,'location','northeastoutside');
          drawnow;
          clear figData
          figData.x = [];
          figData.y = [];
          saveFigure(outDir,sprintf('scatter_lfpPwrVSlfpPwr_%s_%s_%s_%dms-%dms_Run%s',channelNames{channel_i},channelNames{channel2_i},groupName,frEpochs(epoch_i,1),frEpochs(epoch_i,2),runNum), figData, figStruct, figStruct.figTag );
        end
      end
    end
  end
end
  
% lfp p2p amplitude vs lfp p2p amplitude, across channels
if isfield(plotSwitch,'lfpPeakToPeakAcrossChannels') && plotSwitch.lfpPeakToPeakAcrossChannels && channel_i < length(lfpChannels)
  for epoch_i = 1:size(frEpochs,1)
    for channel_i = 1:length(channelNames)
      for channel2_i = channel_i+1:length(lfpChannels)
        for group_i = 1:length(analysisGroups.stimulusLabelGroups.groups)
          fh = figure();
          group = analysisGroups.stimulusLabelGroups.groups{group_i};
          groupName = analysisGroups.stimulusLabelGroups.names{group_i};
          groupColors = analysisGroups.stimulusLabelGroups.colors{group_i};
          hold on;
          handlesForLegend = gobjects(length(group),1);
          for item_i = 1:length(group)
            if isfield(catInds,group{item_i})
              lfpByTrialCh1 = squeeze(lfpByCategory{catInds.(group{item_i})}(1,channel_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
              lfpByTrialCh2 = squeeze(lfpByCategory{catInds.(group{item_i})}(1,channel2_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
            else
              lfpByTrialCh1 = squeeze(lfpByEvent{imInds.(group{item_i})}(1,channel_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
              lfpByTrialCh2 = squeeze(lfpByEvent{imInds.(group{item_i})}(1,channel2_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
            end
            h = scatter(max(lfpByTrialCh1,[],2)-min(lfpByTrialCh1,[],2),max(lfpByTrialCh2,[],2)-min(lfpByTrialCh2,[],2),36,groupColors(item_i,:),'filled');
            handlesForLegend(item_i) = h;
          end
          xlabel(sprintf('LFP Amplitude, %s (uV)',channelNames{channel_i}));
          ylabel(sprintf('LFP Amplitude, %s (uV)',channelNames{channel2_i}));
          title(sprintf('LFP Amplitude, %s vs. %s, %d ms - %d ms',channelNames{channel_i},channelNames{channel2_i}, frEpochs(epoch_i,1), frEpochs(epoch_i,2)))
          legend(handlesForLegend,group,'location','northeastoutside');
          drawnow;
          clear figData
          figData.x = [];
          figData.y = [];
          saveFigure(outDir,sprintf('scatter_lfpP2PVSlfpP2P_%s_%s_%s_%dms-%dms_Run%s',channelNames{channel_i},channelNames{channel2_i},groupName,frEpochs(epoch_i,1),frEpochs(epoch_i,2),runNum), figData, figStruct, figStruct.figTag );
        end
      end
    end
  end
end
  
% lfp latency vs lfp latency, across channels, defined as max dot product with mean
if isfield(plotSwitch,'lfpLatencyShiftAcrossChannels') && plotSwitch.lfpLatencyShiftAcrossChannels && channel_i < length(lfpChannels)
  for channel_i = 1:length(channelNames)
    for channel2_i = channel_i+1:length(lfpChannels)
      for group_i = 1:length(analysisGroups.stimulusLabelGroups.groups)
        fh = figure();
        group = analysisGroups.stimulusLabelGroups.groups{group_i};
        groupName = analysisGroups.stimulusLabelGroups.names{group_i};
        groupColors = analysisGroups.stimulusLabelGroups.colors{group_i};
        hold on;
        handlesForLegend = gobjects(length(group),1);
        for item_i = 1:length(group)
          if isfield(catInds,group{item_i})
            lfpByTrialCh1 = squeeze(lfpByCategory{catInds.(group{item_i})}(1,channel_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
            lfpByTrialCh2 = squeeze(lfpByCategory{catInds.(group{item_i})}(1,channel2_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
          else
            lfpByTrialCh1 = squeeze(lfpByEvent{imInds.(group{item_i})}(1,channel_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
            lfpByTrialCh2 = squeeze(lfpByEvent{imInds.(group{item_i})}(1,channel2_i,:,lfpPaddedBy+1+psthPre+frEpochs(epoch_i,1):lfpPaddedBy+1+psthPre+frEpochs(epoch_i,2)));
          end
          lfpMeanSubByTrialCh1 = lfpByTrialCh1 - mean(lfpByTrialCh1,2);
          lfpTrialAveCh1 = mean(lfpMeanSubByTrialCh1,1);
          lfpMeanSubByTrialCh2 = lfpByTrialCh2 - mean(lfpByTrialCh2,2);
          lfpTrialAveCh2 = mean(lfpMeanSubByTrialCh2,1);
          shiftsCh1 = zeros(size(lfpByTrialCh1,1),1);
          shiftsCh2 = zeros(size(lfpByTrialCh1,1),1);
          for trial_i = 1:size(lfpByTrialCh1,1)
            bestShiftCh1 = -50;
            bestMatchCh1 = 0;
            for shift = -50:50 %magic number; todo replace with variable
              shiftedTrialLfp = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel_i,trial_i,lfpPaddedBy+1+psthPre+frCalcOn+shift:lfpPaddedBy+1+psthPre+frCalcOff+shift));
              shiftedTrialLfp = shiftedTrialLfp - mean(shiftedTrialLfp);
              match = lfpTrialAveCh1*shiftedTrialLfp;
              if match > bestMatchCh1
                bestMatchCh1 = match;
                bestShiftCh1 = shift;
              end
              shiftsCh1(trial_i) = bestShiftCh1;
            end
            bestShiftCh2 = -50;
            bestMatchCh2 = 0;
            for shift = -50:50 %magic number; todo replace with variable
              shiftedTrialLfp = squeeze(lfpByCategory{categorySlimInds(cat_i)}(1,channel2_i,trial_i,lfpPaddedBy+1+psthPre+frCalcOn+shift:lfpPaddedBy+1+psthPre+frCalcOff+shift));
              shiftedTrialLfp = shiftedTrialLfp - mean(shiftedTrialLfp);
              match = lfpTrialAveCh2*shiftedTrialLfp;
              if match > bestMatchCh2
                bestMatchCh2 = match;
                bestShiftCh2 = shift;
              end
              shiftsCh2(trial_i) = bestShiftCh2;
            end
          end
          h = scatter(shiftsCh1,shiftsCh2,36,groupColors(item_i,:),'filled');
          handlesForLegend(item_i) = h;
        end
        xlabel(sprintf('LFP latency shift from category mean, %s (ms)',channelNames{channel_i}));
        ylabel(sprintf('LFP latency shift from category mean, %s (ms)',channelNames{channel2_i}));
        title(sprintf('LFP Latency shift from category mean, %s vs. %s, %d ms - %d ms',channelNames{channel_i},channelNames{channel2_i}, frCalcOn, frCalcOff))
        legend(handlesForLegend,group,'location','northeastoutside');
        drawnow;
        clear figData
        figData.x = [];
        figData.y = [];
        saveFigure(outDir,sprintf('scatter_lfpShiftVSlfpShift_%s_%s__%s_Run%s',channelNames{channel_i},channelNames{channel2_i},groupName,runNum), figData, figStruct, figStruct.figTag );
      end
    end
  end
end
% todo
if isfield(plotSwitch,'singleTrialAnalogInByCategory') && plotSwitch.singleTrialAnalogInByCategory
  for group_i = 1:length(analysisGroups.analogInSingleTrialsByCategory.groups)
    group = analysisGroups.analogInSingleTrialsByCategory.groups{group_i};
    groupName = analysisGroups.analogInSingleTrialsByCategory.names{group_i};
    groupChannels = analysisGroups.analogInSingleTrialsByCategory.channels{group_i};
    for channel_i = 1:length(groupChannels)
      channel = groupChannels(channel_i);
      fh = figure();
      figData = cell(length(group),1);
      for item_i = 1:length(group)
        subplot(length(group),1,item_i)
        hold on
        ydata = zeros(50,length(times));
        if isfield(catInds,group{item_i})
          for i = 1:min(50,size(analogInByCategory{catInds.(group{item_i})},3))
            plot(times, squeeze(analogInByCategory{catInds.(group{item_i})}(:,channel,i,lfpPaddedBy+1:end-lfpPaddedBy)));
            ydata(i,:) = squeeze(analogInByCategory{catInds.(group{item_i})}(:,channel,i,lfpPaddedBy+1:end-lfpPaddedBy));
          end
        else
          for i = 1:min(50,size(analogInByEvent{imInds.(group{item_i})},3))
            plot(times, squeeze(analogInByEvent{catInds.(group{item_i})}(:,channel,i,lfpPaddedBy+1:end-lfpPaddedBy)));
            ydata(i,:) = squeeze(analogInByEvent{catInds.(group{item_i})}(:,channel,i,lfpPaddedBy+1:end-lfpPaddedBy));
          end
        end
        h = get(gca,'ylim');
        plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
        hold off
        title(sprintf('single trials, %s, %s', group{item_i},analogInChannelNames{channel}), 'FontSize',18);
        xlabel('time after stimulus (ms)', 'FontSize',18);
        ylabel(analogInChannelUnits{channel}, 'FontSize',18);
        set(gca,'fontsize',18);
        xlim([min(times) max(times)]);
        clear tmp
        tmp.y = ydata;
        tmp.x = times;
        figData{item_i} = tmp;
      end
      drawnow;
      saveFigure(outDir,sprintf('%s_singleTrials_%s_Run%s',analogInChannelNames{channel},groupName,runNum), figData, figStruct, figStruct.figTag );
    end
  end
end

% other to-do: one epoch for each peak of trial ave evoked potential, ideally set automatically
% todo: across-peak latency within and across channels
% todo: spike burstiness vs. lfp bumpiness; can do as p2p evoked/psth; or (low freq power)/(high freq power)
calcSwitchNames = {'evoked', 'induced'};
tfCalcSwitchTitleSuffixes = {'',', induced'}; % appended to titles
tfCalcSwitchFnameSuffixes = {'','_induced'}; % appended to filenames
calcSwitches = [calcSwitch.evokedSpectra, calcSwitch.inducedSpectra];
for calc_i = 1:length(calcSwitches)
  if calcSwitches(calc_i)
    if strcmp(calcSwitchNames{calc_i},'induced')
      if calcSwitch.spikeTimes
        disp('induced TF not yet implemented for spike times; change to spike bins using calcSwitch.spikeTimes = 0');
        continue
      end
      disp('switching spikes and lfps to induced for spectrum computation');
      spikesByCategoryBinnedEvokedTmp = spikesByCategoryBinned;
      spikesByCategoryBinned = spikesByCategoryBinnedInduced;
      lfpByCategoryEvokedTmp = lfpByCategory;
      lfpByCategory = lfpByCategoryInduced;
    end
    % single trial evoked lfps; better as multi-channel subplot?
    % todo: programatize on category
    if isfield(plotSwitch,'singleTrialLfpByCategory') && plotSwitch.singleTrialLfpByCategory
      for group_i = 1:length(analysisGroups.lfpSingleTrialsByCategory.groups)
        group = analysisGroups.lfpSingleTrialsByCategory.groups{group_i};
        groupName = analysisGroups.lfpSingleTrialsByCategory.names{group_i};
        for channel_i = 1:length(channelNames)
          fh = figure();
          figData = cell(length(group),1);
          for item_i = 1:length(group)
            subplot(length(group),1,item_i)
            hold on
            ydata = zeros(50,length(times));
            if isfield(catInds,group{item_i})
              for i = 1:min(50,size(lfpByCategory{catInds.(group{item_i})},3))
                plot(times, squeeze(lfpByCategory{catInds.(group{item_i})}(1,channel_i,i,lfpPaddedBy+1:end-lfpPaddedBy)));
                ydata(i,:) = squeeze(lfpByCategory{catInds.(group{item_i})}(1,channel_i,i,lfpPaddedBy+1:end-lfpPaddedBy));
              end
            else
              for i = 1:min(50,size(lfpByEvent{imInds.(group{item_i})},3))
                plot(times, squeeze(lfpByEvent{imInds.(group{item_i})}(1,channel_i,i,lfpPaddedBy+1:end-lfpPaddedBy)));
                ydata(i,:) = squeeze(lfpByEvent{imInds.(group{item_i})}(1,channel_i,i,lfpPaddedBy+1:end-lfpPaddedBy));
              end
            end
            h = get(gca,'ylim');
            plot([0, psthImDur],[h(1)+0.05*(h(2)-h(1)), h(1)+0.05*(h(2)-h(1))],'color','k','linewidth',3);
            hold off
            title(sprintf('single trial LFPs, %s, %s%s', group{item_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
            xlabel('time after stimulus (ms)', 'FontSize',18);
            ylabel('voltage (uV)', 'FontSize',18);
            set(gca,'fontsize',18);
            xlim([min(times) max(times)]);
            clear tmp
            tmp.y = ydata;
            tmp.x = times;
            figData{item_i} = tmp;
          end
          drawnow;
          saveFigure(outDir,sprintf('Evoked_singleTrials_%s_%s%s_Run%s',channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
        end
      end
    end
    
    
    if isfield(plotSwitch,'lfpSpectraByCategory') && plotSwitch.lfpSpectraByCategory
      for group_i = 1:length(analysisGroups.spectraByCategory.groups)
        group = analysisGroups.spectraByCategory.groups{group_i};
        groupName = analysisGroups.spectraByCategory.names{group_i};
        groupColors = analysisGroups.spectraByCategory.colors{group_i};
        for channel_i = 1:length(channelNames)
          for item_i = 1:length(group)
            [S,F,E] = mtspectrumc(squeeze(lfpByCategory{catInds.(group{item_i})}(1,channel_i,:,:))', chronuxParams);
            if item_i == 1
              spectra = zeros(length(group),length(S));
              freqs = zeros(length(group),length(F));
              errs = zeros(length(group),length(E),2);
            end
            spectra(item_i,:) = S';
            freqs(item_i,:) = 1000*F';
            errs(item_i,:,1) = E(1,:)'-S;
            errs(item_i,:,2) = S-E(2,:)';
          end
          fh = figure();
          lineProps.width = 3;
          lineProps.col = analysisGroups.spectraByCategory.colors{group_i};
          mseb(freqs,spectra,errs,lineProps);
          hold on
          % need to plot individually as well, to get point markers
          handles = gobjects(length(group),1);
          for item_i = 1:length(group)
            h = plot(freqs(item_i,:),spectra(item_i,:),'linewidth',3,'linestyle','none','marker','o','color',groupColors{item_i});
            handles(item_i) = h;
          end
          set(gca,'yScale','log');
          title(sprintf('%s evoked power spectra, %s%s',channelNames{channel_i},groupName,tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
          xlabel('frequency (Hz)', 'FontSize',18);
          ylabel('voltage, log(uV)', 'FontSize',18);
          legend(handles,group);
          set(gca,'fontsize',18);
          clear figData
          figData.y = spectra;
          figData.x = freqs;
          figData.e = errs;
          drawnow;
          saveFigure(outDir,sprintf('spectrum_%s_%s_LFP%s_Run%s',groupName,channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
          
          %loglog spectrum with power law fit
          specModel = fitlm(log(1000*mean(freqs(:,6:end),1)), log(mean(spectra(:,6:end),1)));
          m = specModel.Coefficients.Estimate(2);
          y0 = specModel.Coefficients.Estimate(1);
          fh = figure();
          lineProps.width = 3;
          lineProps.col = analysisGroups.spectraByCategory.colors{group_i};
          mseb(freqs,spectra,errs,lineProps);
          hold on
          % need to plot individually as well, to get point markers
          handles = gobjects(length(group),1);
          for item_i = 1:length(group)
            h = plot(freqs(item_i,:),spectra(item_i,:),'linewidth',3,'linestyle','none','marker','o','color',groupColors{item_i});
            handles(item_i) = h;
          end
          h = plot(mean(freqs(:,6:end),1),exp(m*log(1000*mean(freqs(:,6:end),1)) + y0), 'k--');
          set(gca,'yScale','log');
          set(gca,'xScale','log');
          xlabel('frequency (Hz)', 'FontSize',18);
          ylabel('voltage (uV)', 'FontSize',18);
          title(sprintf('%s evoked power spectra, %s%s',channelNames{channel_i},groupName,tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
          legend(vertcat(handles,h),vertcat(group,{sprintf('fit,m = %.3f',m)}));
          set(gca,'fontsize',18);
          clear figData
          figData.y = spectra;
          figData.x = freqs;
          figData.e = errs;
          drawnow;
          saveFigure(outDir,sprintf('spectrum_loglog_%s_%s_LFP%s_Run%s',groupName,channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
        
          if calcSwitch.trialMeanSpectra
            for item_i = 1:length(group)
              [S,F,E] = mtspectrumc(squeeze(mean(lfpByCategory{catInds.(group{item_i})}(1,channel_i,:,:),3)), chronuxParams);
              if item_i == 1
                spectra = zeros(length(group),length(S));
                freqs = zeros(length(group),length(F));
                errs = zeros(length(group),length(E),2);
              end
              spectra(item_i,:) = S'; %note: explicit transposition not technically necessary, since matlab will transpose automatically for singleton dimensions
              freqs(item_i,:) = 1000*F';
              errs(item_i,:,1) = E(1,:)'-S;
              errs(item_i,:,2) = S-E(2,:)';
            end
            fh = figure();
            lineProps.width = 3;
            lineProps.col = analysisGroups.spectraByCategory.colors{group_i};
            mseb(freqs,spectra,errs,lineProps);
            hold on
            % need to plot individually as well, to get point markers
            handles = gobjects(length(group),1);
            for item_i = 1:length(group)
              h = plot(freqs(item_i,:),spectra(item_i,:),'linewidth',2,'linestyle','none','marker','o','color',groupColors{item_i});
              handles(item_i) = h;
            end
            set(gca,'yScale','log');
            title(sprintf('%s mean-evoked power spectra, %s%s',channelNames{channel_i},groupName,tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
            xlabel('frequency (Hz)', 'FontSize',18);
            ylabel('voltage, log(uV)', 'FontSize',18);
            legend(handles,group);
            set(gca,'fontsize',18);
            clear figData
            figData.y = spectra;
            figData.x = freqs;
            figData.e = errs;
            drawnow;
            saveFigure(outDir,sprintf('spectrum_mean_%s_%s_LFP%s_Run%s',groupName,channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
            
            fh = figure();
            lineProps.width = 3;
            lineProps.col = analysisGroups.spectraByCategory.colors{group_i};
            mseb(freqs,spectra,errs,lineProps);
            hold on
            % need to plot individually as well, to get point markers
            handles = gobjects(length(group),1);
            for item_i = 1:length(group)
              h = plot(freqs(item_i,:),spectra(item_i,:),'linewidth',3,'linestyle','none','marker','o','color',groupColors{item_i});
              handles(item_i) = h;
            end
            set(gca,'yScale','log');
            set(gca,'xScale','log');
            title(sprintf('%s mean-evoked power spectra, %s%s',channelNames{channel_i},groupName,tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
            specModel = fitlm(log(mean(freqs(:,6:end),1)), log(mean(spectra(:,6:end),1)));
            m = specModel.Coefficients.Estimate(2);
            y0 = specModel.Coefficients.Estimate(1);
            h = plot(mean(freqs(:,6:end),1),exp(m*log(mean(freqs(:,6:end),1)) + y0), 'k--');
            xlabel('frequency (Hz)', 'FontSize',18);
            ylabel('voltage (uV)', 'FontSize',18);
            legend(vertcat(handles,h),vertcat(group,{sprintf('fit,m = %.3f',m)}));
            set(gca,'fontsize',18);
            clear figData
            figData.y = spectra;
            figData.x = freqs;
            figData.e = errs;
            drawnow;
            saveFigure(outDir,sprintf('spectrum__mean_loglog_%s_%s_LFP%s_Run%s',groupName,channelNames{channel_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
          end
        end
      end
    end
    % time domain autocorrelation
    if (isfield(plotSwitch,'lfpAutocorrelTfByItem') && plotSwitch.lfpAutocorrelTfByItem) || (isfield(plotSwitch,'lfpAutocorrelByItem') && plotSwitch.lfpAutocorrelByItem)
      for group_i = 1:length(analysisGroups.spectraByCategory.groups)
        group = analysisGroups.spectraByCategory.groups{group_i};
        groupName = analysisGroups.spectraByCategory.names{group_i};
        groupColors = analysisGroups.spectraByCategory.colors{group_i};
        for channel_i = 1:length(channelNames)
          for item_i = 1:length(group)
            if isfield(catInds,group{item_i})
              lfpByItem = lfpByCategory;
              itemNum = catInds.(group{item_i});
            else
              lfpByItem = lfpByEvent;
              itemNum = imInds.(group{item_i});
            end
            [Cgram, C, shiftList, confCgram, confC] = correlogram(squeeze(lfpByItem{itemNum}(1,channel_i,:,:)),...
              squeeze(lfpByItem{itemNum}(1,channel_i,:,:)), correlParams);
            Cgram = Cgram(lfpPaddedBy+1:end-lfpPaddedBy,lfpPaddedBy+1:end-lfpPaddedBy);
            confCgram = confCgram(lfpPaddedBy+1:end-lfpPaddedBy,lfpPaddedBy+1:end-lfpPaddedBy);
            if item_i == 1
              spectra = zeros(length(group),length(C));
              specErrs = zeros(length(group),length(C),2); % support for asymmetric error bars, as from jacknife
              shifts =  zeros(length(group),length(C));
            end
            spectra(item_i,:) = C';
            specErrs(item_i,:,:) = confC;
            shifts(item_i,:) = shiftList';
            t = -psthPre:psthImDur+psthPost;
            if isfield(plotSwitch,'lfpAutocorrelTfByItem') && plotSwitch.lfpAutocorrelTfByItem
              fh = figure();
              imagesc(t,t,Cgram); axis xy
              xlabel('Time (ms)');
              ylabel('Time(ms)');
              c = colorbar();
              ylabel(c,'Correlation');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              line(xlim,[0 0],'Color',[0.8,0.8,0.9],'LineWidth',4);
              line(xlim,[psthImDur psthImDur],'Color',[0.8,0.8,0.9],'LineWidth',4);
              title(sprintf('%s LFP autocorrelation, %s%s',channelNames{channel_i},group{item_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = t;
              figData.y = 1;
              figData.z = Cgram;
              drawnow;
              saveFigure(outDir,sprintf('autocorrel_TF_%s_LFP_%s%s_Run%s',channelNames{channel_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
            end
          end
          if isfield(plotSwitch,'lfpAutocorrelByItem') && plotSwitch.lfpAutocorrelByItem
            fh = figure();
            lineProps.width = 3;
            lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
            mseb(shifts,spectra,specErrs,lineProps);
            legend(group);
            xlabel('time (ms)');
            ylabel('correlation');
            title(sprintf('%s LFP autocorrelation %s',channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = shifts;
            figData.y = spectra;
            figData.e = specErrs;
            drawnow;
            saveFigure(outDir,sprintf('autocorrel_%s_LFP_%s%s_Run%s',channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
          end
        end
      end
    end
    
    
    if isfield(plotSwitch,'spikeSpectraByCategory') && plotSwitch.spikeSpectraByCategory
      for group_i = 1:length(analysisGroups.spectraByCategory.groups)
        group = analysisGroups.spectraByCategory.groups{group_i};
        groupName = analysisGroups.spectraByCategory.names{group_i};
        groupColors = analysisGroups.spectraByCategory.colors{group_i};
        for channel_i = 1:length(channelNames)
          for unit_i = 1:length(channelUnitNames{channel_i})
            if length(channelUnitNames{channel_i}) == 2 && unit_i == 1 %skip unsorted if no isolated units defined
              continue
            end
            for item_i = 1:length(group)
              if calcSwitch.spikeTimes
                [S,F,~,E] = mtspectrumpt(spikesByCategoryForTF{catInds.(group{item_i})}{channel_i}{unit_i}, chronuxParams); %matlab note: '~' is placeholder for unneeded output
              else
                [S,F,~,E] = mtspectrumpb(spikesByCategoryBinned{catInds.(group{item_i})}{channel_i}{unit_i}', chronuxParams);
              end
              if item_i == 1
                spectra = zeros(length(group),length(S));
                freqs = zeros(length(group),length(F));
                errs = zeros(length(group),length(E),2);
              end
              spectra(item_i,:) = S';
              freqs(item_i,:) = 1000*F';
              errs(item_i,:,1) = E(1,:)'-S;
              errs(item_i,:,2) = S-E(2,:)';
            end
            fh = figure();
            lineProps.width = 3;
            lineProps.col = analysisGroups.spectraByCategory.colors{group_i};
            mseb(freqs,spectra,errs,lineProps);
            hold on
            % plot individually as well, to get point markers
            handles = gobjects(length(group),1);
            for item_i = 1:length(group)
              h = plot(freqs(item_i,:),spectra(item_i,:),'linewidth',3,'linestyle','none','marker','o','color',groupColors{item_i});
              handles(item_i) = h;
            end
            set(gca,'yScale','log');
            title(sprintf('%s %s evoked power spectra, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},groupName,tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
            xlabel('frequency (Hz)', 'FontSize',18);
            ylabel('power (spks/s)', 'FontSize',18);
            legend(handles,group);
            set(gca,'fontsize',18);
            clear figData
            figData.y = spectra;
            figData.x = freqs;
            figData.e = errs;
            drawnow;
            saveFigure(outDir,sprintf('spectrum_%s_%s_%s%s_Run%s',groupName,channelNames{channel_i},channelUnitNames{channel_i}{unit_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
            
            fh = figure();
            lineProps.width = 3;
            lineProps.col = analysisGroups.spectraByCategory.colors{group_i};
            mseb(freqs,spectra,errs,lineProps);
            hold on
            % plot individually as well, to get point markers
            handles = gobjects(length(group),1);
            for item_i = 1:length(group)
              h = plot(freqs(item_i,:),spectra(item_i,:),'linewidth',3,'linestyle','none','marker','o','color',groupColors{item_i});
              handles(item_i) = h;
            end
            set(gca,'yScale','log');
            set(gca,'xScale','log');
            specModel = fitlm(log(mean(freqs(:,6:end),1)), log(mean(spectra(:,6:end),1)));
            m = specModel.Coefficients.Estimate(2);
            y0 = specModel.Coefficients.Estimate(1);
            h = plot(mean(freqs(:,6:end),1),exp(m*log(mean(freqs(:,6:end),1)) + y0), 'k--');
            xlabel('frequency (Hz)', 'FontSize',18);
            ylabel('power (spks/s)', 'FontSize',18);
            title(sprintf('%s %s evoked power spectra, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},groupName,tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
            legend(vertcat(handles,h),vertcat(group,{sprintf('fit,m = %.3f',m)}));
            set(gca,'fontsize',18);
            clear figData
            figData.y = spectra;
            figData.x = freqs;
            figData.e = errs;
            saveFigure(outDir,sprintf('spectrum_loglog_%s_%s_%s%s_Run%s',groupName,channelNames{channel_i},channelUnitNames{channel_i}{unit_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
            
            if calcSwitch.trialMeanSpectra
              for item_i = 1:length(group)
                if calcSwitch.spikeTimes
                  [S,F,~,E] = mtspectrumpt(allSpikesByCategoryForTF{catInds.(group{item_i})}{channel_i}{unit_i}, chronuxParams); %matlab note: '~' is placeholder for unneeded output
                else
                  [S,F,~,E] = mtspectrumpb(mean(spikesByCategoryBinned{catInds.(group{item_i})}{channel_i}{unit_i},1)', chronuxParams);
                end
                if item_i == 1
                  spectra = zeros(length(group),length(S));
                  freqs = zeros(length(group),length(F));
                  errs = zeros(length(group),length(E),2);
                end
                spectra(item_i,:) = S';
                freqs(item_i,:) = 1000*F';
                errs(item_i,:,1) = E(1,:)'-S;
                errs(item_i,:,2) = S-E(2,:)';
              end
              fh = figure();
              lineProps.width = 3;
              lineProps.col = analysisGroups.spectraByCategory.colors{group_i};
              mseb(freqs,spectra,errs,lineProps);
              hold on
              % need to plot individually as well, to get point markers
              handles = gobjects(length(group),1);
              for item_i = 1:length(group)
                h = plot(freqs(item_i,:),spectra(item_i,:),'linewidth',3,'linestyle','none','marker','o','color',groupColors{item_i});
                handles(item_i) = h;
              end
              set(gca,'yScale','log');
              title(sprintf('%s %s mean-evoked power spectra, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},groupName,tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
              xlabel('frequency (Hz)', 'FontSize',18);
              ylabel('power (spks/s)', 'FontSize',18);
              legend(handles,group);
              set(gca,'fontsize',18);
              clear figData
              figData.y = spectra;
              figData.x = freqs;
              figData.e = errs;
              drawnow;
              saveFigure(outDir,sprintf('spectrum_mean_%s_%s_%s%s_Run%s',groupName,channelNames{channel_i},channelUnitNames{channel_i}{unit_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              
              fh = figure();
              lineProps.width = 3;
              lineProps.col = analysisGroups.spectraByCategory.colors{group_i};
              mseb(freqs,spectra,errs,lineProps);
              hold on
              % need to plot individually as well, to get point markers
              handles = gobjects(length(group),1);
              for item_i = 1:length(group)
                h = plot(freqs(item_i,:),spectra(item_i,:),'linewidth',3,'linestyle','none','marker','o','color',groupColors{item_i});
                handles(item_i) = h;
              end
              set(gca,'yScale','log');
              set(gca,'xScale','log');
              specModel = fitlm(log(mean(freqs(:,6:end),1)), log(mean(spectra(:,6:end),1)));
              m = specModel.Coefficients.Estimate(2);
              y0 = specModel.Coefficients.Estimate(1);
              h = plot(mean(freqs(:,6:end),1),exp(m*log(mean(freqs(:,6:end),1)) + y0), 'k--');
              xlabel('frequency (Hz)', 'FontSize',18);
              ylabel('power (spks/s)', 'FontSize',18);
              title(sprintf('%s %s mean-evoked power spectra, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},groupName,tfCalcSwitchTitleSuffixes{calc_i}), 'FontSize',18);
              legend(vertcat(handles,h),vertcat(group,{sprintf('fit,m = %.3f',m)}));
              set(gca,'fontsize',18);
              clear figData
              figData.y = spectra;
              figData.x = freqs;
              figData.e = errs;
              saveFigure(outDir,sprintf('spectrum_mean_loglog_%s_%s_%s%s_Run%s',groupName,channelNames{channel_i},channelUnitNames{channel_i}{unit_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
            end
          end
        end
      end
      
      % time domain autocorrelation
      if (isfield(plotSwitch,'spikeAutocorrelTfByItem') && plotSwitch.spikeAutocorrelTfByItem) || (isfield(plotSwitch,'spikeAutocorrelByItem') && plotSwitch.spikeAutocorrelByItem)
        for group_i = 1:length(analysisGroups.spectraByCategory.groups)
          group = analysisGroups.spectraByCategory.groups{group_i};
          groupName = analysisGroups.spectraByCategory.names{group_i};
          groupColors = analysisGroups.spectraByCategory.colors{group_i};
          for channel_i = 1:length(channelNames)
            if calcSwitch.spikeTimes
              disp('Requested time-domain spike autocorrelation analysis for non-binned spikes. Not implemented. Skipping');
            else
              for item_i = 1:length(group)
                if isfield(catInds,group{item_i})
                  spikesByItemBinned = spikesByCategoryBinned;
                  itemNum = catInds.(group{item_i});
                else
                  spikesByItemBinned = spikesByEventBinned;
                  itemNum = imInds.(group{item_i});
                end
                correlParams.useSmoothing = [1 1];
                [Cgram, C, shiftList, confCgram, confC] = correlogram(spikesByItemBinned{itemNum}{channel_i}{unit_i},...
                  spikesByItemBinned{itemNum}{channel_i}{unit_i}, correlParams);
                correlParams.useSmoothing = [0 0];
                Cgram = Cgram(lfpPaddedBy+1:end-lfpPaddedBy,lfpPaddedBy+1:end-lfpPaddedBy);
                confCgram = confCgram(lfpPaddedBy+1:end-lfpPaddedBy,lfpPaddedBy+1:end-lfpPaddedBy);
                if item_i == 1
                  spectra = zeros(length(group),length(C));
                  specErrs = zeros(length(group),length(C),2); % support for asymmetric error bars, as from jacknife
                  shifts =  zeros(length(group),length(C));
                end
                spectra(item_i,:) = C';
                specErrs(item_i,:,:) = confC;
                shifts(item_i,:) = shiftList';
                t = -psthPre:psthImDur+psthPost;
                if isfield(plotSwitch,'spikeAutocorrelTfByItem') && plotSwitch.spikeAutocorrelTfByItem
                  fh = figure();
                  imagesc(t,t,Cgram); axis xy
                  xlabel('Time (ms)');
                  ylabel('Time (ms)');
                  c = colorbar();
                  ylabel(c,'Correlation');
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  line(xlim,[0 0],'Color',[0.8,0.8,0.9],'LineWidth',4);
                  line(xlim,[psthImDur psthImDur],'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title(sprintf('%s %s autocorrelation, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},group{item_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                  clear figData
                  figData.x = t;
                  figData.y = t;
                  figData.z = Cgram;
                  drawnow;
                  saveFigure(outDir,sprintf('autocorrel_TF_%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                end
              end
              if isfield(plotSwitch,'spikeAutocorrelByItem') && plotSwitch.spikeAutocorrelByItem
                fh = figure();
                lineProps.width = 3;
                lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
                mseb(shifts,spectra,specErrs,lineProps);
                legend(group);
                xlabel('time (ms)');
                ylabel('correlation');
                title(sprintf('%s %s autocorrelation %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = shifts;
                figData.y = spectra;
                figData.e = specErrs;
                drawnow;
                saveFigure(outDir,sprintf('autocorrel_%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              end
            end
          end
        end
      end
    end
    % clean up temporary variables and restore stable variables
    if strcmp(calcSwitchNames{calc_i},'induced')
      spikesByCategoryBinned = spikesByCategoryBinnedEvokedTmp;
      lfpByCategory =  lfpByCategoryEvokedTmp;
    end
  end
end

%% image-wise time-frequency plots, spikes and lfp
% todo: put in option for just preferred image
tfCalcSwitchNames = {'evokedImageTF', 'inducedImageTF'};
tfCalcSwitchTitleSuffixes = {'',', induced'}; % appended to titles
tfCalcSwitchFnameSuffixes = {'','_induced'}; % appended to filenames
tfCalcSwitches = [calcSwitch.evokedImageTF, calcSwitch.inducedImageTF];
for calc_i = 1:length(tfCalcSwitches)
  if tfCalcSwitches(calc_i)
    if strcmp(tfCalcSwitchNames{calc_i},'inducedImageTF')
      if calcSwitch.spikeTimes
        disp('induced TF not yet implemented for spike times; change to spike bins using calcSwitch.spikeTimes = 0');
        continue
      end
      disp('switching spikes and lfps to induced for imagewise tf computation');
      spikesByEventBinnedEvokedTmp = spikesByEventBinned;
      spikesByEventBinned = spikesByEventBinnedInduced;
      lfpByEventEvokedTmp = lfpByEvent;
      lfpByEvent = lfpByEventInduced;
    end
    if isfield(plotSwitch,'spikeSpectraTfByImage') && plotSwitch.spikeSpectraTfByImage
      for channel_i = 1:length(channelNames)
        for image_i = 1:length(eventLabels)
          % todo: put in dB conversion and specgramrowave option
          for unit_i =1:length(channelUnitNames{channel_i})
            if length(channelUnitNames{channel_i}) == 2 && unit_i ==1 %if no isolated spike, don't separate unsorted and MUA
              continue
            end
            if calcSwitch.spikeTimes
              [S,t,f,~,Serr]=mtspecgrampt(spikesByEventForTF{image_i}{channel_i}{unit_i},movingWin,chronuxParams); %last optional param not used: fscorr
            else
              [S,t,f,~,Serr]=mtspecgrampb(spikesByEventBinned{image_i}{channel_i}{unit_i}',movingWin,chronuxParams); %last optional param not used: fscorr
            end
            %todo: add Serr plot, significance plot, effect size plot, or similar
            %todo: add support for image calc groups, incl 'pref' wildcard
            t = t - lfpAlignParams.msPreAlign;
            f = 1000*f;
            fh = figure();
            imagesc(t,f,S'); %TODO: fix to match cat version, which is correct!!
            set(gca,'Ydir','normal'); % also need 'Yscale,'log'??
            xlabel('Time'); % todo: units
            ylabel('Frequency (Hz)');
            c = colorbar();
            ylabel(c,'Power'); % todo: check units
            hold on
            draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
            draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
            title(sprintf('%s %s Time-Frequency, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},eventLabels{image_i},tfCalcSwitchTitleSuffixes{calc_i}));
            clear figData
            figData.x = t;
            figData.y = f;
            figData.z = S';
            drawnow;
            saveFigure(outDir,sprintf('TF_%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},eventLabels{image_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
            
            if calcSwitch.trialMeanSpectra
              if calcSwitch.spikeTimes
                [S,t,f,~,Serr]=mtspecgrampt(allSpikesByImageForTF{image_i}{channel_i}{unit_i},movingWin,chronuxParams); %last optional param not used: fscorr
              else
                [S,t,f,~,Serr]=mtspecgrampb(mean(spikesByEventBinned{image_i}{channel_i}{unit_i},1)',movingWin,chronuxParams); %last optional param not used: fscorr
              end
              %todo: add Serr plot, significance plot, effect size plot, or similar
              %todo: add support for image calc groups, incl 'pref' wildcard
              t = t - lfpAlignParams.msPreAlign;
              f = 1000*f;
              fh = figure();
              imagesc(t,f,S'); %TODO: fix to match cat version, which is correct!!
              set(gca,'Ydir','normal'); % also need 'Yscale,'log'??
              xlabel('Time'); % todo: units
              ylabel('Frequency (Hz)');
              c = colorbar();
              ylabel(c,'Power'); % todo: check units
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              title(sprintf('%s %s Time-Frequency, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},eventLabels{image_i},tfCalcSwitchTitleSuffixes{calc_i}));
              clear figData
              figData.x = t;
              figData.y = f;
              figData.z = S';
              drawnow;
              saveFigure(outDir,sprintf('TF_%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},eventLabels{image_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
            end
          end
        end
      end
    end
    if isfield(plotSwitch,'lfpSpectraTfByImage') && plotSwitch.lfpSpectraTfByImage
      for channel_i = 1:length(channelNames)
        %todo: add Serr plot, significance plot, effect size plot, or similar
        %todo: add support for image calc groups, incl 'pref' wildcard
        [S,t,f,Serr]=mtspecgramc(squeeze(lfpByEvent{image_i}(1,channel_i,:,:))',movingWin,chronuxParams);
        t = t - lfpAlignParams.msPreAlign;
        f = 1000*f;
        fh = figure();
        imagesc(t,f,S');
        set(gca,'Ydir','normal'); % also need 'Yscale,'log'??
        xlabel('Time'); % todo: units
        ylabel('Frequency (Hz)');
        c = colorbar();
        ylabel(c,'Power'); % todo: check units
        hold on
        draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
        draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
        title(sprintf('%s mean LFP Time-Frequency, %s%s',channelNames{channel_i},eventLabels{image_i},tfCalcSwitchTitleSuffixes{calc_i}));
        clear figData
        figData.x = t;
        figData.y = f;
        figData.z = S;
        drawnow;
        saveFigure(outDir,sprintf('TF_mean_LFP_%s_%s%s_Run%s',channelNames{channel_i},eventLabels{image_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
      end
    end
    if strcmp(tfCalcSwitchNames{calc_i},'inducedImageTF')
      spikesByEventBinned = spikesByEventBinnedEvokedTmp;
      lfpByEvent = lfpByEventEvokedTmp;
    end
  end
end

%% category-wise time-frequency plots, spikes
tfCalcSwitchNames = {'evokedCatSpikeTF', 'inducedCatSpikeTF'};
tfCalcSwitchTitleSuffixes = {'',', induced'}; % appended to titles
tfCalcSwitchFnameSuffixes = {'','_induced'}; % appended to filenames
tfCalcSwitches = [calcSwitch.evokedCatSpikeTF, calcSwitch.inducedCatSpikeTF];
for calc_i = 1:length(tfCalcSwitches)
  if tfCalcSwitches(calc_i)
    if strcmp(tfCalcSwitchNames{calc_i},'inducedCatSpikeTF')
      if calcSwitch.spikeTimes
        disp('induced TF not yet implemented for spike times; change to spike bins using calcSwitch.spikeTimes = 0');
        continue
      end
      disp('switching spikes and lfps to induced for categorywise coupling computation');
      spikesByCategoryBinnedEvokedTmp = spikesByCategoryBinned;
      spikesByCategoryBinned = spikesByCategoryBinnedInduced;
      lfpByCategoryEvokedTmp = lfpByCategory;
      lfpByCategory = lfpByCategoryInduced;
    end
    
    if isfield(plotSwitch,'tfSpectraByCategory') && plotSwitch.tfSpectraByCategory %plots spectra individually by category/channel, and as nChannels x nCategories in group subplot
      for group_i = 1:length(analysisGroups.tfSpectraByCategory.groups)
        group = analysisGroups.tfSpectraByCategory.groups{group_i};
        groupName = analysisGroups.tfSpectraByCategory.names{group_i};
        for channel_i = 1:length(channelNames)
          if channel_i == 1
            muaTfFig = figure();
          end
          for item_i = 1:length(group) %note that cat_i here refers to an index into the group
            % spikes
            for unit_i = 1:length(channelUnitNames{channel_i})
              if length(channelUnitNames{channel_i}) == 2 && unit_i == 1
                continue
              end
              if calcSwitch.spikeTimes
                [S,t,f,~,Serr] = mtspecgrampt(spikesByCategoryForTF{catInds.(group{item_i})}{channel_i}{unit_i},movingWin,chronuxParams); %last optional param not used: fscorr
              else
                [S,t,f,~,Serr] = mtspecgrampb(spikesByCategoryBinned{catInds.(group{item_i})}{channel_i}{unit_i}',movingWin,chronuxParams); %last optional param not used: fscorr
              end
              t = t - lfpAlignParams.msPreAlign;
              f = 1000*f;
              totalPower = sum(S,2)';
              fh = figure();
              if specgramRowAve
                for i = 1:size(S,2)
                  S(:,i) = S(:,i)/mean(S(:,i));
                end
                imagesc(t,f,S'); axis xy; c = colorbar();
                ylabel(c,'Row-Normalized Power');
              else
                S = 10*log10(S);
                imagesc(t,f,S'); axis xy; c = colorbar();
                ylabel(c,'Power (dB)');
              end
              xlabel('Time (ms)');
              ylabel('Frequency (Hz)');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              yyaxis right
              plot(t,totalPower,'Color',[0.8,0.8,0.9],'LineWidth',4);
              ylabel('Integrated Power');
              hold off
              title(sprintf('%s %s Time-Frequency, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},categoryList{catInds.(group{item_i})},tfCalcSwitchTitleSuffixes{calc_i}));
              clear figData
              figData.x = t;
              figData.y = f;
              figData.z = S;
              drawnow;
              saveFigure(outDir,sprintf('TF_%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},categoryList{catInds.(group{item_i})},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              
              % contribute mua to shared figure
              if unit_i == length(channelUnitNames{channel_i})
                figure(muaTfFig);
                subplot(length(channelNames),length(group),length(group)*(channel_i-1)+item_i);
                forTitle = sprintf('%s MUA TF, %s%s',categoryList{catInds.(group{item_i})},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i});
                imagesc(t,f,S'); axis xy; c = colorbar(); ylabel(c,'Power (dB)'); %todo: fix unit when row-normalizing
                hold on
                yyaxis right
                plot(t,totalPower,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title(sprintf('%s LFP Time-Frequency, %s%s',channelNames{channel_i},eventLabels{image_i},tfCalcSwitchTitleSuffixes{calc_i}));
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                title(forTitle);
              end
              
              if calcSwitch.trialMeanSpectra
                if calcSwitch.spikeTimes
                  [S,t,f,~,Serr] = mtspecgrampt(allSpikesByCategoryForTF{catInds.(group{item_i})}{channel_i}{unit_i},movingWin,chronuxParams); %last optional param not used: fscorr
                else
                  [S,t,f,~,Serr] = mtspecgrampb(mean(spikesByCategoryBinned{catInds.(group{item_i})}{channel_i}{unit_i},1)',movingWin,chronuxParams); %last optional param not used: fscorr
                end
                t = t - lfpAlignParams.msPreAlign;
                f = 1000*f;
                totalPower = sum(S,2)';
                fh = figure();
                if specgramRowAve
                  for i = 1:size(S,2)
                    S(:,i) = S(:,i)/mean(S(:,i));
                  end
                  imagesc(t,f,S'); axis xy; c = colorbar();
                  ylabel(c,'Row-Normalized Power');
                else
                  S = 10*log10(S);
                  imagesc(t,f,S'); axis xy; c = colorbar();
                  ylabel(c,'Power (dB)');
                end
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                yyaxis right
                plot(t,totalPower,'Color',[0.8,0.8,0.9],'LineWidth',4);
                ylabel('Integrated Power');
                hold off
                title(sprintf('%s %s mean Time-Frequency, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},categoryList{catInds.(group{item_i})},tfCalcSwitchTitleSuffixes{calc_i}));
                clear figData
                figData.x = t;
                figData.y = f;
                figData.z = S;
                drawnow;
                saveFigure(outDir,sprintf('TF_mean_%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},categoryList{catInds.(group{item_i})},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              end
            end
            
            if calcSwitch.trialMeanSpectra
              [S,t,f]=mtspecgramc(squeeze(mean(lfpByCategory{catInds.(group{item_i})}(1,channel_i,:,:),3)),movingWin,chronuxParams);
              t = t - lfpAlignParams.msPreAlign;
              f = 1000*f;
              totalPower = sum(S,2)';
              
              fh = figure();
              if specgramRowAve
                for i = 1:size(S,2)
                  S(:,i) = S(:,i)/mean(S(:,i));
                end
                imagesc(t,f,S'); axis xy; c = colorbar();
                ylabel(c,'Row-Normalized Power');
              else
                S = 10*log10(S);
                imagesc(t,f,S'); axis xy; c = colorbar();
                ylabel(c,'Power (dB)');
              end
              xlabel('Time (ms)');
              ylabel('Frequency (Hz)');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              yyaxis right
              plot(t,totalPower,'Color',[0.8,0.8,0.9],'LineWidth',4);
              ylabel('Integrated Power');
              hold off
              title(sprintf('%s mean LFP Time-Frequency, %s%s',channelNames{channel_i},categoryList{catInds.(group{item_i})},tfCalcSwitchTitleSuffixes{calc_i}));
              clear figData
              figData.x = t;
              figData.y = f;
              figData.z = S';
              drawnow;
              saveFigure(outDir,sprintf('TF_mean_%s_LFP_%s%s_Run%s',channelNames{channel_i},categoryList{catInds.(group{item_i})},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
            end
          end
          drawnow;
          if channel_i == length(lfpChannels) && item_i == length(group)
            figure(muaTfFig);
            saveFigure(outDir,sprintf('TF_MUA_%s_%s%s.fig',channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i}), figData, figStruct, figStruct.figTag );
          end
        end
      end
    end
    if strcmp(tfCalcSwitchNames{calc_i},'inducedCatTF')
      spikesByCategoryBinned = spikesByCategoryBinnedEvokedTmp;
      lfpByCategory = lfpByCategoryEvokedTmp;
    end
  end
end

%% LFP Catagory comparison
tfCalcSwitchNames = {'evokedCatLFPTF', 'inducedCatLFPTF'};
tfCalcSwitchTitleSuffixes = {'',', induced'}; % appended to titles
tfCalcSwitchFnameSuffixes = {'','_induced'}; % appended to filenames
tfCalcSwitches = [calcSwitch.evokedCatLFPTF, calcSwitch.inducedCatLFPTF];
for calc_i = 1:length(tfCalcSwitches)
  if tfCalcSwitches(calc_i)
    if strcmp(tfCalcSwitchNames{calc_i},'inducedCatLFPTF')
      if calcSwitch.spikeTimes
        disp('induced TF not yet implemented for spike times; change to spike bins using calcSwitch.spikeTimes = 0');
        continue
      end
      disp('switching spikes and lfps to induced for categorywise coupling computation');
      spikesByCategoryBinnedEvokedTmp = spikesByCategoryBinned;
      spikesByCategoryBinned = spikesByCategoryBinnedInduced;
      lfpByCategoryEvokedTmp = lfpByCategory;
      lfpByCategory = lfpByCategoryInduced;
    end
    
    if isfield(plotSwitch,'tfSpectraByCategory') && plotSwitch.tfSpectraByCategory %plots spectra individually by category/channel, and as nChannels x nCategories in group subplot
      for group_i = 1:length(analysisGroups.tfSpectraByCategory.groups)
        group = analysisGroups.tfSpectraByCategory.groups{group_i};
        groupName = analysisGroups.tfSpectraByCategory.names{group_i};
        for channel_i = 1:length(channelNames)
          forTitle = sprintf('LFP Time Frequency Analysis, %s, Run%s',channelNames{channel_i},runNum);
          if channel_i == 1
            lfpTfFig = figure('Name',forTitle,'NumberTitle','off');
          end
          for item_i = 1:length(group) %note that cat_i here refers to an index into the group
            % LFP
            [S,t,f] = mtspecgramc(squeeze(lfpByCategory{catInds.(group{item_i})}(1,channel_i,:,:))',movingWin,chronuxParams);
            t = t - lfpAlignParams.msPreAlign;
            f = 1000*f;
            totalPower = sum(S,2)';
            
            % contribute to shared figure
            figure(lfpTfFig);
            subplot(length(channelNames),length(group),length(group)*(channel_i-1)+item_i);
            imagesc(t,f,S'); axis xy; c = colorbar(); ylabel(c,'Power (dB)'); %todo: fix unit when row-normalizing
            hold on
            draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
            draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
            yyaxis right
            plot(t,totalPower,'Color',[0.8,0.8,0.9],'LineWidth',4);
            xlabel('Time (ms)');
            yyaxis left
            ylabel('Frequency (Hz)');
            title(forTitle);
            
            if calcSwitch.trialMeanSpectra
              [S,t,f]=mtspecgramc(squeeze(mean(lfpByCategory{catInds.(group{item_i})}(1,channel_i,:,:),3)),movingWin,chronuxParams);
              t = t - lfpAlignParams.msPreAlign;
              f = 1000*f;
              totalPower = sum(S,2)';
              
              fh = figure();
              if specgramRowAve
                for i = 1:size(S,2)
                  S(:,i) = S(:,i)/mean(S(:,i));
                end
                imagesc(t,f,S'); axis xy; c = colorbar();
                ylabel(c,'Row-Normalized Power');
              else
                S = 10*log10(S);
                imagesc(t,f,S'); axis xy; c = colorbar();
                ylabel(c,'Power (dB)');
              end
              xlabel('Time (ms)');
              ylabel('Frequency (Hz)');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              yyaxis right
              plot(t,totalPower,'Color',[0.8,0.8,0.9],'LineWidth',4);
              ylabel('Integrated Power');
              hold off
              title(sprintf('%s mean LFP Time-Frequency, %s%s',channelNames{channel_i},categoryList{catInds.(group{item_i})},tfCalcSwitchTitleSuffixes{calc_i}));
              clear figData
              figData.x = t;
              figData.y = f;
              figData.z = S';
              drawnow;
              saveFigure(outDir,sprintf('TF_mean_%s_LFP_%s%s_Run%s',channelNames{channel_i},categoryList{catInds.(group{item_i})},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
            end
          end
          drawnow;
          if channel_i == length(lfpChannels) && item_i == length(group)
            figure(lfpTfFig);
            saveFigure(outDir,sprintf('TF_LFP_%s_%s%s.fig',channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i}), figData, figStruct, figStruct.figTag );
          end
        end
      end
    end
    if strcmp(tfCalcSwitchNames{calc_i},'inducedCatTF')
      spikesByCategoryBinned = spikesByCategoryBinnedEvokedTmp;
      lfpByCategory = lfpByCategoryEvokedTmp;
    end
  end
end

%% coupling across channels and modalities (units/unsorted/mua, fields)
tfCalcSwitchNames = {'evokedCoupling', 'inducedCoupling'}; % todo: make separate switch for coherence?
tfCalcSwitchTitleSuffixes = {'',', induced'}; % appended to titles
tfCalcSwitchFnameSuffixes = {'','_induced'}; % appended to filenames
tfCalcSwitches = [calcSwitch.evokedCoupling, calcSwitch.inducedCoupling];
for calc_i = 1:length(tfCalcSwitches)
  if tfCalcSwitches(calc_i)
    if strcmp(tfCalcSwitchNames{calc_i},'inducedCoupling')
      if calcSwitch.spikeTimes
        disp('induced TF not yet implemented for spike times; change to spike bins using calcSwitch.spikeTimes = 0');
        continue
      end
      disp('switching spikes and lfps to induced for categorywise coupling computation');
      spikesByCategoryBinnedEvokedTmp = spikesByCategoryBinned;
      spikesByCategoryBinned = spikesByCategoryBinnedInduced;
      lfpByCategoryEvokedTmp = lfpByCategory;
      lfpByCategory = lfpByCategoryInduced;
    end
    
    %spike -- field time-domain coupling
    for group_i = 1:length(analysisGroups.coherenceByCategory.groups) %todo: make own analysis group
      group = analysisGroups.coherenceByCategory.groups{group_i};
      groupName = analysisGroups.coherenceByCategory.names{group_i};
      for tmp_channel_i = 1:length(channelNames)
        for tmp_channel2_i = tmp_channel_i:length(lfpChannels) %make sure we only calculate once for each pair
          for order_i = 1:2
            if order_i == 2 && tmp_channel_i == tmp_channel2_i %intrachannel, then two orders are equivalent
              continue
            end
            if order_i == 1 % tmp_channel2_i spike -- tmp_channel_i field
              channel_i = tmp_channel_i;
              channel2_i = tmp_channel2_i;
            else
              channel_i = tmp_channel2_i; % tmp_channel2_i spike -- tmp_channel_i field
              channel2_i = tmp_channel_i;
            end
            % channel2_i spike -- channel_i field (order set above)
            for unit_i = 1:length(channelUnitNames{channel2_i})
              if length(channelUnitNames{channel2_i}) == 2 && unit_i == 1 %skip unsorted if no isolated units; just do mua
                continue
              end
              for item_i = 1:length(group)
                if isfield(catInds,group{item_i})
                  lfpByItem = lfpByCategory;
                  spikesByItemBinned = spikesByCategoryBinned;
                  %spikesByItemForTF = spikesByCategoryForTF;
                  itemNum = catInds.(group{item_i});
                else
                  lfpByItem = lfpByEvent;
                  spikesByItemBinned = spikesByEventBinned;
                  %spikesByItemForTF = spikesByEventForTF;
                  itemNum = imInds.(group{item_i});
                end
                if calcSwitch.spikeTimes
                  disp('requested time-domain spike-field analysis. STA lfp not yet implemented.');
                  continue
                else
                  correlParams.useSmoothing = [0 1];
                  [Cgram, C, shiftList, confCgram, confC] = correlogram(squeeze(lfpByItem{itemNum}(1,channel_i,:,:)),...
                    spikesByItemBinned{itemNum}{channel2_i}{unit_i}, correlParams);
                  correlParams.useSmoothing = [0 0];
                end
                Cgram = Cgram(lfpPaddedBy+1:end-lfpPaddedBy,lfpPaddedBy+1:end-lfpPaddedBy);
                confCgram = confCgram(lfpPaddedBy+1:end-lfpPaddedBy,lfpPaddedBy+1:end-lfpPaddedBy);
                if item_i == 1
                  spectra = zeros(length(group),length(C));
                  specErrs = zeros(length(group),length(C),2); % support for asymmetric error bars, as from jacknife
                  shifts =  zeros(length(group),length(C));
                end
                spectra(item_i,:) = C';
                specErrs(item_i,:,:) = confC;
                shifts(item_i,:) = shiftList';
                t = -psthPre:psthImDur+psthPost;
                fh = figure();
                imagesc(t,t,Cgram'); axis xy
                xlabel(sprintf('Time, %s field (ms)', channelNames{channel_i}));
                ylabel(sprintf('Time, %s %s (ms)', channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i}));
                c = colorbar();
                ylabel(c,'Correlation');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                line(xlim,[0 0],'Color',[0.8,0.8,0.9],'LineWidth',4);
                line(xlim,[psthImDur psthImDur],'Color',[0.8,0.8,0.9],'LineWidth',4);
                line(xlim,ylim,'Color','w','LineWidth',0.5);
                title(sprintf('%s %s - %s field correlation, %s%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},group{item_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = t;
                figData.y = 1;
                figData.z = Cgram';
                drawnow;
                saveFigure(outDir,sprintf('correl_TF_%s_%s_-%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              end
              fh = figure();
              lineProps.width = 3;
              lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
              mseb(shifts,spectra,specErrs,lineProps);
              legend(group);
              xlabel('time (ms)');
              ylabel('correlation');
              title(sprintf('%s %s - %s field correlation %s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = shifts;
              figData.y = spectra;
              figData.e = specErrs;
              drawnow;
              saveFigure(outDir,sprintf('correl_%s_%s_%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              if closeFig
                close(fh);
              end
            end
          end
        end
      end
    end
    
    %spike -- field coherency, within and across channel, frequency domain
    for group_i = 1:length(analysisGroups.coherenceByCategory.groups)
      group = analysisGroups.coherenceByCategory.groups{group_i};
      groupName = analysisGroups.coherenceByCategory.names{group_i};
      for tmp_channel_i = 1:length(channelNames)
        for tmp_channel2_i = tmp_channel_i:length(lfpChannels) %make sure we only calculate once for each pair
          for order_i = 1:2
            if order_i == 2 && tmp_channel_i == tmp_channel2_i %intrachannel, then two orders are equivalent
              continue
            end
            if order_i == 1 % tmp_channel2_i spike -- tmp_channel_i field
              channel_i = tmp_channel_i;
              channel2_i = tmp_channel2_i;
            else
              channel_i = tmp_channel2_i; % tmp_channel2_i spike -- tmp_channel_i field
              channel2_i = tmp_channel_i;
            end
            % channel2_i spike -- channel_i field (order set above)
            for unit_i = 1:length(channelUnitNames{channel2_i})
              if length(channelUnitNames{channel2_i}) == 2 && unit_i == 1 %skip unsorted if no isolated units; just do mua
                continue
              end
              for item_i = 1:length(group)
                if isfield(catInds,group{item_i})
                  lfpByItem = lfpByCategory;
                  spikesByItemBinned = spikesByCategoryBinned;
                  spikesByItemForTF = spikesByCategoryForTF;
                  itemNum = catInds.(group{item_i});
                else
                  lfpByItem = lfpByEvent;
                  spikesByItemBinned = spikesByEventBinned;
                  spikesByItemForTF = spikesByEventForTF;
                  itemNum = imInds.(group{item_i});
                end
                if calcSwitch.spikeTimes
                  [C,phi,~,~,~,f,~,confC,phistd, Cerr] = callChronuxCoherency(@coherencycpt,squeeze(lfpByItem{itemNum}(1,channel_i,:,:))',...
                    spikesByItemForTF{itemNum}{channel2_i}{unit_i},chronuxParams); %can have time grid as additional arg
                else
                  [C,phi,~,~,~,f,~,confC,phistd, Cerr] = callChronuxCoherency(@coherencycpb,squeeze(lfpByItem{itemNum}(1,channel_i,:,:))',...
                    spikesByItemBinned{itemNum}{channel2_i}{unit_i}',chronuxParams);
                end
                if item_i == 1
                  spectra = zeros(length(group),length(C));
                  specErrs = zeros(length(group),length(C),2); % support for asymmetric error bars, as from jacknife
                  phases = zeros(length(group),length(C));
                  phaseErrs = zeros(length(group),length(C));
                  freqs=  zeros(length(group),length(C));
                end
                spectra(item_i,:) = C';
                phases(item_i,:) = phi';
                if calcSwitch.useJacknife
                  specErrs(item_i,:,1) = Cerr(1,:)-C;
                  specErrs(item_i,:,2) = C-Cerr(2,:);
                else
                  specErrs(item_i,:,:) = confC;
                end
                phaseErrs(item_i,:) = phistd';
                freqs(item_i,:) = 1000*f';
              end
              % note: this correction maintains consistency with X-Y coherence, where X comes first in the call to coherency, but so we can say X=spike,
              % Y=field, while chronux requires us to put field first in the function call
              phases = -1*phases;
              
              fh = figure();
              lineProps.width = 3;
              lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
              mseb(freqs,spectra,specErrs,lineProps);
              legend(group);
              xlabel('frequency (Hz)');
              ylabel('coherency');
              title(sprintf('%s %s - %s field coherence%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = freqs;
              figData.y = spectra;
              figData.e = specErrs;
              drawnow;
              saveFigure(outDir,sprintf('coh_%s_%s_%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              if closeFig
                close(fh);
              end
              
              %phases %todo: convert from std to sterr?
              fh = figure();
              lineProps.width = 3;
              lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
              mseb(freqs,phases,phaseErrs,lineProps);
              legend(group);
              xlabel('frequency (Hz)');
              ylabel('phase');
              title(sprintf('%s %s - %s field phase%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              drawnow;
              saveFigure(outDir,sprintf('phase_%s_%s-%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              if closeFig
                close(fh);
              end
              
              if isfield(plotSwitch,'couplingPhasesUnwrapped') && plotSwitch.couplingPhasesUnwrapped
                %replot the phases with the implicit mod pi operation undone (so easy to see slope across frequencies
                fh = figure();
                lineProps.width = 3;
                lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
                [unwrappedPhases,piMultiples] = unwrapPhases(phases);
                mseb(freqs,unwrappedPhases,phaseErrs,lineProps);
                xlabel('frequency (Hz)');
                ylabel('phase');
                hold on
                lineHandles = findobj(gca,'Type','line');  %note: these two lines exclude the pi multiples from legend
                legendHandles = lineHandles(1:size(phases,1));
                plot(xlim(),repmat(pi*piMultiples,2,1),'color','k');
                legend(legendHandles,group);
                title(sprintf('%s %s - %s field phase%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                drawnow;
                saveFigure(outDir,sprintf('phase_unwrapped_%s_%s-%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
              end
              
              if isfield(plotSwitch,'couplingPhasesAsOffsets') && plotSwitch.couplingPhasesAsOffsets
                %replot the unwrapped phases as latency differences
                fh = figure();
                lineProps.width = 3;
                lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
                [unwrappedPhases,piMultiples] = unwrapPhases(phases);
                mseb(freqs,1000*(1/(2*pi))*unwrappedPhases./freqs,1000*(1/(2*pi))*phaseErrs./freqs,lineProps);
                xlabel('frequency (Hz)');
                ylabel('offset (ms)');
                hold on
                lineHandles = findobj(gca,'Type','line');  %note: these two lines exclude the pi multiples from legend
                legendHandles = lineHandles(1:size(phases,1));
                plot(freqs(1,:),1000*(1/(2*pi))*kron(piMultiples',1./freqs(1,:)),'color','k');
                legend(legendHandles,group);
                title(sprintf('%s %s - %s field offset%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                drawnow;
                saveFigure(outDir,sprintf('phase_unwrapped_latency_%s_%s-%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
              end
              
              if isfield(plotSwitch,'couplingPhasesPolar') && plotSwitch.couplingPhasesPolar
                %replot the unwrapped phases as latency differences
                fh = figure();
                [unwrappedPhases,~] = unwrapPhases(phases);
                legendHandles = gobjects(length(group),1);
                for item_i = 1:length(group)
                  ah = polarplot(unwrappedPhases(item_i,:),freqs(item_i,:),'color',analysisGroups.coherenceByCategory.colors{group_i}{item_i},'linewidth',3);
                  hold on
                  polarplot(unwrappedPhases(item_i,:)+phaseErrs(item_i,:),freqs(item_i,:),'color',analysisGroups.coherenceByCategory.colors{group_i}{item_i},'linewidth',1);
                  polarplot(unwrappedPhases(item_i,:)-phaseErrs(item_i,:),freqs(item_i,:),'color',analysisGroups.coherenceByCategory.colors{group_i}{item_i},'linewidth',1);
                  legendHandles(item_i) = ah;
                end
                legend(legendHandles,group);
                title(sprintf('%s %s - %s field phases%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                drawnow;
                saveFigure(outDir,sprintf('phase_unwrapped_polar_%s_%s-%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
              end
              
              if calcSwitch.meanEvokedTF
                for item_i = 1:length(group)
                  if isfield(catInds,group{item_i})
                    lfpByItem = lfpByCategory;
                    spikesByItemBinned = spikesByCategoryBinned;
                    if calcSwitch.spikeTimes
                      allSpikesByItemForTF = allSpikesByCategoryForTF;
                    end
                    itemNum = catInds.(group{item_i});
                  else
                    lfpByItem = lfpByEvent;
                    spikesByItemBinned = spikesByEventBinned;
                    if calcSwitch.spikeTimes
                      allSpikesByItemForTF = allSpikesByImageForTF;
                    end
                    itemNum = imInds.(group{item_i});
                  end
                  if calcSwitch.spikeTimes
                    [C,phi,~,~,~,f,~,confC,phistd, Cerr] = callChronuxCoherency(@coherencycpt,squeeze(mean(lfpByItem{itemNum}(1,channel_i,:,:),3)),...
                      allSpikesByItemForTF{itemNum}{channel2_i}{unit_i},chronuxParams); %can have time grid as additional arg
                  else
                    [C,phi,~,~,~,f,~,confC,phistd, Cerr] = callChronuxCoherency(@coherencycpb,squeeze(mean(lfpByItem{itemNum}(1,channel_i,:,:),3)),...
                      mean(spikesByItemBinned{itemNum}{channel2_i}{unit_i},1)',chronuxParams);
                  end
                  if item_i == 1
                    spectra = zeros(length(group),length(C));
                    specErrs = zeros(length(group),length(C),2); % support for asymmetric error bars, as from jacknife
                    phases = zeros(length(group),length(C));
                    phaseErrs = zeros(length(group),length(C));
                    freqs=  zeros(length(group),length(C));
                  end
                  spectra(item_i,:) = C';
                  phases(item_i,:) = phi';
                  if calcSwitch.useJacknife
                    specErrs(item_i,:,1) = Cerr(1,:)-C;
                    specErrs(item_i,:,2) = C-Cerr(2,:);
                  else
                    specErrs(item_i,:,:) = confC;
                  end
                  phaseErrs(item_i,:) = phistd';
                  freqs(item_i,:) = 1000*f';
                end
                % note: this correction maintains consistency with X-Y coherence, where X comes first in the call to coherency, but so we can say X=spike,
                % Y=field, while chronux requires us to put field first in the function call
                phases = -1*phases;
                
                fh = figure();
                lineProps.width = 3;
                lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
                mseb(freqs,spectra,specErrs,lineProps);
                legend(group);
                xlabel('frequency (Hz)');
                ylabel('coherency');
                title(sprintf('Trial mean %s %s - %s field coherence%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = freqs;
                figData.y = spectra;
                figData.e = specErrs;
                drawnow;
                saveFigure(outDir,sprintf('coh_mean_%s_%s_%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
                
                %phases %todo: convert from std to sterr?
                fh = figure();
                lineProps.width = 3;
                lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
                mseb(freqs,phases,phaseErrs,lineProps);
                legend(group);
                xlabel('frequency (Hz)');
                ylabel('phase');
                title(sprintf('Trial mean %s %s - %s field phase%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                drawnow;
                saveFigure(outDir,sprintf('phase_mean_%s_%s-%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
                
                if isfield(plotSwitch,'couplingPhasesUnwrapped') && plotSwitch.couplingPhasesUnwrapped
                  %replot the phases with the implicit mod pi operation undone (so easy to see slope across frequencies
                  fh = figure();
                  lineProps.width = 3;
                  lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
                  [unwrappedPhases,piMultipliers] = unwrapPhases(phases);
                  mseb(freqs,unwrappedPhases,phaseErrs,lineProps);
                  legend(group);
                  xlabel('frequency (Hz)');
                  ylabel('phase');
                  hold on
                  line(repmat(xlim(),length(piMultipliers),1)',repmat(piMultipliers',1,2)','color','k');
                  title(sprintf('Trial mean %s %s - %s field phase%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                  drawnow;
                  saveFigure(outDir,sprintf('phase_mean_unwrapped_%s_%s-%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                  if closeFig
                    close(fh);
                  end
                end
                
                if isfield(plotSwitch,'couplingPhasesAsOffsets') && plotSwitch.couplingPhasesAsOffsets
                  %replot the unwrapped phases as latency differences
                  fh = figure();
                  lineProps.width = 3;
                  lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
                  [unwrappedPhases,piMultipliers] = unwrapPhases(phases);
                  mseb(freqs,(1/(2*pi))*unwrappedPhases./freqs,phaseErrs,lineProps);
                  legend(group);
                  xlabel('frequency (Hz)');
                  ylabel('offset (ms)');
                  hold on
                  plot(repmat(xlim(),length(piMultipliers),1),[(1/(2*pi))*piMultipliers'/freqs(1,1), (1/(2*pi))*piMultipliers'/freqs(1,end)],'color','k');
                  title(sprintf('Trial mean %s %s - %s field offsets%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                  drawnow;
                  saveFigure(outDir,sprintf('phase_mean_unwrapped_latency_%s_%s-%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                  if closeFig
                    close(fh);
                  end
                end
              end
            end
          end
        end
      end
    end
    
    
    %time-frequency spike -- field coherency
    for group_i = 1:length(analysisGroups.coherenceByCategory.groups)
      group = analysisGroups.coherenceByCategory.groups{group_i};
      groupName = analysisGroups.coherenceByCategory.names{group_i};
      for tmp_channel_i = 1:length(channelNames)
        for tmp_channel2_i = tmp_channel_i:length(lfpChannels) %make sure we only calculate once for each pair
          for order_i = 1:2
            if order_i == 2 && tmp_channel_i == tmp_channel2_i %intrachannel, then two orders are equivalent
              continue
            end
            if order_i == 1 % tmp_channel2_i spike -- tmp_channel_i field
              channel_i = tmp_channel_i;
              channel2_i = tmp_channel2_i;
            else
              channel_i = tmp_channel_i; % tmp_channel2_i spike -- tmp_channel_i field
              channel2_i = tmp_channel2_i;
            end
            % channel2_i spike -- channel_i field (order set above)
            for unit_i = 1:length(channelUnitNames{channel2_i})
              if length(channelUnitNames{channel2_i}) == 2 && unit_i == 1 %skip unsorted if no isolated units; just do mua
                continue
              end
              for item_i = 1:length(group)
                if isfield(catInds,group{item_i})
                  lfpByItem = lfpByCategory;
                  spikesByItemBinned = spikesByCategoryBinned;
                  spikesByItemForTF = spikesByCategoryForTF;
                  itemNum = catInds.(group{item_i});
                else
                  lfpByItem = lfpByEvent;
                  spikesByItemBinned = spikesByEventBinned;
                  spikesByItemForTF = spikesByEventForTF;
                  itemNum = imInds.(group{item_i});
                end
                if calcSwitch.spikeTimes
                  [C,phi,~,~,~,t,f,~,confC,phistd,Cerr] = callChronuxCoherency(@cohgramcpt,squeeze(lfpByItem{itemNum}(1,channel_i,:,:))',...
                    spikesByItemForTF{itemNum}{channel2_i}{unit_i}',movingWin,chronuxParams);
                else
                  [C,phi,~,~,~,t,f,~,confC,phistd,Cerr] = callChronuxCoherency(@cohgramcpb,squeeze(lfpByItem{itemNum}(1,channel_i,:,:))',...
                    spikesByItemBinned{itemNum}{channel2_i}{unit_i}',movingWin,chronuxParams);
                end
                if ~calcSwitch.useJacknife
                  Cerr = ones(2,size(C,1),size(C,2));
                  Cerr(1,:,:) = C+confC; %todo: confirm that this works
                  Cerr(2,:,:) = C-confC;
                end
                
                t = t - lfpAlignParams.msPreAlign;
                f = 1000*f;
                % note: this correction maintains consistency with X-Y coherence, where X comes first in the call to coherency, but so we can say X=spike,
                % Y=field, while chronux requires us to put field first in the function call
                phi = -1*phi;
                phaseErrs = zeros(size(Cerr));
                phaseErrs(1,:,:) = phi+phistd;
                phaseErrs(2,:,:) = phi-phistd;
                
                fh = figure();
                imagesc(t,f,C'); axis xy
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'Coherency');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title(sprintf('%s %s - %s field coherence, %s%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},group{item_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = t;
                figData.y = f;
                figData.z = C';
                drawnow;
                saveFigure(outDir,sprintf('coh_TF_%s_%s_-%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
                
                fh = figure();
                imagesc(t,f,phi'); axis xy
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'phase');
                phasemap;
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title(sprintf('%s %s - %s phase, face%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = t;
                figData.y = f;
                figData.z = phi';
                drawnow;
                saveFigure(outDir,sprintf('phase_TF_%s_%s-%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
                
                
                if isfield(plotSwitch,'tfErrs') && plotSwitch.tfErrs
                  fh = figure();
                  subplot(2,3,1);
                  imagesc(t,f,C'); axis xy
                  xlabel('Time (ms)');
                  ylabel('Frequency (Hz)');
                  c = colorbar();
                  ylabel(c,'Coherency');
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title(sprintf('%s %s - %s field coherence, %s%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},group{item_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                  
                  subplot(2,3,2);
                  imagesc(t,f,squeeze(Cerr(1,:,:))'); axis xy
                  xlabel('Time (ms)');
                  ylabel('Frequency (Hz)');
                  c = colorbar();
                  ylabel(c,'Coherency');
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title('coherence, upper confidence bound','FontSize',18);
                  
                  subplot(2,3,3);
                  imagesc(t,f,squeeze(Cerr(2,:,:))'); axis xy
                  xlabel('Time (ms)');
                  ylabel('Frequency (Hz)');
                  c = colorbar();
                  ylabel(c,'Coherency');
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title('coherence, lower confidence bound','FontSize',18);
                  
                  
                  ax = subplot(2,3,4);
                  imagesc(t,f,phi'); axis xy
                  xlabel('Time (ms)');
                  ylabel('Frequency (Hz)');
                  c = colorbar();
                  ylabel(c,'phase');
                  m = phasemap();
                  colormap(ax,m);
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title(sprintf('%s %s - %s LFP phase, face%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                  
                  ax = subplot(2,3,5);
                  imagesc(t,f,squeeze(phistd(2,:,:))'); axis xy
                  xlabel('Time (ms)');
                  ylabel('Frequency (Hz)');
                  c = colorbar();
                  ylabel(c,'phase');
                  m = phasemap();
                  colormap(ax,m);
                  ylabel(c,'phase std');
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title('phase, upper confidence bound','FontSize',18);
                  
                  ax = subplot(2,3,6);
                  imagesc(t,f,squeeze(phistd(2,:,:))'); axis xy
                  xlabel('Time (ms)');
                  ylabel('Frequency (Hz)');
                  c = colorbar();
                  ylabel(c,'phase');
                  m = phasemap();
                  colormap(ax,m);
                  ylabel(c,'phase std');
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title('phase, lower confidence bound','FontSize',18);
                  
                  saveFigure(outDir,sprintf('coh_TF_errs_%s_%s_-%s_LFP_%s%s_Run%s',channelNames{channel2_i},channelUnitNames{channel2_i}{unit_i},channelNames{channel_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                  if closeFig
                    close(fh);
                  end
                end
              end
            end
          end
        end
      end
    end
    
    % field -- field time domain coupling
    for group_i = 1:length(analysisGroups.coherenceByCategory.groups) %todo: make own analysis group
      group = analysisGroups.coherenceByCategory.groups{group_i};
      groupName = analysisGroups.coherenceByCategory.names{group_i};
      for channel_i = 1:length(channelNames)
        for channel2_i = channel_i:length(lfpChannels) %these two lines make sure we only calculate once for each pair
          if channel2_i > channel_i                    %avoids figuring out matlab behavior when channel_i == length(lfpChannels)
            for item_i = 1:length(group)
              if isfield(catInds,group{item_i})
                lfpByItem = lfpByCategory;
                itemNum = catInds.(group{item_i});
              else
                lfpByItem = lfpByEvent;
                itemNum = imInds.(group{item_i});
              end
              [Cgram, C, shiftList, confCgram, confC] = correlogram(squeeze(lfpByItem{itemNum}(1,channel_i,:,:)),...
                squeeze(lfpByItem{itemNum}(1,channel2_i,:,:)), correlParams);
              Cgram = Cgram(lfpPaddedBy+1:end-lfpPaddedBy,lfpPaddedBy+1:end-lfpPaddedBy);
              confCgram = confCgram(lfpPaddedBy+1:end-lfpPaddedBy,lfpPaddedBy+1:end-lfpPaddedBy);
              if item_i == 1
                spectra = zeros(length(group),length(C));
                specErrs = zeros(length(group),length(C),2); % support for asymmetric error bars, as from jacknife
                shifts =  zeros(length(group),length(C));
              end
              spectra(item_i,:) = C';
              specErrs(item_i,:,:) = confC;
              shifts(item_i,:) = shiftList';
              t = -psthPre:psthImDur+psthPost;
              fh = figure();
              imagesc(t,t,Cgram'); axis xy
              xlabel(sprintf('Time, %s field (ms)', channelNames{channel_i}));
              ylabel(sprintf('Time, %s field (ms)', channelNames{channel2_i}));
              c = colorbar();
              ylabel(c,'Correlation');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              line(xlim,[0 0],'Color',[0.8,0.8,0.9],'LineWidth',4);
              line(xlim,[psthImDur psthImDur],'Color',[0.8,0.8,0.9],'LineWidth',4);
              line(xlim,ylim,'Color','w','LineWidth',0.5);
              title(sprintf('%s field - %s field correlation, %s%s',channelNames{channel_i},channelNames{channel2_i},group{item_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = t;
              figData.y = 1;
              figData.z = Cgram';
              drawnow;
              saveFigure(outDir,sprintf('correl_TF_%s_LFP-%s_LFP_%s%s_Run%s',channelNames{channel_i},channelNames{channel2_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              if closeFig
                close(fh);
              end
            end
            fh = figure();
            lineProps.width = 3;
            lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
            mseb(shifts,spectra,specErrs,lineProps);
            legend(group);
            xlabel('time (ms)');
            ylabel('correlation');
            title(sprintf('%s field - %s field correlation %s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = shifts;
            figData.y = spectra;
            figData.e = specErrs;
            drawnow;
            saveFigure(outDir,sprintf('correl_%s_LFP_%s_LFP_%s%s_Run%s',channelNames{channel_i},channelNames{channel2_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
            if closeFig
              close(fh);
            end
          end
        end
      end
    end
    
    
    % field -- field coherency, frequency domain
    for group_i = 1:length(analysisGroups.coherenceByCategory.groups)
      group = analysisGroups.coherenceByCategory.groups{group_i};
      groupName = analysisGroups.coherenceByCategory.names{group_i};
      for channel_i = 1:length(channelNames)
        for channel2_i = channel_i:length(lfpChannels) %these two lines make sure we only calculate once for each pair
          if channel2_i > channel_i
            for item_i = 1:length(group)
              if isfield(catInds,group{item_i})
                lfpByItem = lfpByCategory;
                itemNum = catInds.(group{item_i});
              else
                lfpByItem = lfpByEvent;
                itemNum = imInds.(group{item_i});
              end
              
              [C,phi,~,~,~,f,confC,phistd, Cerr] = callChronuxCoherency(@coherencyc,squeeze(lfpByItem{itemNum}(1,channel_i,:,:))',...
                squeeze(lfpByItem{itemNum}(1,channel2_i,:,:))',chronuxParams);
              if item_i == 1
                spectra = zeros(length(group),length(C));
                specErrs = zeros(length(group),length(C),2); % support for asymmetric error bars, as from jacknife
                phases = zeros(length(group),length(C));
                phaseErrs = zeros(length(group),length(C));
                freqs=  zeros(length(group),length(C));
              end
              spectra(item_i,:) = C';
              phases(item_i,:) = phi';
              if calcSwitch.useJacknife
                specErrs(item_i,:,1) = Cerr(1,:)-C;
                specErrs(item_i,:,2) = C-Cerr(2,:);
              else
                specErrs(item_i,:,:) = confC;
              end
              phaseErrs(item_i,:) = phistd';
              freqs(item_i,:) = 1000*f';
            end
            
            fh = figure();
            lineProps.width = 3;
            lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
            mseb(freqs,spectra,specErrs,lineProps);
            legend(group);
            xlabel('frequency (Hz)');
            ylabel('coherency');
            title(sprintf('%s field - %s field coherence%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            clear figData
            figData.x = freqs;
            figData.y = spectra;
            figData.e = specErrs;
            drawnow;
            saveFigure(outDir,sprintf('coh_%s_LFP_%s_LFP_%s%s_Run%s',channelNames{channel_i},channelNames{channel2_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
            if closeFig
              close(fh);
            end
            
            %phases %todo: convert from std to sterr?
            fh = figure();
            lineProps.width = 3;
            lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
            mseb(freqs,phases,phaseErrs,lineProps);
            legend(group);
            xlabel('frequency (Hz)');
            ylabel('phase');
            title(sprintf('%s field - %s field phase%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
            drawnow;
            saveFigure(outDir,sprintf('phase_%s_LFP-%s_LFP_%s%s_Run%s',channelNames{channel_i},channelNames{channel2_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
            if closeFig
              close(fh);
            end
            
            if calcSwitch.meanEvokedTF
              for item_i = 1:length(group)
                if isfield(catInds,group{item_i})
                  lfpByItem = lfpByCategory;
                  itemNum = catInds.(group{item_i});
                else
                  lfpByItem = lfpByEvent;
                  itemNum = imInds.(group{item_i});
                end
                [C,phi,~,~,~,f,confC,phistd, Cerr] = callChronuxCoherency(@coherencyc,squeeze(mean(lfpByItem{itemNum}(1,channel_i,:,:),3)),...
                  squeeze(mean(lfpByItem{itemNum}(1,channel2_i,:,:),3)),chronuxParams);
                if item_i == 1
                  spectra = zeros(length(group),length(C));
                  specErrs = zeros(length(group),length(C),2); % support for asymmetric error bars, as from jacknife
                  phases = zeros(length(group),length(C));
                  phaseErrs = zeros(length(group),length(C));
                  freqs=  zeros(length(group),length(C));
                end
                spectra(item_i,:) = C';
                phases(item_i,:) = phi';
                if calcSwitch.useJacknife
                  specErrs(item_i,:,1) = Cerr(1,:)-C;
                  specErrs(item_i,:,2) = C-Cerr(2,:);
                else
                  specErrs(item_i,:,:) = confC;
                end
                phaseErrs(item_i,:) = phistd';
                freqs(item_i,:) = 1000*f';
              end
              
              fh = figure();
              lineProps.width = 3;
              lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
              mseb(freqs,spectra,specErrs,lineProps);
              legend(group);
              xlabel('frequency (Hz)');
              ylabel('coherency');
              title(sprintf('Trial mean %s field - %s field coherence%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = freqs;
              figData.y = spectra;
              figData.e = specErrs;
              drawnow;
              saveFigure(outDir,sprintf('coh_mean_%s_LFP-%s_LFP_%s%s_Run%s',channelNames{channel_i},channelNames{channel2_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              if closeFig
                close(fh);
              end
              
              %phases %todo: convert from std to sterr?
              fh = figure();
              lineProps.width = 3;
              lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
              mseb(freqs,phases,phaseErrs,lineProps);
              legend(group);
              xlabel('frequency (Hz)');
              ylabel('phase');
              title(sprintf('Trial mean %s field - %s field phase%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              drawnow;
              saveFigure(outDir,sprintf('phase_mean_%s_LFP-%s_LFP_%s%s_Run%s',channelNames{channel_i},channelNames{channel2_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              if closeFig
                close(fh);
              end
            end
          end
        end
      end
    end
    
    % time-frequency field --field coherency
    for group_i = 1:length(analysisGroups.coherenceByCategory.groups)
      group = analysisGroups.coherenceByCategory.groups{group_i};
      groupName = analysisGroups.coherenceByCategory.names{group_i};
      for channel_i = 1:length(channelNames)
        for channel2_i = channel_i:length(lfpChannels) %these two lines make sure we only calculate once for each pair
          if channel2_i > channel_i
            for item_i = 1:length(group)
              if isfield(catInds,group{item_i})
                lfpByItem = lfpByCategory;
                itemNum = catInds.(group{item_i});
              else
                lfpByItem = lfpByEvent;
                itemNum = imInds.(group{item_i});
              end
              [C,phi,~,~,~,t,f,confC,phistd,Cerr]=callChronuxCoherency(@cohgramc,squeeze(lfpByItem{itemNum}(1,channel_i,:,:))',...
                squeeze(lfpByItem{itemNum}(1,channel2_i,:,:))',movingWin,chronuxParams);
              if ~calcSwitch.useJacknife
                Cerr = ones(2,size(C,1),size(C,2));
                Cerr(1,:,:) = C+confC; %todo: confirm that this works
                Cerr(2,:,:) = C-confC;
              end
              t = t - lfpAlignParams.msPreAlign;
              f = 1000*f;
              phaseErrs = zeros(size(Cerr));
              phaseErrs(1,:,:) = phi+phistd;
              phaseErrs(2,:,:) = phi-phistd;
              
              fh = figure();
              imagesc(t,f,C'); axis xy
              xlabel('Time (ms)');
              ylabel('Frequency (Hz)');
              c = colorbar();
              ylabel(c,'Coherency');
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              title(sprintf('%s field - %s field coherence, %s%s',channelNames{channel_i},channelNames{channel2_i},group{item_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = t;
              figData.y = f;
              figData.z = C';
              drawnow;
              saveFigure(outDir,sprintf('coh_TF_%s_LFP-%s_LFP_%s%s_Run%s',channelNames{channel_i},channelNames{channel2_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              if closeFig
                close(fh);
              end
              
              fh = figure();
              imagesc(t,f,phi'); axis xy
              xlabel('Time (ms)');
              ylabel('Frequency (Hz)');
              c = colorbar();
              ylabel(c,'phase');
              phasemap;
              hold on
              draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
              draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
              title(sprintf('%s field - %s field phase%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = t;
              figData.y = f;
              figData.z = phi';
              drawnow;
              saveFigure(outDir,sprintf('phase_TF_%s_LFP-%s_LFP_%s%s_Run%s',channelNames{channel_i},channelNames{channel2_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              if closeFig
                close(fh);
              end
              
              
              if isfield(plotSwitch,'tfErrs') && plotSwitch.tfErrs
                fh = figure();
                subplot(2,3,1);
                imagesc(t,f,C'); axis xy
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'Coherency');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title(sprintf('%s field - %s field coherence, %s%s',channelNames{channel_i},channelNames{channel2_i},group{item_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                
                subplot(2,3,2);
                imagesc(t,f,squeeze(Cerr(1,:,:))'); axis xy
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'Coherency');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title('coherence, upper confidence bound','FontSize',18);
                
                subplot(2,3,3);
                imagesc(t,f,squeeze(Cerr(2,:,:))'); axis xy
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'Coherency');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title('coherence, lower confidence bound','FontSize',18);
                
                
                ax = subplot(2,3,4);
                imagesc(t,f,phi'); axis xy
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'phase');
                m = phasemap();
                colormap(ax,m);
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title(sprintf('%s field - %s field phase, face%s',channelNames{channel_i},channelNames{channel2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                
                ax =subplot(2,3,5);
                imagesc(t,f,squeeze(phistd(1,:,:))'); axis xy
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'phase');
                m = phasemap();
                colormap(ax,m);
                ylabel(c,'phase std');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title('phase, upper confidence bound','FontSize',18);
                
                ax =subplot(2,3,6);
                imagesc(t,f,squeeze(phistd(2,:,:))'); axis xy
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'phase');
                m = phasemap();
                colormap(ax,m);
                ylabel(c,'phase std');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title('phase, lower confidence bound','FontSize',18);
                
                saveFigure(outDir,sprintf('coh_TF_errs_%s_LFP-%s_LFP_%s%s_Run%s',channelNames{channel_i},channelNames{channel2_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
              end
            end
          end
        end
      end
    end
    
    %spike -- spike time-domain coupling
    for group_i = 1:length(analysisGroups.coherenceByCategory.groups) %todo: make own analysis group
      group = analysisGroups.coherenceByCategory.groups{group_i};
      groupName = analysisGroups.coherenceByCategory.names{group_i};
      for channel_i = 1:length(channelNames)
        for channel2_i = channel_i:length(lfpChannels) %make sure we only calculate once for each pair
          for unit_i = 1:length(channelUnitNames{channel_i})
            if length(channelUnitNames{channel_i}) == 2 && unit_i == 1 %skip unsorted if no isolated unit
              continue
            end
            for unit2_i = 1:length(channelUnitNames{channel2_i})
              if length(channelUnitNames{channel2_i}) == 2 && unit2_i == 1 %skip unsorted if no isolated unit
                continue
              end
              %if intrachannel, check that pair is unique and non-self
              if channel_i == channel2_i && unit_i >= unit2_i
                continue
              end
              for item_i = 1:length(group)
                if isfield(catInds,group{item_i})
                  spikesByItemBinned = spikesByCategoryBinned;
                  %spikesByItemForTF = spikesByCategoryForTF;
                  itemNum = catInds.(group{item_i});
                else
                  spikesByItemBinned = spikesByEventBinned;
                  %spikesByItemForTF = spikesByEventForTF;
                  itemNum = imInds.(group{item_i});
                end
                if calcSwitch.spikeTimes
                  disp('requested time-domain spike-field analysis. STA lfp not yet implemented.');
                  continue
                else
                  correlParams.useSmoothing = [1 1];
                  [Cgram, C, shiftList, confCgram, confC] = correlogram(spikesByItemBinned{itemNum}{channel_i}{unit_i},...
                    spikesByItemBinned{itemNum}{channel2_i}{unit2_i}, correlParams);
                  correlParams.useSmoothing = [0 0];
                end
                Cgram = Cgram(lfpPaddedBy+1:end-lfpPaddedBy,lfpPaddedBy+1:end-lfpPaddedBy);
                confCgram = confCgram(lfpPaddedBy+1:end-lfpPaddedBy,lfpPaddedBy+1:end-lfpPaddedBy);
                if item_i == 1
                  spectra = zeros(length(group),length(C));
                  specErrs = zeros(length(group),length(C),2); % support for asymmetric error bars, as from jacknife
                  shifts =  zeros(length(group),length(C));
                end
                spectra(item_i,:) = C';
                specErrs(item_i,:,:) = confC;
                shifts(item_i,:) = shiftList';
                t = -psthPre:psthImDur+psthPost;
                fh = figure();
                imagesc(t,t,Cgram'); axis xy
                xlabel(sprintf('Time, %s %s (ms)', channelNames{channel_i},channelUnitNames{channel_i}{unit_i}));
                ylabel(sprintf('Time, %s %s (ms)', channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i}));
                c = colorbar();
                ylabel(c,'Correlation');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                line(xlim,[0 0],'Color',[0.8,0.8,0.9],'LineWidth',4);
                line(xlim,[psthImDur psthImDur],'Color',[0.8,0.8,0.9],'LineWidth',4);
                line(xlim,ylim,'Color','w','LineWidth',0.5);
                title(sprintf('%s %s - %s %s correlation, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},group{item_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = t;
                figData.y = 1;
                figData.z = Cgram';
                drawnow;
                saveFigure(outDir,sprintf('correl_TF_%s_%s_-%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
              end
              fh = figure();
              lineProps.width = 3;
              lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
              mseb(shifts,spectra,specErrs,lineProps);
              legend(group);
              xlabel('time (ms)');
              ylabel('correlation');
              title(sprintf('%s %s - %s %s correlation %s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = shifts;
              figData.y = spectra;
              figData.e = specErrs;
              drawnow;
              saveFigure(outDir,sprintf('correl_%s_%s_%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel_i},channelUnitNames{channel2_i}{unit2_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              if closeFig
                close(fh);
              end
            end
          end
        end
      end
    end
    
    % spike --spike coherency, frequency domain
    for group_i = 1:length(analysisGroups.coherenceByCategory.groups)
      group = analysisGroups.coherenceByCategory.groups{group_i};
      groupName = analysisGroups.coherenceByCategory.names{group_i};
      for channel_i = 1:length(channelNames)
        for channel2_i = channel_i:length(lfpChannels) %make sure we only calculate once for each pair
          for unit_i = 1:length(channelUnitNames{channel_i})
            if length(channelUnitNames{channel_i}) == 2 && unit_i == 1
              continue
            end
            for unit2_i = 1:length(channelUnitNames{channel2_i})
              if length(channelUnitNames{channel2_i}) == 2 && unit2_i == 1
                continue
              end
              %if intrachannel, check that pair is unique and non-self
              if channel_i == channel2_i && unit_i >= unit2_i
                continue
              end
              for item_i = 1:length(group)
                if isfield(catInds,group{item_i})
                  spikesByItemBinned = spikesByCategoryBinned;
                  spikesByItemForTF = spikesByCategoryForTF;
                  itemNum = catInds.(group{item_i});
                else
                  spikesByItemBinned = spikesByEventBinned;
                  spikesByItemForTF = spikesByEventForTF;
                  itemNum = imInds.(group{item_i});
                end
                if calcSwitch.spikeTimes
                  [C,phi,~,~,~,f,~,confC,phistd, Cerr] = callChronuxCoherency(@coherencypt,spikesByItemForTF{itemNum}{channel_i}{unit_i},...
                    spikesByItemForTF{itemNum}{channel2_i}{unit2_i},chronuxParams); %can have time grid as additional arg
                else
                  [C,phi,~,~,~,f,~,confC,phistd, Cerr] = callChronuxCoherency(@coherencypb,spikesByItemBinned{itemNum}{channel_i}{unit_i}',...
                    spikesByItemBinned{itemNum}{channel2_i}{unit2_i}',chronuxParams);
                end
                if item_i == 1
                  spectra = zeros(length(group),length(C));
                  specErrs = zeros(length(group),length(C),2); % support for asymmetric error bars, as from jacknife
                  phases = zeros(length(group),length(C));
                  phaseErrs = zeros(length(group),length(C));
                  freqs=  zeros(length(group),length(C));
                end
                spectra(item_i,:) = C';
                phases(item_i,:) = phi';
                if calcSwitch.useJacknife
                  specErrs(item_i,:,1) = Cerr(1,:)-C;
                  specErrs(item_i,:,2) = C-Cerr(2,:);
                else
                  specErrs(item_i,:,:) = confC;
                end
                phaseErrs(item_i,:) = phistd';
                freqs(item_i,:) = 1000*f';
              end
              
              fh = figure();
              lineProps.width = 3;
              lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
              mseb(freqs,spectra,specErrs);
              legend(group);
              xlabel('frequency (Hz)');
              ylabel('coherency');
              title(sprintf('%s %s - %s %s coherence%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              clear figData
              figData.x = freqs;
              figData.y = spectra;
              figData.e = specErrs;
              drawnow;
              saveFigure(outDir,sprintf('coh_%s_%s_%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              if closeFig
                close(fh);
              end
              
              %phases %todo: convert from std to sterr?
              fh = figure();
              lineProps.width = 3;
              lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
              mseb(freqs,phases,phaseErrs);
              legend(group);
              xlabel('frequency (Hz)');
              ylabel('phase');
              title(sprintf('%s %s - %s %s phase%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
              drawnow;
              saveFigure(outDir,sprintf('phase_%s_%s-%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
              if closeFig
                close(fh);
              end
              
              if calcSwitch.meanEvokedTF
                for item_i = 1:length(group)
                  if isfield(catInds,group{item_i})
                    spikesByItemBinned = spikesByCategoryBinned;
                    if calcSwitch.spikeTimes
                      allSpikesByItemForTF = allSpikesByCategoryForTF;
                    end
                    itemNum = catInds.(group{item_i});
                  else
                    spikesByItemBinned = spikesByEventBinned;
                    if calcSwitch.spikeTimes
                      allSpikesByItemForTF = allSpikesByImageForTF;
                    end
                    itemNum = imInds.(group{item_i});
                  end
                  if calcSwitch.spikeTimes
                    [C,phi,~,~,~,f,~,confC,phistd, Cerr] = callChronuxCoherency(@coherencypt,allSpikesByItemForTF{itemNum}{channel_i}{unit_i},...
                      allSpikesByItemForTF{itemNum}{channel2_i}{unit2_i},chronuxParams); %can have time grid as additional arg
                  else
                    [C,phi,~,~,~,f,~,confC,phistd, Cerr] = callChronuxCoherency(@coherencypb,squeeze(mean(spikesByItemBinned{itemNum}{channel_i}{unit_i},3)),...
                      squeeze(mean(spikesByItemBinned{itemNum}{channel2_i}{unit2_i},3)),chronuxParams);
                  end
                  if item_i == 1
                    spectra = zeros(length(group),length(C));
                    specErrs = zeros(length(group),length(C),2); % support for asymmetric error bars, as from jacknife
                    phases = zeros(length(group),length(C));
                    phaseErrs = zeros(length(group),length(C));
                    freqs=  zeros(length(group),length(C));
                  end
                  spectra(item_i,:) = C';
                  phases(item_i,:) = phi';
                  if calcSwitch.useJacknife
                    specErrs(item_i,:,1) = Cerr(1,:)-C;
                    specErrs(item_i,:,2) = C-Cerr(2,:);
                  else
                    specErrs(item_i,:,:) = confC;
                  end
                  phaseErrs(item_i,:) = phistd';
                  freqs(item_i,:) = 1000*f';
                end
                
                fh = figure();
                lineProps.width = 3;
                lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
                mseb(freqs,spectra,specErrs);
                legend(group);
                xlabel('frequency (Hz)');
                ylabel('coherency');
                title(sprintf('Trial mean %s %s - %s %s coherence%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = freqs;
                figData.y = spectra;
                figData.e = specErrs;
                drawnow;
                saveFigure(outDir,sprintf('coh_mean_%s_%s_%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
                
                %phases %todo: convert from std to sterr?
                fh = figure();
                lineProps.width = 3;
                lineProps.col = analysisGroups.coherenceByCategory.colors{group_i};
                mseb(freqs,phases,phaseErrs);
                legend(group);
                xlabel('frequency (Hz)');
                ylabel('phase');
                title(sprintf('Trial mean %s %s - %s %s phase%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                drawnow;
                saveFigure(outDir,sprintf('phase_mean_%s_%s-%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},groupName,tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
              end
            end
          end
        end
      end
    end
    
    %time-frequency spike --spike coherency
    for group_i = 1:length(analysisGroups.coherenceByCategory.groups)
      group = analysisGroups.coherenceByCategory.groups{group_i};
      groupName = analysisGroups.coherenceByCategory.names{group_i};
      for channel_i = 1:length(channelNames)
        for channel2_i = channel_i:length(lfpChannels) %make sure we only calculate once for each pair
          for unit_i = 1:length(channelUnitNames{channel_i})
            if length(channelUnitNames{channel_i}) == 2 && unit_i == 1
              continue
            end
            for unit2_i = 1:length(channelUnitNames{channel2_i})
              if length(channelUnitNames{channel2_i}) == 2 && unit2_i == 1
                continue
              end
              %if intrachannel, check that pair is unique and non-self
              if channel_i == channel2_i && unit_i >= unit2_i
                continue
              end
              for item_i = 1:length(group)
                if isfield(catInds,group{item_i})
                  spikesByItemBinned = spikesByCategoryBinned;
                  spikesByItemForTF = spikesByCategoryForTF;
                  itemNum = catInds.(group{item_i});
                else
                  spikesByItemBinned = spikesByEventBinned;
                  spikesByItemForTF = spikesByEventForTF;
                  itemNum = imInds.(group{item_i});
                end
                if calcSwitch.spikeTimes
                  [C,phi,~,~,~,t,f,~,confC,phistd]=callChronuxCoherency(@cohgrampt,spikesByItemForTF{itemNum}{channel_i}{unit_i},...
                    spikesByItemForTF{itemNum}{channel2_i}{unit2_i}, movingWin,chronuxParams);
                else
                  [C,phi,~,~,~,t,f,~,confC,phistd]=callChronuxCoherency(@cohgrampb,spikesByItemBinned{itemNum}{channel_i}{unit_i}',...
                    spikesByItemBinned{itemNum}{channel2_i}{unit2_i}', movingWin,chronuxParams);
                end
                if ~calcSwitch.useJacknife
                  Cerr = ones(2,size(C,1),size(C,2));
                  Cerr(1,:,:) = C+confC; %todo: confirm that this works
                  Cerr(2,:,:) = C-confC;
                end
                t = t - lfpAlignParams.msPreAlign;
                f = 1000*f;
                phaseErrs = zeros(size(Cerr));
                phaseErrs(1,:,:) = phi+phistd;
                phaseErrs(2,:,:) = phi-phistd;
                
                fh = figure();
                imagesc(t,f,C'); axis xy
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'Coherency');
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title(sprintf('%s %s - %s %s coherence, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},group{item_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = t;
                figData.y = f;
                figData.z = C';
                drawnow;
                saveFigure(outDir,sprintf('coh_TF_%s_%s_-%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
                
                fh = figure();
                imagesc(t,f,phi'); axis xy
                xlabel('Time (ms)');
                ylabel('Frequency (Hz)');
                c = colorbar();
                ylabel(c,'phase');
                phasemap;
                hold on
                draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                title(sprintf('%s %s - %s %s phase, face%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                clear figData
                figData.x = t;
                figData.y = f;
                figData.z = phi';
                drawnow;
                saveFigure(outDir,sprintf('phase_TF_%s_%s-%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                if closeFig
                  close(fh);
                end
                
                
                if isfield(plotSwitch,'tfErrs') && plotSwitch.tfErrs
                  fh = figure();
                  subplot(2,3,1);
                  imagesc(t,f,C'); axis xy
                  xlabel('Time (ms)');
                  ylabel('Frequency (Hz)');
                  c = colorbar();
                  ylabel(c,'Coherency');
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title(sprintf('%s %s - %s %s coherence, %s%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},group{item_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                  
                  subplot(2,3,2);
                  imagesc(t,f,squeeze(Cerr(1,:,:))'); axis xy
                  xlabel('Time (ms)');
                  ylabel('Frequency (Hz)');
                  c = colorbar();
                  ylabel(c,'Coherency');
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title('coherence, upper confidence bound','FontSize',18);
                  
                  subplot(2,3,3);
                  imagesc(t,f,squeeze(Cerr(2,:,:))'); axis xy
                  xlabel('Time (ms)');
                  ylabel('Frequency (Hz)');
                  c = colorbar();
                  ylabel(c,'Coherency');
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title('coherence, lower confidence bound','FontSize',18);
                  
                  
                  ax = subplot(2,3,4);
                  imagesc(t,f,phi'); axis xy
                  xlabel('Time (ms)');
                  ylabel('Frequency (Hz)');
                  c = colorbar();
                  ylabel(c,'phase');
                  m = phasemap();
                  colormap(ax,m);
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title(sprintf('%s %s - %s %s phase, face%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},tfCalcSwitchTitleSuffixes{calc_i}),'FontSize',18);
                  
                  ax = subplot(2,3,5);
                  imagesc(t,f,squeeze(phistd(1,:,:))'); axis xy
                  xlabel('Time (ms)');
                  ylabel('Frequency (Hz)');
                  c = colorbar();
                  ylabel(c,'phase');
                  m = phasemap();
                  colormap(ax,m);
                  ylabel(c,'phase std');
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title('phase, upper confidence bound','FontSize',18);
                  
                  ax = subplot(2,3,6);
                  imagesc(t,f,squeeze(phistd(2,:,:))'); axis xy
                  xlabel('Time (ms)');
                  ylabel('Frequency (Hz)');
                  c = colorbar();
                  ylabel(c,'phase');
                  m = phasemap();
                  colormap(ax,m);
                  ylabel(c,'phase std');
                  hold on
                  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
                  title('phase, lower confidence bound','FontSize',18);
                  saveFigure(outDir,sprintf('coh_TF_errs_%s_%s_-%s_%s_%s%s_Run%s',channelNames{channel_i},channelUnitNames{channel_i}{unit_i},channelNames{channel2_i},channelUnitNames{channel2_i}{unit2_i},group{item_i},tfCalcSwitchFnameSuffixes{calc_i},runNum), figData, figStruct, figStruct.figTag );
                  if closeFig
                    close(fh);
                  end
                end
              end
            end
          end
        end
      end
    end
    
    % clean up temporary variables and restore stable variables
    if strcmp(tfCalcSwitchNames{calc_i},'inducedCoupling')
      spikesByCategoryBinned = spikesByCategoryBinnedEvokedTmp;
      lfpByCategory =  lfpByCategoryEvokedTmp;
    end
  end
end

errorMsgRun = 'None';
save(analysisOutFilename, 'errorMsgRun', '-append')

end
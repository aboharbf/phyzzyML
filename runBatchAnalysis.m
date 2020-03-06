function [analysisOutFilename] = runBatchAnalysis(inputs)
%%Unpack Inputs
analysisLog = struct();
load(inputs{1},'outputDir')   % Path to spikes.
spikeDataBank = saveSpikeDataBank([], [], 'load', outputDir); %Spike variables
load(inputs{1})   % Non spike variables, including path to spikes.

%Overwrite switches with what is currently in file
load(fullfile(outputDir, 'batchAnalysisParams.mat'));

%% Analyses
runList = fields(spikeDataBank);

if ~isfield(spikeDataBank.(runList{end}), 'stimPresCount')
  spikeDataBank = stimulusStatistics(spikeDataBank);
  %Save modified struct.
  saveEnv()
end

if plotSwitch.stimPresCount
  disp('plot stimPresCount')
  %Plot as an image the 'daysSinceLastPres' matrix, stimuli on Y, date on
  %X.will likely require a figure callback for making Y axis labels
  %legible.
  %At the top of the plot, I will have a single vector with 'days since
  %last recording' to highlight days after long breaks.
end

if 1%~exist('unitCounts','var')
  [trueCellInd, trueCellInfo, unitCounts, resultTable, nullCells] = trueCellCount(cellCountParams.batchRunxls, cellCountParams.recordingLogxls);
  for run_ind = 1:length(runList)
    sessionName = extractBetween(runList{run_ind}, 2, length(runList{run_ind}));
    dataInd = (strcmp(sessionName,trueCellInfo(:,1)));
    spikeDataBank.(runList{run_ind}).gridHoles  = trueCellInfo(dataInd, 2);
    spikeDataBank.(runList{run_ind}).recDepth  = trueCellInfo(dataInd, 3);
  end
%   saveEnv()
end

plotPies = 0;
if plotPies
  unitCount = sum(resultTable{1}{3,1});
  taskModUnits = sum(resultTable{1}{3,2});
  socSigUnits = sum(resultTable{1}{3,3});
  
  nonTaskMod = (unitCount-taskModUnits)/unitCount;
  TaskMod = (taskModUnits-socSigUnits)/unitCount;
  socSig = (socSigUnits)/unitCount;
  sigFrac = [nonTaskMod, TaskMod, socSig];
  sigTitle = {sprintf('Non Task Modulated (n = %d, %s%%)', nonTaskMod*unitCount, num2str(round(nonTaskMod,3)*100)),...
    sprintf('Task Modulated (n = %d, %s%%)', round(TaskMod*unitCount), num2str(round(TaskMod,3)*100)),...
    sprintf('Socially activated (n = %d, %s%%)', socSig*unitCount, num2str(round(socSig,3)*100)),...
    };
  patchColor = {[0.25, 0.25, 0.25]; [0, 0.4470, 0.7410]; 	[0.6350, 0.0780, 0.1840]};
  
  for fig_i = 1:length(sigFrac)
    figure('NumberTitle', 'off','units','normalized','outerposition',[0 0 params.plotSize]);
    h = pie(sigFrac(1:fig_i), sigTitle(1:fig_i));
    for patch_i = 1:fig_i
      obj_i = ((patch_i - 1) * 2) + 1;
      h(obj_i).FaceColor = patchColor{patch_i};
    end
    saveFigure([], num2str(fig_i), [], 0, 1, 0, {''})
  end
end

meanPSTHParams.tTestTable = resultTable;

% Remove repeated Runs
if ~calcSwitch.excludeRepeats && isfield(analysisLog, 'repeatsExcluded')
  error('Repeats already excluded in this spikeDataBank, delete and restart or change parameters')
elseif calcSwitch.excludeRepeats && ~isfield(analysisLog, 'repeatsExcluded')
  % Exclude repeated recordings at the same site.
  for run_ind = 1:length(runList)
    sessionName = extractBetween(runList{run_ind}, 2, length(runList{run_ind}));
    validInd = trueCellInd(strcmp(sessionName,trueCellInfo(:,1)));
    if sum(validInd) == 0
      % Remove entire field if all channels are repeated recordings.
      spikeDataBank = rmfield(spikeDataBank,(runList{run_ind}));
    else
      for event_ind = 1:length(spikeDataBank.(runList{run_ind}).spikesByEvent)
        % Remove individual channel info where one of the channels is
        % recording new units.
        spikeDataBank.(runList{run_ind}).spikesByEvent{event_ind} = spikeDataBank.(runList{run_ind}).spikesByEvent{event_ind}(validInd);
      end
      spikeDataBank.(runList{run_ind}).gridHoles(validInd) = spikeDataBank.(runList{run_ind}).gridHoles(validInd);
      spikeDataBank.(runList{run_ind}).recDepth(validInd) = spikeDataBank.(runList{run_ind}).recDepth(validInd);
    end
  end
  %Save modified struct.
  runListIndex = unitCounts.Exclude;
  analysisLog.repeatsExcluded = 1;
  saveEnv()
else
  runListIndex = unitCounts.nonExclude;
end

%% Combine PSTH across all runs for a particular stimulus.
%This will crash if the PSTHs aren't the same length.
if plotSwitch.meanPSTH 
  [stimPSTH, meanPSTHStruct] = meanPSTH(spikeDataBank, meanPSTHParams);
  saveEnv()
end

%% Combine PSTH across all runs for a particular event.

if plotSwitch.subEventPSTH %&& ~exist('meanPSTHStruct','var')
  [stimPSTH, meanPSTHStruct] = meanPSTH(spikeDataBank, meanPSTHParams);
  saveEnv()
end

%% Gather information on frame rates
if plotSwitch.frameFiringRates %&& ~exist('frameFiringStruct','var')
  [spikeDataBank, frameFiringStruct] = frameFiringRates(spikeDataBank, frameFiringParams);
  saveEnv()
end

%% Check whether the novelty of the runs
if plotSwitch.novelty
  assert(logical(exist('meanPSTHStruct','var')), 'Must run w/ meanPSTH enabled for novelty analysis');
  spikeDataBank = noveltyAnalysis(spikeDataBank, stimPSTH, meanPSTHStruct,[], noveltyParams);
end

%% Perform sliding window ANOVA and Omega calculations
if plotSwitch.slidingWindowANOVA % && ~isfield(spikeDataBank, [Some new field generated by function])
  spikeDataBank = slidingWindowTest(spikeDataBank, slidingTestParams);
  saveEnv()
end

end

%% Functions

function spikeDataBank = stimulusStatistics(spikeDataBank)
% Stimuli Presentation count and 'Novelty' Related Information.
%Code below creates a single large vector of stimuli used, and uses this to
%create individual vectors containing which viewing of the stimulus this
%represent (i.e. 'this run represents the 10th viewing of X.avi'). It also
%appends a dateTime vector to each structure related to how long since the
%last recording day.

runList = fields(spikeDataBank);

%extract the eventIDs field, generate a cell array of unique stimuli
allStimuliVec = struct2cell(structfun(@(x) x.eventIDs, spikeDataBank,'UniformOutput', 0));
allStimuliVec = unique(vertcat(allStimuliVec{:}));

%Produce matrix (N stim * M runs) which gives 0 for non present stim, count of stim presentation otherwise.
stimLogicalArray = zeros(length(allStimuliVec),length(runList));
for run_ind = 1:length(runList)
  stimLogicalArray(:,run_ind) = ismember(allStimuliVec,spikeDataBank.(runList{run_ind}).eventIDs);
end
csStimLogicalArray = cumsum(stimLogicalArray,2);
csStimLogicalArray(~stimLogicalArray) = 0;

% When was a stimulus first seen? Index of runList where first presentation took place.
firstStimPresInd = zeros(length(allStimuliVec),1);
for stim_ind = 1:length(allStimuliVec)
  firstStimPresInd(stim_ind) = find(stimLogicalArray(stim_ind,:),1,'first');
end

% Append a dateTime to each field with the time in days since the last
% recording. Add the relevant slice of the larger csStimLogicalArray.
for run_ind = 1:length(runList)
  spikeDataBank.(runList{run_ind}).stimPresArray = csStimLogicalArray(:,run_ind);
  spikeDataBank.(runList{run_ind}).dateTime = datetime(extractBetween(spikeDataBank.(runList{run_ind}).dateSubject,1,8),'InputFormat','yyyyMMdd');     %Generate 'daysSinceLastRec' for each field.
end

% find unique recording dates, and the distance between them. Add these
% to the spikeDataBank
allDateTimeVec = struct2cell(structfun(@(x) x.dateTime, spikeDataBank,'UniformOutput', 0));
allDateTimeVec = [allDateTimeVec{:}]';
uniqueDateTimeVec = unique(allDateTimeVec);
daysSinceLastRec = [1000; days(diff(uniqueDateTimeVec))];
for run_ind = 1:length(runList)
  spikeDataBank.(runList{run_ind}).daysSinceLastRec = daysSinceLastRec(spikeDataBank.(runList{run_ind}).dateTime == uniqueDateTimeVec);
end

% use the dateTime and stimulus presentation matrix to find out how long
% in days takes place before a particular showing of a stimulus.
daysSinceLastPres = zeros(size(stimLogicalArray));
for stim_ind = 1:size(stimLogicalArray,1)
  presentationInd = logical(stimLogicalArray(stim_ind,:)); %When was the stim shown
  daysSinceLastPres(stim_ind,presentationInd) = [1000; days(diff(allDateTimeVec(presentationInd)))]; %Duration between those dates in days
end

%Return days since last presentation and the stimulus presentation count to spikeDataBank, arranged in way that matches stim table.
for run_ind = 1:size(stimLogicalArray,2)
  [~, big2SmallInd] = ismember(spikeDataBank.(runList{run_ind}).eventIDs,allStimuliVec);
  spikeDataBank.(runList{run_ind}).daysSinceLastPres = daysSinceLastPres(big2SmallInd,run_ind);
  spikeDataBank.(runList{run_ind}).stimPresCount = csStimLogicalArray(big2SmallInd,run_ind);
end

end

function [stimPSTH, meanPSTHStruct] = meanPSTH(spikeDataBank, params)
% Function which combines stimulus presentations across all runs in the spikeDataBank.
% Inputs include spikeDataBank and list of parameters.
disp('Starting mean PSTH Analysis...');
exportFig = params.exportFig;

% Rebuild variables
% extract the eventIDs field, generate a cell array of unique stimuli
allStimuliVec = struct2cell(structfun(@(x) x.eventIDs, spikeDataBank,'UniformOutput', 0));
allStimuliVec = unique(vertcat(allStimuliVec{:}));

% Generate grid for indexing into individual runs and extracting relevant
% PSTHes.
runList = fields(spikeDataBank);
small2BigInd = zeros(length(allStimuliVec), length(runList));
for run_ind = 1:length(runList)
  [~, small2BigInd(:,run_ind)] = ismember(allStimuliVec, spikeDataBank.(runList{run_ind}).eventIDs);
end

stimPresCounts = sum(logical(small2BigInd),2);

% If you only want stim above a certain count, remove here.
if params.plotTopStim
  keepInd = stimPresCounts >= params.topStimPresThreshold;
  stimPresCounts = stimPresCounts(keepInd);
  small2BigInd = small2BigInd(keepInd,:);
  allStimuliVec = allStimuliVec(keepInd,:);
end

animParam.stimParamsFilename = params.stimParamsFilename;
animParam.plotLabels = {'animSocialInteraction', 'animControl'};
[tmp, ~, ~] = plotIndex(allStimuliVec, animParam);
animInd = logical(sum(tmp,2));

%Generate the figure directory
dirTags = {'without Rewards','fixAligned','Normalized','Pres Threshold'};
dirTagSwitch = [params.removeRewardEpoch, params.fixAlign, params.normalize, params.plotTopStim];
params.outputDir = [params.outputDir strjoin(dirTags(logical(dirTagSwitch)),' - ')];

if params.stimInclude == 1
  stimPresCounts = stimPresCounts(animInd);
  small2BigInd = small2BigInd(animInd,:);
  allStimuliVec = allStimuliVec(animInd,:);
  params.outputDir = [params.outputDir '- AnimOnly'];
elseif params.stimInclude == 2
  stimPresCounts = stimPresCounts(~animInd);
  small2BigInd = small2BigInd(~animInd,:);
  allStimuliVec = allStimuliVec(~animInd,:);
  params.outputDir = [params.outputDir '- NoAnim'];
end

[plotMat, briefStimList, params] = plotIndex(allStimuliVec, params);

if params.normalize
  normTag = ' - normalized';
else
  normTag = '';
end

if ~exist(params.outputDir, 'dir')
  mkdir(params.outputDir)
end

% Step 1 - Concatonate PSTHs - iterate across stimuli, grab PSTHes from each run, store
groupingType = {'Unsorted', 'Units', 'MUA'};
dataType = {'PSTH', 'PSTH Err', 'presCount','daysSinceLastPres', 'daysSinceLastRec', 'Run of Day', 'Grid Hole', 'Recording Depth', 'Run Ind'};
% {'Run of Day', 'Grid Hole', 'Run of Day', 'Recording Depth'};

meanPSTHStruct = struct();
meanPSTHStruct.IndStructs{1} = allStimuliVec;
meanPSTHStruct.IndStructs{2} = groupingType;
meanPSTHStruct.IndStructs{3} = dataType;
%stimPSTH{stim,grouping,dataType}
stimPSTH = cell([length(allStimuliVec), length(groupingType), length(dataType)]);
rewardStart = abs(spikeDataBank.(runList{1}).start) + spikeDataBank.(runList{1}).stimDur;

for stim_ind = 1:length(allStimuliVec)
  % Find the runs where the stimulus was present, generate a list of them.
  stimRunIndex = small2BigInd(stim_ind,:);
  psthStimIndex = nonzeros(stimRunIndex);
  psthRunIndex = find(stimRunIndex);
  subRunList = runList(psthRunIndex);
  
  % For all runs containing a particular stimuli, retrieve relevant activity vector in each.
  for subRun_ind = 1:length(subRunList)
    tmpRunStruct = spikeDataBank.(subRunList{subRun_ind});
    for chan_ind = 1:length(tmpRunStruct.psthByImage)
      for unit_ind = 1:length(tmpRunStruct.psthByImage{chan_ind})
        % Retrieve correct PSTH from run
        unitActivity = tmpRunStruct.psthByImage{chan_ind}{unit_ind}(psthStimIndex(subRun_ind),:);
        unitErr = tmpRunStruct.psthErrByImage{chan_ind}{unit_ind}(psthStimIndex(subRun_ind),:);
        
        % Do not include activity which does not meet 1 Hz Threshold.
        if params.rateThreshold
          if mean(unitActivity) < params.rateThreshold
            continue
          end
        end
        
        if params.removeRewardEpoch
          unitActivity = unitActivity(1:rewardStart);
          unitErr = unitErr(1:rewardStart);
          params.psthPost = 0;
        end
        presCount = tmpRunStruct.stimPresCount(psthStimIndex(subRun_ind));
        gridHole = tmpRunStruct.gridHoles;
        depth = tmpRunStruct.recDepth;
        assert(length(gridHole) == length(depth), 'F');
        
        daysSinceLastPres = tmpRunStruct.daysSinceLastPres(psthStimIndex(subRun_ind));
        daysSinceLastRec = tmpRunStruct.daysSinceLastRec;
        % If desired, Z score PSTHs here based on fixation period activity.
        if params.normalize
          tmp = tmpRunStruct.psthByImage{chan_ind}{unit_ind}(:,1:abs(tmpRunStruct.start));
          tmp = reshape(tmp, [size(tmp,1) * size(tmp,2), 1]);
          fixMean = mean(tmp); %Find activity during fixation across all stim.
          fixSD = std(tmp);
          if params.normalize == 1 && fixMean ~= 0 && fixSD ~= 0
            unitActivity = ((unitActivity - fixMean)/fixSD);
            unitErr = ((unitErr - fixMean)/fixSD);
          end
        end
        % Store the generated values into a dataArray
        dataArray = cell(length(dataType),1); % Structure to store and iterate across.
        % dataType = {'PSTH', 'PSTH Err', 'presCount','daysSinceLastPres', 'daysSinceLastRec', 'Run of Day', 'Grid Hole', 'Recording Depth', 'Run Ind'};
        dataArray{1} = unitActivity;
        dataArray{2} = unitErr;
        dataArray{3} = presCount;
        dataArray{4} = daysSinceLastPres;
        dataArray{5} = daysSinceLastRec;
        dataArray{6} = str2double(tmpRunStruct.runNum);
        dataArray{7} = {num2str(gridHole{chan_ind})};
        dataArray{8} = depth{chan_ind};
        dataArray{9} = psthRunIndex(subRun_ind);
                
        % Concatonate to the correct Matrix (index matching phyzzy convention)
        if unit_ind == 1 % Unsorted
          groupInd = 1;
        elseif unit_ind == length(tmpRunStruct.psthByImage{chan_ind}) % MUA
          groupInd = 3;
        else % Unit
          groupInd = 2;
        end
        for data_ind = 1:length(dataArray)
          stimPSTH{stim_ind,groupInd,data_ind} = [stimPSTH{stim_ind,groupInd,data_ind}; dataArray{data_ind}];
        end
      end
    end
  end
end

% Step 2 - Plots
allStimuliNames = cellfun(@(x) extractBetween(x, 1, length(x)-4), allStimuliVec);

stimPresMat = cellfun(@(x) size(x,1),stimPSTH(:,:,end));
meanPSTHStruct.stimPresMat = stimPresMat;
broadLabelInd = logical(plotMat);
dataInd2Plot = 1:2;

% Plot 1 - All Stimuli means in the same plot.
if params.allStimPSTH
  for group_ind = 1:size(stimPSTH,2)
    stimCounts = stimPresMat(:,group_ind);
    if params.traceCountLabel
      allStimuliLabel = cell(length(allStimuliNames),1);
      for stim_ind = 1:length(allStimuliNames)
        allStimuliLabel{stim_ind} = [allStimuliNames{stim_ind} ',n = ' num2str(stimCounts(stim_ind))];
      end
    else
      allStimuliLabel = allStimuliNames;
    end
    groupData = stimPSTH(:, group_ind, 1);
    groupData = cellfun(@(x) mean(x, 1), groupData, 'UniformOutput',0);
    plotData = vertcat(groupData{:});
    
    if params.sortPresCountSort
      [~, newOrder] = sort(stimCounts);
      allStimuliLabel = allStimuliLabel(newOrder);
      plotData = plotData(newOrder, :);
    end
    
    catPSTHTitle = sprintf('All Stimuli, Mean PSTH %s, %s', normTag, groupingType{group_ind});
    h = figure('NumberTitle', 'off', 'Name', catPSTHTitle,'units','normalized','outerposition',[0 0 params.plotSizeAllStimPSTH]);
    [catPSTHAxes, cbh] = plotPSTH(plotData, [], axes(), params, 'color', catPSTHTitle, allStimuliLabel);
    if params.normalize == 1
      cbh.Label.String = 'Signal Change relative to Baseline (%)';
    elseif params.normalize == 2
      cbh.Label.String = 'Z scored relative to fixation';
      ylabel('Normalized Activity (Baseline Z scored)');
    else
      ylabel('Activity (Firing Rate)');
    end
    title(catPSTHTitle);
    xlabel('Time (ms)');
    saveFigure(params.outputDir, ['1.' catPSTHTitle], [], 1, exportFig, 0, {''})
    close(h)
    clear allStimuliLabel
  end
end

% Plot 2 - Catagory Plot - 'All Chasing Stimuli, mean PSTH'
if params.catPSTH
  for broad_ind = 1:length(params.plotLabels)
    % Extract relevant slices of larger matricies
    sliceStimPSTH = stimPSTH(broadLabelInd(:,broad_ind),:,dataInd2Plot);
    sliceStimLabels = briefStimList(broadLabelInd(:,broad_ind));
    sliceStimPresMat = stimPresMat(broadLabelInd(:,broad_ind),:);
    
    % Generate a matrix of the meanPSTHes, with the last row being total
    % means. Add this to the counts as well.
    meanPSTH = cellfun(@(x) mean(x),sliceStimPSTH,'UniformOutput',0);
    catMeanMat = cell(1,size(sliceStimPSTH,2),size(sliceStimPSTH,3));
    catErrMat = cell(1,size(sliceStimPSTH,2),size(sliceStimPSTH,3));
    for group_ind = 1:length(groupingType)
      for data_ind = dataInd2Plot
        tmp = vertcat(sliceStimPSTH{:,group_ind,data_ind});
        catMeanMat{1,group_ind,data_ind} = mean(tmp,1);
        catErrMat{1,group_ind,data_ind} = std(tmp)/sqrt(size(tmp,1));
      end
    end
    meanPSTH = cat(1, meanPSTH, catMeanMat);
    
    if params.includeMeanTrace
      sliceStimPresMat = [sliceStimPresMat; sum(sliceStimPresMat, 1)];
      if iscell(params.plotLabels{broad_ind})
        params.plotLabels{broad_ind} = strjoin(params.plotLabels{broad_ind});
      end
      sliceStimLabels = [sliceStimLabels; params.plotLabels{broad_ind}];
    end
    
    for group_ind = 1:length(groupingType)
      % Prepare labels
      if params.traceCountLabel
        stimLabels = cell(length(sliceStimLabels),1);
        for stim_ind = 1:length(stimLabels)
          stimLabels{stim_ind} = [sliceStimLabels{stim_ind} ',n = ' num2str(sliceStimPresMat(stim_ind,group_ind))];
        end
      else
        stimLabels = sliceStimLabels;
      end
      
      % Extract correct data
      plotData = vertcat(meanPSTH{:,group_ind,1});
      
      % If Sorting data, do so here.
      if params.sortPresCountSort
        [~, newOrder] = sort(sliceStimPresMat(:,group_ind));
        plotLabels = stimLabels(newOrder);
        plotData = plotData(newOrder, :);
      else
        plotLabels = stimLabels;
      end
      
      psthTitle = sprintf('All %s Stimuli - Mean %s%s, %s', params.plotLabels{broad_ind}, dataType{1}, normTag, groupingType{group_ind});
      % Plot the Activity
      h = figure('NumberTitle', 'off', 'Name', psthTitle,'units','normalized','outerposition',[0 0 params.plotSizeCatPSTH]);
      [psthAxes, cbHandle] = plotPSTH(plotData, [], axes(), params, 'color', psthTitle, plotLabels);
      psthAxes.FontSize = 15;
      psthAxes.YLabel.String = 'Stimulus Name';
      title(psthTitle)
      hold on
      if params.normalize == 1
        cbHandle.Label.String = 'Signal Change relative to Baseline (%)';
      elseif params.normalize == 2
        cbHandle.Label.String = 'Z scored relative to fixation';
      end
      cbHandle.Label.FontSize = 12;
      
      if params.plotMeanLine
        % Include a black line which plots the mean trace on a second axis.
        lineprops.col{1} = 'k';
        mseb(1:size(plotData,2), plotData(end,:), catErrMat{end,group_ind,2}, lineprops,1);
        grandMeanAxes = axes('YAxisLocation','right','Color', 'none','xtick',[],'xticklabel',[],'xlim',[0 size(plotData,2)]);%,'ytick',[],'yticklabel',[]);
        grandMeanAxes.YAxis.FontSize = 12;
        linkprop([psthAxes, grandMeanAxes],{'Position'});
        tmp = psthAxes.Position;
        cbHandle.Position(1) = cbHandle.Position(1) + .03;
        psthAxes.Position = tmp;
        
      end
      
      saveFigure(params.outputDir, ['2.' psthTitle], [], 1, exportFig, 0, {''})
      close(h);
    end
    
  end
end

% Plot 3 - Stimuli Plot - 'All chasing 1 PSTHs, sorted by...'
if params.allRunStimPSTH
  % Make a PSTH of each stimulus across all its repetitions.
  sortType = {'Run of Day', 'Grid Hole', 'Recording Depth', 'Run Ind'};
  % {'Run of Day', 'Grid Hole', 'Recording Depth', 'Run Ind'};
  % Must be the same order as the final parts of stimPSTH
  sortMat = cell(size(stimPSTH,1), size(stimPSTH,2), length(sortType));
  sortLabel = cell(size(stimPSTH,1), size(stimPSTH,2), length(sortType));
  
  % Generate the PSTH sorting indicies
  for stim_i = 1:size(stimPSTH,1)
    for group_i = 1:size(stimPSTH,2)
      for sort_i = 1:length(sortType) % the sortings that need processing
        
        if sort_i == 1 || sort_i == 2
          % Store indicies in 1st, labels in 2nd %Grid Holes --> Indicies
          [tmpLabels , sortMat{stim_i, group_i, sort_i}] = sort(stimPSTH{stim_i, group_i,strcmp(dataType, sortType{sort_i})});
          
          % Generate Labeling index
          uniqueLabels = unique(tmpLabels);
          labeledIndex = zeros(length(uniqueLabels),1);
          for uni_i = 1:length(uniqueLabels)
            if sort_i == 1
              labeledIndex(uni_i) = find(tmpLabels == uniqueLabels(uni_i), 1);
            elseif sort_i == 2
              labeledIndex(uni_i) = find(strcmp(tmpLabels, uniqueLabels(uni_i)),1);
            end
          end
          
        elseif sort_i == 3 %Recording Depth --> Indicies
          tmpDepths = [stimPSTH{stim_i, group_i,strcmp(dataType, sortType{sort_i})}];
          [tmpLabels , sortMat{stim_i, group_i, sort_i}] = sort(tmpDepths');
          labeledIndex = 1:5:length(tmpLabels);
          
        elseif sort_i == 4
          sortMat{stim_i, group_i, sort_i} = 1:length(stimPSTH{stim_i, group_i, strcmp(dataType, sortType{sort_i})});
          tmpLabels = runList(stimPSTH{stim_i, group_i,strcmp(dataType, sortType{sort_i})});
          labeledIndex = 1:5:length(tmpLabels);
        end
        
        % Use the labeledIndex to generate the proper label array with 0's
        % everywhere else.
        sortLabelTmp = cell(length(tmpLabels),1);
        for tmp_i = 1:length(labeledIndex)
          if sort_i == 2 || sort_i == 4
            sortLabelTmp{labeledIndex(tmp_i)} = tmpLabels{labeledIndex(tmp_i)};
          else
            sortLabelTmp{labeledIndex(tmp_i)} = tmpLabels(labeledIndex(tmp_i));
          end
        end
        sortLabel{stim_i, group_i, sort_i} = tmpLabels;
        
      end
    end
  end
  
  % Iterate through PSTH, generating plots
  for stim_i = 1:length(stimPSTH)
    for sort_i = 4%1:length(sortType)
      for group_i = 2:size(stimPSTH,2)
        figTitle = sprintf('%s - %s PSTHs, Sorted by %s, %s', allStimuliNames{stim_i}, groupingType{group_i} ,sortType{sort_i}, normTag);
        h = figure('NumberTitle', 'off', 'Name', figTitle,'units','normalized','outerposition',[0 0 params.plotSizeAllRunStimPSTH]);
        sgtitle(figTitle)
        
        sortIndex = sortMat{stim_i, group_i, sort_i};
        
        % Break up PSTHes with many lines into subplots.
        TracesPerPlot = ceil(size(sortIndex,2)/3);
        if TracesPerPlot < 20
          subplot2Plot = 1;
        elseif TracesPerPlot < 45
          subplot2Plot = 2;
        else
          subplot2Plot = 3;
        end
        TracesPerPlot = ceil(size(sortIndex,2)/subplot2Plot);
        
        plotStarts = 1:TracesPerPlot:size(sortIndex,2);
        plotEnds = [plotStarts(2:end)-1, size(sortIndex,2)];
        
        subplotAxes = gobjects(subplot2Plot,1);
        cbHandle = gobjects(subplot2Plot,1);
        for plot_i = 1:subplot2Plot
          plotLabels = sortLabel{stim_i, group_i, sort_i}(plotStarts(plot_i):plotEnds(plot_i));
          plotData = stimPSTH{stim_i,group_i,1}(plotStarts(plot_i):plotEnds(plot_i), :);
          %sortIndex = sortMat{stim_i, group_i, sort_i}(plotStarts(plot_i):plotEnds(plot_i));
          psthAxes = subplot(1,subplot2Plot,plot_i);
          [subplotAxes(plot_i), cbHandle(plot_i)] = plotPSTH(plotData, [], psthAxes, params, 'color', [], plotLabels); %(sortIndex,:)
          cbHandle(plot_i).Label.FontSize = 12;
          if plot_i == subplot2Plot
            if params.normalize == 1
              cbHandle(plot_i).Label.String = 'Signal Change relative to Baseline (%)';
            elseif params.normalize == 2
              cbHandle(plot_i).Label.String = 'Z score relative to fixation';
            end
          else
            delete(cbHandle(plot_i))
          end
          set(gca,'FontSize',10,'TickLength',[.01 .01],'LineWidth',.25);
        end
        linkprop(subplotAxes, 'CLim');
        
        saveFigure(params.outputDir, ['3.' figTitle], [], 1, exportFig, 0, {''})
        close(h)
      end
    end
  end
  
end

% Plot 4 - Line plot with Line per Catagory
if params.lineCatPlot
  % Find the end of the 'stimuli catagory' count.
  catCount = find(strcmp(params.plotLabels, 'scene'));
  if isempty(catCount)
    catCount = find(strcmp(params.plotLabels, 'grooming'));
  end
  singleCatPlotMat = plotMat(:,1:catCount);
  singleCatplotLabels = params.plotLabels(1:catCount);
  singleCatplotLabelsSocialInd = params.plotLabelSocialInd(1:catCount);
  %tmp = distinguishable_colors(catCount);
  plotColors = cell(catCount,1);
  socialColorsHSV = repmat(rgb2hsv(params.socialColor), [sum(singleCatplotLabelsSocialInd),1]);
  nonSocialColorsHSV = repmat(rgb2hsv(params.nonSocialColor), [sum(singleCatplotLabelsSocialInd),1]);
  % HSV values cap at 240
  hsvGradient = 25/240;
  tmp = 1:5;
  hsvGradient = hsvGradient * tmp';
  
  socialColorsHSV(:,3) = socialColorsHSV(:,3) - hsvGradient;
  nonSocialColorsHSV(:,3) = nonSocialColorsHSV(:,3) - hsvGradient;
  plotColorTmp = [socialColorsHSV; nonSocialColorsHSV];
  plotColorTmp = hsv2rgb(plotColorTmp);
  
  for col_i = 1:catCount
    plotColors{col_i} = plotColorTmp(col_i,:);
  end
  
  for group_ind = 1:size(stimPSTH,2)
    catPSTHTitle = sprintf('%s Mean PSTH %s', groupingType{group_ind}, normTag);
    h = figure('NumberTitle', 'off', 'Name', catPSTHTitle,'units','normalized','outerposition',[0 0 params.plotSizeLineCatPlot]);
    hold on
    h.Children.FontSize = 15;
    % Generate line plots w/ error bars.
    lineProps.width = 2;
    lineProps.col = plotColors;
    lineProps.patch.FaceAlpha = '0.5';
    
    [groupData, groupErr] = deal(cell(size(singleCatPlotMat,2),1));
    for cat_ind = 1:size(singleCatPlotMat,2)
      tmp = vertcat(stimPSTH{logical(singleCatPlotMat(:,cat_ind)), group_ind});
      groupData{cat_ind} = mean(tmp,1);
      groupErr{cat_ind} = std(tmp)/sqrt(size(tmp,1));
    end
    
    groupData = vertcat(groupData{:});
    
    if params.fixAlign
      groupMean = mean(mean(groupData(:,500:800)));
      for g_ind = 1:size(groupData,1)
        groupData(g_ind,:) = groupData(g_ind,:) - mean(groupData(g_ind,500:800)) + groupMean;
      end
    end
    
    mseb(-800:(size(groupData,2)-801),groupData, vertcat(groupErr{:}), lineProps);
    xlim([-800 (size(groupData,2)-801)])
    legend(singleCatplotLabels, 'AutoUpdate','off','location', 'northeastoutside');
    line([0 0], ylim(), 'Linewidth',3,'color','k');
    line([2800 2800], ylim(), 'Linewidth',3,'color','k');
    xlabel('Time (ms)');
    ylabel('Normalized Activity (Baseline Z scored)');
    title(catPSTHTitle);
    saveFigure(params.outputDir, ['4.' catPSTHTitle], [], 1, exportFig, 0, {''})
    close(h)
  end
  
end

% Plot 5 - Means across broad catagorizations (like Social vs non Social)
if params.lineBroadCatPlot
  
  % Generate line plots across different labeling schemes.
  socialInd = logical(plotMat(:,strcmp(params.plotLabels,'socialInteraction')));
  agentInd = logical(plotMat(:,strcmp(params.plotLabels,'agents')));
  headTurnInd = logical(plotMat(:,strcmp(params.plotLabels,'headTurning')));
  
  % Agent videos without head turning.
  agentNHInd = agentInd;
  agentNHInd(headTurnInd) = false;
  
  % Social videos, Head turning vs Not
  htSocialIndTmp = socialInd + headTurnInd;
  socialHTInd = false(length(htSocialIndTmp), 1);
  socialHTInd(htSocialIndTmp == 2) = true;
  socialNonHTInd = socialInd;
  socialNonHTInd(headTurnInd) = false;
  
  % Social Agent vs Non-Social Agent
  socialAgentInd = socialInd;
  nonSocialAgentInd = logical(agentInd-socialInd);
  
  figIncInd = {agentInd, socialInd, headTurnInd, socialHTInd, socialAgentInd};
  figExcInd = {~agentInd, ~socialInd, agentNHInd, socialNonHTInd, nonSocialAgentInd};
  figTitInd = {'Agent containing Stimuli Contrast',...
    'Social Interactions Contrast',...
    'Agents engaging in Head Turning Contrast',...
    'Social Interactions with Head turning Contrast',...
    'Agents engaging in Social Interactions Contrast'};
  figLegends = {{'Agent containing Stimuli','Non-Agent containing Stimuli'},...
    {'Social Interaction Stimuli','non-Social Interaction Stimuli'},...
    {'Agents with Head Turning','Agents without Head Turning'},...
    {'Socially Interacting Agents with Head Turning','Socially Interacting Agents without  Head Turning'},...
    {'Agents engaging in Social Interactions ','Agents not engaging in Social Interactions'}};
  lineprops.col = {params.socialColor, params.nonSocialColor};
  
  for fig_ind = 1:length(figIncInd)
    line1Array = stimPSTH(figIncInd{fig_ind},:,1);
    line2Array = stimPSTH(figExcInd{fig_ind},:,1);
    if ~isempty(line1Array) && ~isempty(line2Array)
      for group_ind = 2:size(line1Array,2)
        % Generate Data
        line1Data = vertcat(line1Array{:,group_ind});
        line2Data = vertcat(line2Array{:,group_ind});
        %       sigBins = zeros(size(line1Data,2),1);
        %       for ii = 1:size(line1Data, 2)
        %         sigBins(ii) = ttest2(line1Data(:, ii), line2Data(:, ii));
        %       end
        %       sigBins = find(sigBins);
        lineMean = [mean(line1Data, 1); mean(line2Data, 1)];
        lineErr = [std(line1Data)/sqrt(size(line1Data,1)); std(line2Data)/sqrt(size(line2Data,1))];
        
        if params.fixAlign
          fixMean = mean(mean(lineMean(:,500:800),1));
          for line_i = 1:size(lineMean,1)
            lineMean(line_i,:) = lineMean(line_i,:) - mean(lineMean(line_i,500:800)) + fixMean;
          end
        end
        
        % Prepare figure
        plotTitle = sprintf('%s - %s, %s', ['mean PSTH of ' figTitInd{fig_ind}], groupingType{group_ind}, normTag);
        h = figure('NumberTitle', 'off', 'Name', plotTitle,'units','normalized','outerposition',[0 0 params.plotSizeLineBroadCatPlot]);
        hold on
        h.Children.FontSize = 15;
        title(plotTitle)
        mseb((-800:size(lineMean,2)-801),lineMean, lineErr, lineprops, []);
        legend(figLegends{fig_ind}, 'AutoUpdate','off')
        xlim([-800 (length(lineMean)-801)])
        ylim(ylim());
        line([0 0], ylim(), 'Linewidth',3,'color','k');
        line([2800 2800], ylim(), 'Linewidth',3,'color','k');
        xlabel('Time (ms)')
        ylabel('Normalized Activity')
        
        %       lineMax = max(lineMean,[],1);
        %       for star_i = 1:length(sigBins)
        %         text(sigBins(star_i)-799,(lineMax(sigBins(star_i))*1.05),'*','FontSize',20)
        %       end
        
        saveFigure(params.outputDir, ['5.' plotTitle], [], 1, exportFig, 0, {''})
        close(h)
      end
    end
  end
  
end

end

function spikeDataBank = slidingWindowTest(spikeDataBank,params)
disp('Starting sliding window test...');
% Perform sliding scale ANOVA, calculate Omega at each bin.

if ~exist(params.outputDir, 'dir')
  mkdir(params.outputDir)
end

% Generate variables needed
runList = fields(spikeDataBank);
binSize = params.binSize;
binStep = params.binStep;
target = params.target;
Omega = params.Omega;
stimParamFile = params.stimParamFile;
nestStrCmp = @(x, y) any(strcmp(x, y));

% Step 1 - for each bin (epoch), calculate the rates and counts. save them
% into a larger epochRates structure, which can be stored in the folder for
% this function for subsequent runs. 

if ~exist(fullfile(params.outputDir, 'epochRates.mat'), 'file')
  epochRates = cell(length(runList),1);
  catagoryInd = cell(length(runList),1);
  % Step 1 - generate bin times and spike rates, and proper memberships to groups.
  for run_ind = 1:length(runList)
    runStruct = spikeDataBank.(runList{run_ind});
    starts = (runStruct.start:binStep:(runStruct.end - binSize))';
    ends = (runStruct.start+binSize:binStep:(runStruct.end))';
    spikeDataBank.(runList{run_ind}).epochs = [starts,ends];
    epochRates{run_ind} = cell(length(starts),1);
    for bin_ind = 1:length(starts)
      [epochRates{run_ind}{bin_ind}, ~, ~] = spikeCounter(spikeDataBank.(runList{run_ind}).spikesByEvent, starts(bin_ind), ends(bin_ind));
    end
    for group_ind = 1:length(target)
      catagoryInd{run_ind}(:,group_ind) = cell2mat(cellfun(nestStrCmp, runStruct.eventCategories, repmat(target(group_ind), [length(runStruct.eventCategories),1]), 'UniformOutput',0));
    end
  end
  save(fullfile(params.outputDir, 'epochRates'), 'epochRates', 'catagoryInd');
else
  load(fullfile(params.outputDir, 'epochRates'), 'epochRates', 'catagoryInd');
  assert(length(epochRates) == length(runList), 'epochRates length does not match runList length');
  fprintf('Sliding windows rates already calculated, continuing... \n');
end

% Step 2 - Perform ANOVA, Omega calculation across bins, and store values.
% Save the output in the function folder for later retrieval. 

statsType = {'groupRates', 'nonGroupRates', 'pVec', 'cohensVec', 'nullPVec', 'nullCohensVec', 'sigBins'};

if ~exist(fullfile(params.outputDir, 'allStatsArray.mat'), 'file')
  chanUnitStruct = structfun(@(x) x.psthByImage, spikeDataBank,'UniformOutput',0);
  chanUnitStruct = struct2cell(chanUnitStruct);
  allStatsArray = initNestedCellArray(chanUnitStruct, 'cell', [length(epochRates{1}), length(target), length(statsType)], 3);
  clear chanUnitStruct;
  scrambleCount = params.scrambleCount;
  
  tic
  if license('test','Distrib_Computing_Toolbox')
    
    epochRatesPar = parallel.pool.Constant(epochRates);
    catagoryIndPar =  parallel.pool.Constant(catagoryInd);
%     epochRatesPar.Value = epochRates; for debugging.
%     catagoryIndPar.Value =  catagoryInd;
    parfor run_ind = 1:length(runList)
      % Initialize relevant Structures
      if length(unique(catagoryIndPar.Value{run_ind})) > 1 %Only run tests where you have members in each group.
        binCount = length(epochRatesPar.Value{run_ind});
        for chan_ind = 1:length(epochRatesPar.Value{run_ind}{1})
          for unit_ind = 1:length(epochRatesPar.Value{run_ind}{1}{chan_ind})
            for bin_ind = 1:binCount
              unitResponsePerEvent = epochRatesPar.Value{run_ind}{bin_ind}{chan_ind}{unit_ind};
              catagortyIndSlice = catagoryIndPar.Value{run_ind};
              %unitData{event}.rates = trial*1
              for target_ind = 1:length(target)
                if length(unique(catagortyIndSlice(:,target_ind))) == 1
                  break % If there is only 1 label type, there is no comparison to be made.
                end
                [trialSpikes, trialLabels]  = deal([]);
                % grab the relevant events
                targetInd = catagortyIndSlice(:,target_ind);
                targetSpikes = unitResponsePerEvent(targetInd);
                otherSpikes = unitResponsePerEvent(~targetInd);
                % Initialize relevant vecotrs
                spikeGroups = {targetSpikes otherSpikes};
                spikeGroupLabels ={(target{target_ind}) (['non-' target{target_ind}])};
                % Cluster and reshape the arrays properly
                for group_i = 1:length(spikeGroups)
                  tmp = spikeGroups{group_i};
                  tmp = [tmp{:}];
                  dataVec = vertcat(tmp.rates);
                  labelVec = repmat(spikeGroupLabels(group_i), length(dataVec),1);
                  trialSpikes = vertcat(trialSpikes,dataVec);
                  trialLabels = vertcat(trialLabels, labelVec);
                end
                % Check for social v non-social
                % 3rd ind in 3D mat statsType = {'pVec', 'cohensVec', 'nullPVec', 'nullCohensVec'};
                pop1 = trialSpikes(strcmp(trialLabels,spikeGroupLabels{1}));
                pop2 = trialSpikes(strcmp(trialLabels,spikeGroupLabels{2}));
                allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('groupRates',statsType)} = pop1;
                allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('nonGroupRates',statsType)} = pop2;
                %[allStatsArray{run_ind}{chan_ind}{unit_ind}(bin_ind,target_ind, strcmp('pVec',statsType)), pStatsTable, ~] = anovan(trialSpikes,{trialLabels},'model','interaction','varnames',{'SvNS'}, 'alpha', 0.05,'display','off');
                [allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('sigBins',statsType)}, allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('pVec',statsType)}, ~, pStatsTable] = ttest2(pop1, pop2); %,'model','interaction','varnames',{'SvNS'}, 'alpha', 0.05,'display','off');
                allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('cohensVec',statsType)} = (mean(pop1) - mean(pop2))/pStatsTable.sd;
                [nullPVec, nullCohensVec] = deal(zeros(1,scrambleCount));
                for rand_ind = 1:scrambleCount
                  trialLabelsTmp = trialLabels(randperm(length(trialLabels)));
                  pop1 = trialSpikes(strcmp(trialLabelsTmp,spikeGroupLabels{1}));
                  pop2 = trialSpikes(strcmp(trialLabelsTmp,spikeGroupLabels{2}));
                  [~, nullPVec(rand_ind), ~, pStatsTableTmp] = ttest2(pop1, pop2);
                  nullCohensVec(rand_ind) = (mean(pop1) - mean(pop2))/pStatsTableTmp.sd;
                end
                allStatsArray{run_ind}{chan_ind}{unit_ind}(bin_ind,target_ind, strcmp('nullPVec',statsType)) = {nullPVec};
                allStatsArray{run_ind}{chan_ind}{unit_ind}(bin_ind,target_ind, strcmp('nullCohensVec',statsType)) = {nullCohensVec};
              end
            end
          end
        end
      else
        disp('Skipping')
      end
    end
    
  else
    
    for run_ind = 1:length(runList)
      % Initialize relevant Structures
      if length(unique(catagoryInd{run_ind})) > 1 %Only run tests where you have members in each group.
        chanCount = length(epochRates{run_ind}{1});
        binCount = length(epochRates{run_ind});
        for chan_ind = 1:chanCount
          for unit_ind = 1:length(epochRates{run_ind}{1}{chan_ind})
            for bin_ind = 1:binCount
              unitResponsePerEvent = epochRates{run_ind}{bin_ind}{chan_ind}{unit_ind};
              catagortyIndSlice = catagoryInd{run_ind};
              %unitData{event}.rates = trial*1
              for target_ind = 1:length(target)
                if length(unique(catagortyIndSlice(:,target_ind))) == 1
                  break % If there is only 1 label type, there is no comparison to be made.
                end
                [trialSpikes, trialLabels]  = deal([]);
                % grab the relevant events
                targetInd = catagortyIndSlice(:,target_ind);
                targetSpikes = unitResponsePerEvent(targetInd);
                otherSpikes = unitResponsePerEvent(~targetInd);
                % Initialize relevant vecotrs
                spikeGroups = {targetSpikes otherSpikes};
                spikeGroupLabels ={(target{target_ind}) (['non-' target{target_ind}])};
                % Cluster and reshape the arrays properly
                for group_i = 1:length(spikeGroups)
                  tmp = spikeGroups{group_i};
                  tmp = [tmp{:}];
                  dataVec = vertcat(tmp.rates);
                  labelVec = repmat(spikeGroupLabels(group_i), length(dataVec),1);
                  trialSpikes = vertcat(trialSpikes,dataVec);
                  trialLabels = vertcat(trialLabels, labelVec);
                end
                % Check for social v non-social
                % 3rd ind in 3D mat statsType = {'pVec', 'cohensVec', 'nullPVec', 'nullCohensVec'};
                pop1 = trialSpikes(strcmp(trialLabels,spikeGroupLabels{1}));
                pop2 = trialSpikes(strcmp(trialLabels,spikeGroupLabels{2}));
                allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('groupRates',statsType)} = pop1;
                allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('nonGroupRates',statsType)} = pop2;
                %[allStatsArray{run_ind}{chan_ind}{unit_ind}(bin_ind,target_ind, strcmp('pVec',statsType)), pStatsTable, ~] = anovan(trialSpikes,{trialLabels},'model','interaction','varnames',{'SvNS'}, 'alpha', 0.05,'display','off');
                [allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('sigBins',statsType)}, allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('pVec',statsType)}, ~, pStatsTable] = ttest2(pop1, pop2); %,'model','interaction','varnames',{'SvNS'}, 'alpha', 0.05,'display','off');
                allStatsArray{run_ind}{chan_ind}{unit_ind}{bin_ind,target_ind, strcmp('cohensVec',statsType)} = (mean(pop1) - mean(pop2))/pStatsTable.sd;
                [nullPVec, nullCohensVec] = deal(zeros(1,scrambleCount));
                for rand_ind = 1:scrambleCount
                  trialLabelsTmp = trialLabels(randperm(length(trialLabels)));
                  pop1 = trialSpikes(strcmp(trialLabelsTmp,spikeGroupLabels{1}));
                  pop2 = trialSpikes(strcmp(trialLabelsTmp,spikeGroupLabels{2}));
                  [~, nullPVec(rand_ind), ~, pStatsTableTmp] = ttest2(pop1, pop2);
                  nullCohensVec(rand_ind) = (mean(pop1) - mean(pop2))/pStatsTableTmp.sd;
                end
                allStatsArray{run_ind}{chan_ind}{unit_ind}(bin_ind,target_ind, strcmp('nullPVec',statsType)) = {nullPVec};
                allStatsArray{run_ind}{chan_ind}{unit_ind}(bin_ind,target_ind, strcmp('nullCohensVec',statsType)) = {nullCohensVec};
              end
            end
          end
        end
      else
        disp('Skipping')
      end
    end
    
  end
  save(fullfile(params.outputDir, 'allStatsArray'), 'allStatsArray');
  fprintf('Done in %d hours \n', toc/3600)
else
  fprintf('Loading statsArray from saved data... \n')
  load(fullfile(params.outputDir, 'allStatsArray.mat'), 'allStatsArray');
  assert(length(allStatsArray) == length(runList), 'allStatsArray length does not match runList length');
end

runStruct = spikeDataBank.(runList{1});
starts = (runStruct.start:binStep:(runStruct.end - binSize))';
xlabelPlot = starts(1:8:length(starts));
stimOnInd = find(starts == 0);
stimEndInd = find(starts == 2800);

if params.plotTest
  for run_i = 1:length(allStatsArray)
    statsArray = allStatsArray{run_i};
    % Plot the results for each unit seen
    for chan_ind = 1:length(statsArray)
      for unit_ind = 1:length(statsArray{chan_ind})
        if unit_ind == 1
          ANOVAvarName = ['Ch' num2str(chan_ind) ' Unsorted - SocVsNonSoc'];
        elseif unit_ind == length(statsArray{chan_ind})
          ANOVAvarName = ['Ch' num2str(chan_ind) ' MUA - SocVsNonSoc'];
        else
          ANOVAvarName = ['Ch' num2str(chan_ind) ' U' num2str(unit_ind-1) ' - SocVsNonSoc'];
        end
        % {'groupRates', 'nonGroupRates', 'pVec', 'cohensVec', 'nullPVec', 'nullCohensVec'};
        titleString = sprintf("%s - Sliding Window Group mean comparison - %s", runList{run_i} ,ANOVAvarName);
        fileString = fullfile(spikeDataBank.(runList{run_i}).figDir, [titleString{1}]);
        if ~exist([fileString '.fig'], 'file') || (~exist([fileString '.png'], 'file') && params.exportFig)
          fullfile(spikeDataBank.(runList{run_i}).figDir, titleString)
          groupMeans = cellfun(@(x) mean(x), statsArray{chan_ind}{unit_ind}(:,:,strcmp('groupRates',statsType)));
          if ~(nnz(groupMeans) == 0)
            groupSEMs = cellfun(@(x) (std(x)/sqrt(length(x))), statsArray{chan_ind}{unit_ind}(:,:,strcmp('groupRates',statsType)));
            nonGroupMeans = cellfun(@(x) mean(x), statsArray{chan_ind}{unit_ind}(:,:,strcmp('nonGroupRates',statsType)));
            nonGroupSEMs = cellfun(@(x) (std(x)/sqrt(length(x))), statsArray{chan_ind}{unit_ind}(:,:,strcmp('nonGroupRates',statsType)));
            pVec = cell2mat(statsArray{chan_ind}{unit_ind}(:,:,strcmp('pVec',statsType)));
            %cohensVec = cell2mat(statsArray{chan_ind}{unit_ind}(:,:,strcmp('cohensVec',statsType)));
            if size(pVec,2) ~= length(target)
              tmp = statsArray{chan_ind}{unit_ind}(:,:,strcmp('pVec',statsType));
              emptyInd = cellfun(@(x) isempty(x), tmp);
              tmp(emptyInd) = deal({1});
              pVec = cell2mat(tmp);
            end
            % plot
            h = figure('NumberTitle', 'off', 'Name', titleString,'units','normalized','outerposition',[0 0 params.plotSize]);
            sgtitle(titleString)
            for target_i = 1:length(target)
              subplot(length(target),1,target_i)
              title(target{target_i});
              lineData = [groupMeans(:, target_i), nonGroupMeans(:, target_i)];
              errorData = [groupSEMs(:, target_i), nonGroupSEMs(:, target_i)];
              mseb(1:length(lineData), lineData', errorData');
              legend([target(target_i), ['non-' target{target_i}]], 'location', 'northwest','AutoUpdate','off')
              hold on
              xlabel('Bin Start time (ms)')
              ylabel('Firing Rates (Hz)')
              xlim([0,length(lineData)]);
              xticks(1:8:length(starts))
              xticklabels(xlabelPlot);
              plot([stimOnInd, stimOnInd], ylim(),'k','linewidth',4);
              plot([stimEndInd, stimEndInd], ylim(),'k','linewidth',4);
              
              % Star Significant Bins
              sigBins = find(pVec(:, target_i) < 0.05);
              lineMax = max(lineData,[],2);
              for star_i = 1:length(sigBins)
                text(sigBins(star_i),(lineMax(sigBins(star_i))*1.05),'*','FontSize',20)
              end
            end
            % Save and close
            % (outDir, filename, figData, saveFig, exportFig, saveData, varargin )
            saveFigure(spikeDataBank.(runList{run_i}).figDir, [titleString], 1, 1, params.exportFig, 0, {''});
            close(h)
          end
        end
      end
    end
  end
end

fprintf('Sliding Window Calculations finished... \n')

% Step 3 - Count Results of ANOVA across all units and count stretches of significant
% bins.
totalUnitCount = 0;
totalChannelCount = 0;
groupingType = {'Unsorted', 'Unit', 'MUA'};
% target = {'socialInteraction','agents','interaction'};, pulled from params.
dataType = {'binTotal','sigBins','nullSigBins'};
statType = {'Counts', 'Sig Run Lengths', 'Mean Stretch Lengths', 'Max run lengths', 'Sig Run Starts'};
cohensStatType = {'Bin Value', 'Maximum values'};

% Find Bins in each Epoch of the stimulus.
starts = (runStruct.start:binStep:(runStruct.end - binSize))';
epochNames = {'Early','Middle','Late'};
epochBinTimes = [[find(starts == 50) find(starts == 1025)];...
                [find(starts == 1050) find(starts == 2025)];...
                [find(starts == 2050) find(starts == 2850)]];
epochBinSig = zeros(length(epochNames), length(groupingType) , length(target));

testMat = zeros(length(epochRates{1}), length(target), length(groupingType), length(dataType)); % (binCount, target, group, dataType)
targetRunLens = cell(length(target), length(groupingType), length(statType));
unitSigEpochData = {'Count', 'Run Ch Unit Index'};
unitSigEpoch = cell(length(epochNames), length(groupingType) , length(target), length(unitSigEpochData));
cohensData = cell(length(target), length(groupingType), length(cohensStatType));

for run_ind = 1:length(allStatsArray)
  for chan_ind = 1:length(allStatsArray{run_ind})
    unitCount = length(allStatsArray{run_ind}{chan_ind});
    totalUnitCount = totalUnitCount + (unitCount - 2);
    totalChannelCount = totalChannelCount + 1;
    for unit_ind = 1:unitCount
      % {'groupRates', 'nonGroupRates', 'pVec', 'cohensVec', 'nullPVec', 'nullCohensVec'};
      uPVec = cell2mat(allStatsArray{run_ind}{chan_ind}{unit_ind}(:,:,strcmp('pVec',statsType)));
      if size(uPVec,2) ~= length(target)
        tmp = allStatsArray{run_ind}{chan_ind}{unit_ind}(:,:,strcmp('pVec',statsType));
        emptyInd = cellfun(@(x) isempty(x), tmp);
        tmp(emptyInd) = deal({0});
        uPVec = cell2mat(tmp);
      end
      
      if ~isempty(uPVec)
        % Check for significant bins
        uPVec(uPVec == 0) = nan;
        tmpSigCount = double(uPVec < 0.05);
        tmpNullSigCount = cellfun(@(x) sum(x < 0.05), allStatsArray{run_ind}{chan_ind}{unit_ind}(:,:,strcmp('nullPVec',statsType)));
        cohensVec = allStatsArray{run_ind}{chan_ind}{unit_ind}(:,:,strcmp('cohensVec',statsType));
%         omegaVec(isnan(omegaVec)) = 0;
        % Add to relevant structures
        if unit_ind == 1
          group_i = 1;
        elseif unit_ind == unitCount
          group_i = 3;
        else
          group_i = 2;
        end
        % Save significant bins and null bins correctly.
        %{'binTotal','sigBins','nullSigBins'}
        testMat(:, :, group_i, strcmp(dataType, 'binTotal')) = testMat(:, :, group_i, strcmp(dataType, 'binTotal')) + ~isnan(uPVec);
        testMat(:, :, group_i, strcmp(dataType, 'sigBins')) = testMat(:, :, group_i, strcmp(dataType, 'sigBins')) + tmpSigCount;
        testMat(:, :, group_i, strcmp(dataType, 'nullSigBins')) = testMat(:, :, group_i, strcmp(dataType, 'nullSigBins')) + tmpNullSigCount;
        % Keep count of consecutive bins
        for targ_i = 1:length(target)
          trace = tmpSigCount(:,targ_i);
          trace(isnan(trace)) = 0;
          starts = find(diff([0; trace]) == 1);
          ends = find(diff([trace; 0]) == -1)+1;
          runLengths = ends - starts;
          % {'Counts', 'Sig Run Lengths', 'Sig Run Starts', 'Mean Stretch Lengths', 'Max run lengths'};
          targetRunLens{targ_i, group_i, strcmp(statType, 'Counts')} = [targetRunLens{targ_i, group_i, strcmp(statType, 'Counts')}; sum(trace)];
          targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Lengths')} = [targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Lengths')}; runLengths];
          targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Starts')} = [targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Starts')}; starts];
          targetRunLens{targ_i, group_i, strcmp(statType, 'Mean Stretch Lengths')} = [targetRunLens{targ_i, group_i, strcmp(statType, 'Mean Stretch Lengths')}; mean(runLengths)];
          targetRunLens{targ_i, group_i, strcmp(statType, 'Max run lengths')} = [targetRunLens{targ_i, group_i, strcmp(statType, 'Max run lengths')}; max(runLengths)];
          
          cohensTrace = cell2mat(cohensVec(:,targ_i)');
          cohensData{targ_i, group_i, 1} = [cohensData{targ_i, group_i, 1}; cohensTrace];
          cohensData{targ_i, group_i, 2} = [cohensData{targ_i, group_i, 2}; max(cohensTrace)];
          
          %See if the unit has a significant stretch of bins during each
          %epoch
          if any(runLengths >= 3)
            
            runsOfInterest = runLengths >= 3;
            startsOfInterest = starts(runsOfInterest);
            for epoch_i = 1:length(epochNames)
              % Check where run belongs
              % unitSigEpoch = (epochNames), (groupingType) ,(target));
              startBin = epochBinTimes(epoch_i,1);
              endBin = epochBinTimes(epoch_i,2);
              newSigRuns = ((startsOfInterest >= startBin) + (startsOfInterest <= endBin) == 2);
              if any(newSigRuns)
                %{'Count', 'Run Ch Unit Index'};
                unitSigEpoch{epoch_i, group_i, targ_i, strcmp(unitSigEpochData, 'Count')} = [unitSigEpoch{epoch_i, group_i, targ_i, strcmp(unitSigEpochData, 'Count')};  sum(newSigRuns)];
                unitSigEpoch{epoch_i, group_i, targ_i, strcmp(unitSigEpochData, 'Run Ch Unit Index')} = [unitSigEpoch{epoch_i, group_i, targ_i, strcmp(unitSigEpochData, 'Run Ch Unit Index')};  [run_ind chan_ind unit_ind]];
              end
            end

          end
        end
      end
      
    end
  end
end

% for targ_i = 1:length(target)
%   for group_i = 1:length(groupingType)
%     binLengths = targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Lengths')};
%     binStarts = targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Starts')};
%     longBinInd = binLengths >= 3;
%     binStarts(~longBinInd) = 0; 
%     binLengths(~longBinInd) = 0; 
%     for epoch_i = 1:length(epochNames)
%       startBin = epochBinTimes(epoch_i,1);
%       endBin = epochBinTimes(epoch_i,2);
%       binsOfInterest = ((binStarts >= startBin) + (binStarts <= endBin) == 2);
%       epochBinSig(epoch_i, group_i, targ_i) = sum(binLengths(binsOfInterest));
%     end
%   end
% end

% Go through the unitSigEpoch structure, generate a count of how many
% units, MUA had a significant run of at least 3 bins in the early, middle,
% and late phases of stimulus presentation.
% unitSigEpoch = ('Early,'Mid','Late'), ('Unsort','Unit','MUA') ,(target), ('Count', 'Run Ch Unit Index'));
unitSigEpochSoc = squeeze(unitSigEpoch(:, 2, 1, :));
MUASigEpochSoc = squeeze(unitSigEpoch(:, 3, 1, :));

% totalFractionUnits - how many units have at least 1 bout?
totalUnitswitSigRun = length(unique(vertcat(unitSigEpochSoc{:,2}), 'row'));
totalMUAswitSigRun = length(unique(vertcat(MUASigEpochSoc{:,2}), 'row'));

% Plot Results - using values found in ANOVAmat and targetRunLens
% Fig 1 - 'Social Interactions, Run Length Stats' 4 * 3 grid, 1 subplot per
% group and stat.
statTypePlotLabel = {'Number of Significant Bins', 'Length of Significant Stretches'};
statTypePlotLabel2 = {'Average Significant Stretch Length','Max Significant Stretch Length'};
statTypePlots = {statTypePlotLabel statTypePlotLabel2};
xLabelPlotArray = {'Counts', 'Run Lengths', 'Run Lengths', 'Run Lengths'};

for targ_i = 1:length(target)
  for stat_i = 1:length(statTypePlots)
    statTypePlotLabelTmp = statTypePlots{stat_i};
    figTitle = sprintf('Summary Statistic for Sliding T tests %d - %s', stat_i,  target{targ_i});
    h = figure('NumberTitle', 'off', 'Name', figTitle,'units','normalized','outerposition',[0 0 params.plotSize]);
    sgtitle(figTitle);
    for group_i = 2:length(groupingType)
      for plot_i = 1:length(statTypePlotLabelTmp)
        plotInd = (((group_i - 2) * length(statTypePlotLabelTmp))) + plot_i;
        subplot(length(groupingType)-1, length(statTypePlotLabelTmp), plotInd);
        dataInd = (((stat_i - 1) * length(statTypePlotLabelTmp))) + plot_i;
        plotData = targetRunLens{targ_i, group_i, dataInd};
        histogram(plotData);
        if plot_i == 1
          ylabel(groupingType{group_i})
        end
        xlabel(xLabelPlotArray{dataInd})
        plotTitle = sprintf('%s (90th p = %s)', statTypePlotLabelTmp{plot_i}, num2str(round(prctile(plotData, 90), 2)));
        title(plotTitle);
      end
    end
    saveFigure(params.outputDir, ['1.' figTitle], [], 1, params.exportFig, 0, {''})
    close(h)
  end
end

% Figure 1.5 - The Population Landscape - find runs of 3+ Bins, plot
% frequen

runStruct = spikeDataBank.(runList{1});
starts = (runStruct.start:binStep:(runStruct.end - binSize))';
xlabelPlot = starts(1:8:length(starts));
stimOnInd = find(starts == 0);
stimEndInd = find(starts == 2800);

for targ_i = 1:length(target)
    figTitle = sprintf('Population Summary for Sliding T tests - %s',  target{targ_i});
    h = figure('NumberTitle', 'off', 'Name', figTitle, 'units', 'normalized', 'outerposition', [0 0 params.plotSize]);
    sgtitle(figTitle);
    for group_i = 2:length(groupingType)
      s(group_i) = subplot(2,1,group_i - 1);
      stretchStarts = targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Starts')};
      stretchLen = targetRunLens{targ_i, group_i, strcmp(statType, 'Sig Run Lengths')};
      
      keepInd = stretchLen>=3;
      stretchStarts = stretchStarts(keepInd);
      stretchLen = stretchLen(keepInd);
      
      histoBins = zeros(1, 163);
      %histoBins = [];
      for stretch_i = 1:length(stretchStarts)
        histoBins(stretchStarts(stretch_i):stretchStarts(stretch_i) + stretchLen(stretch_i)) = histoBins(stretchStarts(stretch_i):stretchStarts(stretch_i) + stretchLen(stretch_i))  + 1;
        %histoBins = [histoBins, stretchStarts(stretch_i):stretchStarts(stretch_i) + stretchLen(stretch_i) - 1];
      end
      plot(histoBins)
      title(groupingType{group_i})
      xlim([0 161])
      hold on
      xlabel('Bin Start time (ms)')
      ylabel('Frequency (Hz)')
      xlim([0,length(histoBins)]);
      xticks(1:8:length(starts))
      xticklabels(xlabelPlot);
      if group_i == 3
        linkaxes([s(2), s(3)],'xy');
      end
      plot([stimOnInd, stimOnInd], ylim(),'k','linewidth',4);
      plot([stimEndInd, stimEndInd], ylim(),'k','linewidth',4);
    end
    saveFigure(params.outputDir, ['2.' figTitle], [], 1, params.exportFig, 0, {''})
    close(h)
end

cohensStatTypePlotLabel = ["Cohen's D Values", "Maximum Cohen's D Values"];
% Figure 2 - the mean Omega Curve.
for targ_i = 1:length(target)
  figTitle = sprintf("Summary Statistic for Cohen's D - %s", target{targ_i});
  h = figure('NumberTitle', 'off', 'Name', figTitle,'units','normalized','outerposition',[0 0 params.plotSize]);
  sgtitle(figTitle);
  for group_i = 2:length(groupingType)
    for plot_i = 1:length(cohensStatType)
      plotInd = ((group_i - 2) * length(cohensStatType)) + plot_i;
      subplot(length(groupingType)-1, length(cohensStatType), plotInd);
      plotData = cohensData{targ_i, group_i, plot_i};
      plotData = [plotData(:)];
      histogram(plotData);
      if plot_i == 1
        ylabel(groupingType{group_i})
      end
      xlabel(cohensStatType{plot_i});
      plotTitle = sprintf('%s (90th p = %s)', cohensStatTypePlotLabel{plot_i}, num2str(round(prctile(plotData, 90), 2)));
      title(plotTitle);
    end
  end
  saveFigure(params.outputDir, strjoin(["2.", figTitle]), [], 1, params.exportFig, 0, {''})
  close(h)
end

end

function [spikeDataBank, frameFiringStruct]  = frameFiringRates(spikeDataBank,params)
% Generates means and distributions based on what is being looking at.
disp('Starting Frame Firing Rate Analysis...');

% Generate variables needed
runList = fields(spikeDataBank);

if ~exist(params.outputDir, 'dir')
  mkdir(params.outputDir)
end

% Step 1 - generate bin times and spike rates.
if ~isfield(spikeDataBank.(runList{end}), 'frameFiringRates')
  for run_ind = 1:length(runList)
    runStruct = spikeDataBank.(runList{run_ind});
    if length(unique(cellfun('length',runStruct.attendedObjData.frameStartInd))) == 1 % If all the videos are the same number of frames
      starts = runStruct.attendedObjData.frameStartInd{1}'+params.delay;
      ends = runStruct.attendedObjData.frameEndInd{1}'+params.delay;
      spikeDataBank.(runList{run_ind}).frameTimes = [starts,ends];
      spikeDataBank.(runList{run_ind}).frameFiringRates = cell(length(starts),1);
      frameFiringRatesTmp = cell(length(starts),1);
      for frame_ind = 1:length(starts)
        % spikeDataBank.run.framingFiringRates{epoch/bin}{channel}{unit}{stim}
        [frameFiringRatesTmp{frame_ind}, ~, ~] = spikeCounter(spikeDataBank.(runList{run_ind}).spikesByEvent, starts(frame_ind), ends(frame_ind));
      end
    else
      error('Different stimuli contain different number of frames. Implement different spikeCounting method.')
    end
    % Step 1.5 - Rearrange frameFiringRatesTmp to match attendedObjVect data {channel}{unit}{stimuli}(trials*frames)
    frameCount = length(frameFiringRatesTmp);
    trialCount = length(frameFiringRatesTmp{1}{1}{1}{1}.counts);
    frameFiringRates = initNestedCellArray(frameFiringRatesTmp{1},'zeros',[trialCount, frameCount]);
    for chan_ind = 1:length(frameFiringRatesTmp{1})
      for unit_ind = 1:length(frameFiringRatesTmp{1}{chan_ind})
        frameFiringRatesUnitTmp = initNestedCellArray(runStruct.attendedObjData.attendedObjVect, 'zeros', [size(runStruct.attendedObjData.attendedObjVect{1})],1);
        for frame_ind = 1:length(frameFiringRatesTmp)
          for stim_ind = 1:length(frameFiringRatesTmp{frame_ind}{chan_ind}{unit_ind})
            if frame_ind == 1
              frameFiringRatesUnitTmp{stim_ind} = zeros(size(runStruct.attendedObjData.attendedObjVect{stim_ind}));
            end
            if params.useRates
            frameFiringRatesUnitTmp{stim_ind}(:,frame_ind) = frameFiringRatesTmp{frame_ind}{chan_ind}{unit_ind}{stim_ind}.rates;
            else
              frameFiringRatesUnitTmp{stim_ind}(:,frame_ind) = frameFiringRatesTmp{frame_ind}{chan_ind}{unit_ind}{stim_ind}.counts;
            end
          end
        end
        frameFiringRates{chan_ind}{unit_ind} = frameFiringRatesUnitTmp;
      end
    end
    spikeDataBank.(runList{run_ind}).frameFiringRates = frameFiringRates;
    disp(run_ind)
  end
  saveSpikeDataBank(spikeDataBank, 2, 'save',fileparts(params.outputDir));
else
  fprintf('Frame firing rates already calculated., continuing... \n');
end

% Step 2 - Go through each run, sum frame counts/rates associated with the
% specific objects being attended.
objList = spikeDataBank.(runList{1}).attendedObjData.objList;
objListBroad = {'Face','Body','Hand','GazeFollow','Object','Bkg'};

if params.broadLabels
  objListPlot = objListBroad;
else
  objListPlot = objList;
end

% Step 3 - Generate vectors with data from all runs, organized as follows
% from ObjFrameFiringRates{channel}{unit}{stim}{obj}{dataType} to
% objFrameFiringRatesTotal{stim}{obj}{groupingType}{dataType}
groupingType = {'Unsorted','Units','MUA'};
dataType = {'Raw','Mean','Run Ind'};

% Pooling across runs requires vector of all stim, extract the eventIDs field, generate a cell array of unique stimuli
allStimuliVec = struct2cell(structfun(@(x) x.eventIDs, spikeDataBank,'UniformOutput', 0));
[allStimuliVec, ~, ic] = unique(vertcat(allStimuliVec{:}));
allStimuliVecCounts = accumarray(ic,1);

objFrameFiringRatesTotal = cell([length(allStimuliVec),length(groupingType), length(objListPlot), length(dataType)]);
%objFrameFiringRatesTotal{stim_ind, group_ind, obj_ind, data_ind};

channelUnitMat = [];
for run_ind = 1:length(runList)
  frameFiringRates = spikeDataBank.(runList{run_ind}).frameFiringRates;
  attendedObjVect = spikeDataBank.(runList{run_ind}).attendedObjData.attendedObjVect;
  eventIDs = spikeDataBank.(runList{run_ind}).eventIDs;
  chanCount = 0;
  unitCount = 0;
  
  % Per Stimuli Plot
  for channel_ind = 1:length(frameFiringRates)
    chanCount = chanCount + 1;
    for unit_ind = 1:length(frameFiringRates{channel_ind})
      unitCount = unitCount + 1;
      if unit_ind == 1
        groupInd = 1;  %Unsorted
      elseif unit_ind == length(frameFiringRates{channel_ind}) % MUA
        groupInd = 3;  %MUA
      else
        groupInd = 2;  %Unit
      end
      for stim_ind = 1:length(attendedObjVect)
        stimStoreInd = strcmp(spikeDataBank.(runList{run_ind}).eventIDs{stim_ind}, allStimuliVec);
        stimAttendedObj = attendedObjVect{stim_ind};
        spikeStruct = frameFiringRates{channel_ind}{unit_ind}{stim_ind};
        %ObjFrameFiringRatesStim = initNestedCellArray([length(objListPlot), length(dataType)],'zeros');
        for obj_ind = 1:length(objList)
          objRates = spikeStruct(strcmp(stimAttendedObj, objList{obj_ind}));
          
          % Identify correct object group, when using broad labels.
          if params.broadLabels
            switch find(strcmp(objList,objList{obj_ind}))
              case {1,7};         objStoreInd = 1;
              case {2,8};         objStoreInd = 2;
              case {3,4,9,10};    objStoreInd = 3;
              case {5,11};        objStoreInd = 4;
              case {6,12};        objStoreInd = 5;
              case 13;            objStoreInd = 6;
            end
          else
            objStoreInd = obj_ind;
          end
          if ~isempty(objRates)
            %ObjFrameFiringRatesStim{objStoreInd}{1} = [ObjFrameFiringRatesStim{objStoreInd}{1}; objRates];
            objFrameFiringRatesTotal{stimStoreInd, groupInd, objStoreInd, 1} = [objFrameFiringRatesTotal{stimStoreInd, groupInd, objStoreInd, 1}; objRates];
            objFrameFiringRatesTotal{stimStoreInd, groupInd, objStoreInd, 2} = [objFrameFiringRatesTotal{stimStoreInd, groupInd, objStoreInd, 2}; mean(objRates)];
            objFrameFiringRatesTotal{stimStoreInd, groupInd, objStoreInd, 3} = [objFrameFiringRatesTotal{stimStoreInd, groupInd, objStoreInd, 3}; run_ind];
          end
        end
      end
    end
  end
  % Structure is ObjFrameFiringRates{stim}{obj}{channel}{unit}{dataType}
  channelUnitMat = [channelUnitMat; [run_ind, chanCount, unitCount]];
end
disp('Finished calculating objFrameFiringRatesTotal...')

% Generate the correct labels depending on whether broad labels are used
% for objects.
% Broad Label - Swaps labels of individual stimuli for broader catagories,
% as defined in stimParamFile and the parameters.
tmp = load(params.stimParamsFilename);
for event_i = 1:length(tmp.paramArray)
  totalEventIDs{event_i} = tmp.paramArray{event_i}{1}; %per the stimParamFile spec, this is the event ID
end
totalEventIDs = totalEventIDs';
[~, paramSortVec] = ismember(allStimuliVec, totalEventIDs);
paramArray = tmp.paramArray(paramSortVec);

%Iterate across stimuli and assign new labels.
allStimuliVecBroad = cell(length(allStimuliVec),1);
for label_ind = 1:length(allStimuliVec)
  paramStimSet = paramArray{label_ind};
  allStimuliVecBroad{label_ind} = paramStimSet{ismember(paramStimSet,params.broadLabelPool)};
end
[uniqueStimLabels, ~, stimIndex] = unique(allStimuliVecBroad);

% Find out how many of each of the newly labeled stimuli there are
% uniqueStimLabelCounts = cell2mat(arrayfun(@(x)length(find(stimIndex == x)), unique(stimIndex), 'Uniform', false));

% Step 5 - Generate figures
% Produces 72 images
for stim_ind = 1:length(uniqueStimLabels)
  for group_ind = 1:length(groupingType)
    figTitle = sprintf('%s, Per Object Rates, %s - %s', uniqueStimLabels{stim_ind}, dataType{2}, groupingType{group_ind});
    figFilePath = fullfile(params.outputDir, figTitle);
    if ~exist([figFilePath, '.fig'],'file')
      % Fig 1 - 'Chasing, Per Object Rates Means - MUA'
      stimMatInd = strcmp(uniqueStimLabels(stim_ind), allStimuliVecBroad);
      stimObjData = objFrameFiringRatesTotal(stimMatInd, group_ind, :,  2);
      plotArray = cell(length(objListPlot),2);
      removeInd = false(size(stimObjData,3),1);
      for obj_ind = 1:size(stimObjData,3)
        plotArray{obj_ind,1} = vertcat(stimObjData{:,:,obj_ind});
        removeInd(obj_ind) = isempty(plotArray{obj_ind,1});
        plotArray{obj_ind,2} = objListPlot{obj_ind};
      end
      plotArray(removeInd,:) = [];
      h = figure('Name', figTitle, 'NumberTitle', 'off');
      sgtitle(figTitle)
      hold on
      uniqueObjLabels = plotArray(:,2);
      [tmpValues, tmpLabels] = deal([]);
      for obj_ind = 1:length(uniqueObjLabels)
        subplot(2,length(uniqueObjLabels),obj_ind+length(uniqueObjLabels))
        singleObjMeans = plotArray{obj_ind,1};
        histogram(singleObjMeans, 20)
        title(sprintf('%s (\x03bc  = %s)',uniqueObjLabels{obj_ind}, num2str(round(mean(singleObjMeans), 2)) ));
        tmpValues = [tmpValues; singleObjMeans];
        tmpLabels = [tmpLabels; repmat(uniqueObjLabels(obj_ind) ,[length(singleObjMeans),1])];
      end
      subplot(2,length(uniqueObjLabels),[1:length(uniqueObjLabels)])
      boxplot(tmpValues, tmpLabels, 'Symbol', 'o')
      linkaxes(h.Children(2:length(h.Children)-1),'xy')
      savefig(fullfile(params.outputDir, figTitle))
      close(h)
    end 
  end
end

% Produces 36 Images
for obj_ind = 1:length(objListPlot)
  for group_ind = 1:length(groupingType)
    % Fig 2 - 'Mean Face Rates - MUA, All Stimuli'
    figTitle = sprintf('%s %s rates, %s, All Stimuli', dataType{2}, objListPlot{obj_ind}, groupingType{group_ind});
    figFilePath = fullfile(params.outputDir, figTitle);    
    if ~exist([figFilePath, '.fig'],'file')
      % Collect relevant data for plot
      stimObjData = objFrameFiringRatesTotal(:, group_ind, obj_ind,  2);
      plotArray = cell(length(uniqueStimLabels),2);
      removeInd = false(length(uniqueStimLabels),1);
      for stim_i = 1:length(uniqueStimLabels)
        stimMatInd = strcmp(uniqueStimLabels(stim_i), allStimuliVecBroad);
        plotArray{stim_i,1} = vertcat(stimObjData{stimMatInd});
        removeInd(stim_i) = isempty(plotArray{stim_i,1});
        plotArray{stim_i,2} = uniqueStimLabels{stim_i};
      end
      plotArray(removeInd,:) = [];
      
      if ~isempty(plotArray)
        % Set up grid of subplots correctly
        h = figure('Name', figTitle, 'NumberTitle', 'off', 'units','normalized','outerposition',[0 0 1 1]);
        sgtitle(figTitle)
        hold on
        stim2Plot = size(plotArray,1);
        [tmpValues, tmpLabels] = deal([]);
        for stim_ind = 1:stim2Plot
          subplot(2,stim2Plot,stim_ind+stim2Plot)
          stimObjMeans = plotArray{stim_ind,1};
          histogram(stimObjMeans, 20);
          title(sprintf('%s (\x03bc  = %s)',plotArray{stim_ind,2},num2str(round(mean(stimObjMeans), 2))));
          tmpValues = [tmpValues; stimObjMeans];
          tmpLabels = [tmpLabels; repmat(plotArray(stim_ind,2) ,[length(stimObjMeans),1])];
        end
        subplot(2,stim2Plot,1:stim2Plot)
        boxplot(tmpValues, tmpLabels, 'Symbol', 'o')
        linkaxes(h.Children(2:length(h.Children)-1),'xy')
        savefig(fullfile(params.outputDir, figTitle))
        close(h)
      end
    end
  end
end

% Convert to Arrays of run indices into tables
frameFiringStruct.channelUnitMat = channelUnitMat;
frameFiringStruct.objList = objListPlot;
frameFiringStruct.groupList = groupingType;
frameFiringStruct.stimList = uniqueStimLabels;
frameFiringStruct.dataList = dataType(dataInd2Plot);

end

function  [spikeDataBank, stimPSTH] = noveltyAnalysis(spikeDataBank, stimPSTH, meanPSTHStruct, frameFiringStruct, params)
% Function seeks to analyze whether most active PSTHes are enriched from
% novel runs.
% Inputs:
% - spikeDataBank
% - stimPSTH - the output of meanPSTH, which has indexing information in
% meanPSTHStruct.IndStructs.
disp('Starting Novelty Analysis...');

if ~exist(params.outputDir, 'dir')
  mkdir(params.outputDir)
end

stimuliFileNames = meanPSTHStruct.IndStructs{1};
groupTypes = meanPSTHStruct.IndStructs{2};
dataTypes = meanPSTHStruct.IndStructs{3};
stimuliNames = cellfun(@(x) extractBetween(x, 1, length(x)-4), stimuliFileNames);
stimPresMat = meanPSTHStruct.stimPresMat;
epochs = {[300 800], [860 3380], [3680 4100]};
epochNames = {'Fixation','Stim Presentation','Reward'};
dataMetric = {'Max','Avg'};
lineColors = {'k','b'};
figLegend = {'Metric Trace', 'Recording pause > 7d','Presentation Pause > 7d'};
stimPSTHSortMat = cell(size(stimPresMat,1), size(stimPresMat,2), length(epochNames));

for stim_i = 1:length(stimuliNames)
  stimData = stimPSTH(stim_i,:,:);
  for group_i = 1:length(groupTypes)
    psthData = stimData{:,group_i,strcmp(dataTypes, 'PSTH')};
    presCountData = stimData{:,group_i,strcmp(dataTypes, 'presCount')};
    recData = stimData{:,group_i,strcmp(dataTypes, 'daysSinceLastRec')};
    daysPresData = stimData{:,group_i,strcmp(dataTypes, 'daysSinceLastPres')};
    
    % vertical lines denoting long time since last Rec
    recLine = find(recData > 30);
    presLine = find(daysPresData(2:end) > 30);
    if ~isempty(recLine) && ~isempty(presLine)
      recLine = recLine(logical([1;diff(recLine)>1]));
      tmp1 = ones(length(recLine),1);
      tmp2 = ones(length(presLine),1) * 2;
      lineMat = [[recLine;presLine],[tmp1;tmp2]];
      legendInd = [1, length(tmp1)+1, length(lineMat)+1];
    else
      lineMat = [];
      legendInd = [];
    end
    
    
    metaPlotTitle = sprintf('%s - PSTH Values - %s', stimuliNames{stim_i},  groupTypes{group_i});
    h = figure('Name', metaPlotTitle, 'NumberTitle', 'off', 'units','normalized','outerposition',[0 0 1 1]);
    sgtitle(metaPlotTitle);
      
    for epoch_i = 1:length(epochs)
      % Prepare data types of Interest
      epoch = epochs{epoch_i};
      psthMax = max(psthData(:,epoch(1):epoch(2)),[],2);
      psthMean = mean(psthData(:,epoch(1):epoch(2)),2);
      dataMat = [psthMax, psthMean];
      
      % Save sorted matricies
      [~, maxSort] = sort(psthMax);
      [~, meanSort] = sort(psthMean);
      stimPSTHSortMat{stim_i, group_i, epoch_i} = [meanSort, maxSort];

      for data_i = 1:length(dataMetric)
        subPlotTitle = sprintf('%s - %s', dataMetric{data_i}, epochNames{epoch_i});
        plotInd = epoch_i + ((data_i-1) * 3);
        subplot(length(dataMetric), length(epochs), plotInd)
        
        lineArray = gobjects(size(lineMat,1)+1,1);
        lineArray(1) = plot(1:length(presCountData), dataMat(:,data_i));
        title(subPlotTitle);
        xlim([1, length(presCountData)])
        xlabel('Presentation Count')
        ylabel('Rates');
        tickLabels = xticklabels();
        newLabel = cell(length(tickLabels),1);
        for tick_i = 1:length(tickLabels)
          tickVal = round(str2num(tickLabels{tick_i}));
          newLabel{tick_i} = presCountData(tickVal);
        end
        xticklabels(newLabel);
        hold on
        if ~isempty(recLine) && ~isempty(presLine)
          for line_i = 1:size(lineMat,1)
            lineArray(line_i+1) = plot([lineMat(line_i,1), lineMat(line_i,1)], ylim(), 'Color', lineColors{lineMat(line_i,2)}, 'lineWidth',3);
          end
        end
        legend(lineArray(legendInd), figLegend,'AutoUpdate','off')
      end
    end
    savefig(h, fullfile(params.outputDir, metaPlotTitle));
    close(h);
  end
  
end

% Figure 1 - iterate through stimuli, grabbing mean and max values across
% presentations in order, plot them on a single line. Mark dates where long
% recording dates were taken with vertical lines.


end
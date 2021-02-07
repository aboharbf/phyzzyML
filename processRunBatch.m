function [] = processRunBatch(varargin)
% applies the processRun function to a set of analysisParamFiles created
% with differente dateSubj and run variables
%   Inputs:
%   - varargin can have the following forms:
%     - one argument of the form paramBuilderFunctionName as a string i.e.
%       'buildAnalysisParamFileSocVid'. User will be prompted to select
%       directory with files.
%     - two arguments, the first as described above, the second of the form 
%       {{'param1', 'arg1'}; {'param2','arg2'}}, which will replace
%       parameters defined in the analysisParamFile and the paramTable.

[~, machine] = system('hostname');
machine = machine(~isspace(machine));

switch machine
  case 'Skytech_FA'
    outputVolumeBatch = 'D:\DataAnalysis';                                            % The output folder for analyses performed.
    dataLog = 'D:\EphysData\Data\analysisParamTable.xlsx';                            % Only used to find recording log, used to overwrite params.
    eventDataPath = 'C:\Users\aboha\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories\eventData.mat';
    subEventBatchStructPath = sprintf('%s/subEventBatchStruct.mat',outputVolumeBatch);
  case 'homeDesktopWork'
    outputVolumeBatch = 'H:\Analyzed';                                            % The output folder for analyses performed.
    dataLog = 'H:\EphysData\Data\analysisParamTable.xlsx';                        % Only used to find recording log, used to overwrite params.
    eventDataPath = 'C:\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories\eventData.mat';
    subEventBatchStructPath = sprintf('%s/subEventBatchStruct.mat',outputVolumeBatch);
  otherwise
    error('Matching machine not found')
end

replaceAnalysisOut = 0;                                                       % This generates an excel file at the end based on previous analyses. Don't use when running a new.
usePreprocessed = 1;                                                          % uses preprocessed version of Phyzzy, only do when changing plotSwitch or calcSwitch and nothing else.
runParallel = 1;                                                              % Use parfor loop to go through processRun. Can't be debugged within the loop.
debugNoTry = 1;                                                               % Allows for easier debugging of the non-parallel loop.

%% Load Appropriate variables and paths
addpath('buildAnalysisParamFileLib');
addpath(genpath('dependencies'));
if ~replaceAnalysisOut
  if nargin == 1 || nargin == 2
    %If there is only 1 file, it loads the analysisParamFile and composes a
    %list from all the data files in the ephysVolume.
    vararginInputs = varargin;
    analysisParamFile = varargin{1};
    load(feval(analysisParamFile));
    varargin = vararginInputs;
    runListFolder = ephysVolume;
    runList = buildRunList(runListFolder, 'nev');    
  elseif nargin >= 3 || nargin == 0
    disp('Must have 1  or 2 inputs.')
    return
  end
  
  %% Create all of the appropriate AnalysisParamFiles as a cell array.
  [analysisParamFileList, analysisParamFileName, outDirList] = deal(cell(0));
  meta_ind = 1;
  
  % Load in a page from an excel sheet in the data directory.
  if exist(dataLog, 'file')
    paramTable = readtable(dataLog,'ReadRowNames', true, 'PreserveVariableNames', true,'Format','auto');
  end
  
  for dateSubj_i = 1:length(runList)
    dateSubject = runList{dateSubj_i}{1};
    for run_i = 1:length(runList{dateSubj_i}{2})
      runNum = runList{dateSubj_i}{2}{run_i};
      
      % Look through paramTable (if available) and replace variables
      % Load variables from paramTable, if the row is present.
      if exist('paramTable', 'var') && any(strcmp(paramTable.Properties.RowNames, {[dateSubject runNum]}))
        paramTableRow = paramTable([dateSubject runNum], :);
        paramTableVars = paramTable.Properties.VariableNames;
        for param_i = 1:width(paramTableRow)
          if isa(paramTableRow.(paramTableVars{param_i}), 'double') && ~isnan(paramTableRow.(paramTableVars{param_i}))
            eval(sprintf('%s = %d;', paramTableVars{param_i}, paramTableRow.(paramTableVars{param_i})))
          elseif isa(paramTableRow.(paramTableVars{param_i}), 'cell')
            eval(sprintf('%s = %s;', paramTableVars{param_i}, paramTableRow.(paramTableVars{param_i}){1}))
          end
        end
      end
      
      % Evaluates each row of 2nd argument input as a variable argument pair.
      if nargin == 2
        for arg_i = 1:size(varargin{2},1)
          eval(sprintf('%s = %s;', varargin{2}{arg_i, 1}, num2str(varargin{2}{arg_i, 2})));
        end
      end
      
      % If Channels have been changed, redefine all the relevant channel
      % names
      if exist('paramTableVars', 'var') && any(strcmp(paramTableVars, 'channels2Read'))
        ephysParams.spikeChannels = channels2Read; %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
        ephysParams.lfpChannels = channels2Read;
        ephysParams.channelNames = arrayfun(@(x) {sprintf('Ch%d', x)}, channels2Read);
        ephysParams.lfpChannelScaleBy = repmat(8191/32764, [length(channels2Read),1]); %converts raw values to microvolts
        ephysParams.unitsToUnsort = cell(length(channels2Read),1); %these units will be re-grouped with u0
        ephysParams.unitsToDiscard = cell(length(channels2Read),1); %these units will be considered noise and discarded
      end
      
      % Generate appropriate Paths
      outputVolume = outputVolumeBatch;
      analogInFilename = sprintf('%s/%s/%s%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);   %#ok
      [lfpFilename, photodiodeFilename, lineNoiseTriggerFilename] = deal(sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum));
      spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
      taskFilename = sprintf('%s/%s/%s%s.bhv2',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance
      [outDir, stimSyncParams.outDir, stimSyncParams.outDir, ephysParams.outDir]  = deal(sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum));
      analysisParamFilename = strcat(outDir,analysisParamFilenameStem);
      preprocessedDataFilename = strcat(outDir,preprocessedDataFilenameStem);                     %#ok      
      
      % Generate Directories
      if ~exist(outDir,'dir')
        mkdir(outDir);
      end
      
      % In case difference logfile is being used.
      if ~logical(exist(taskFilename,'file'))
        [A, B, C] = fileparts(taskFilename);
        switch C
          case '.mat'
            taskFilename = [A '/' B '.bhv2'];
          case '.bhv2'
            taskFilename = [A '/' B '.mat'];
        end
      end
      
      % Autodetecting spike channels - If the file is parsed, retrieve channels present
      parsedFolderName = sprintf('%s/%s/%s%s_parsed',ephysVolume,dateSubject,dateSubject,runNum);
      if exist(parsedFolderName,'dir') == 7
        [ephysParams.spikeChannels, ephysParams.lfpChannels, ephysParams.channelNames] = autoDetectChannels(parsedFolderName);
      end

      % Save files
      save(analysisParamFilename);
      outDirList{meta_ind} = outDir;
      analysisParamFileList{meta_ind} = analysisParamFilename;
      analysisParamFileName{meta_ind} = [dateSubject runNum];
      meta_ind = meta_ind + 1;
    end
  end
  outDirList = outDirList';
  analysisParamFileList = analysisParamFileList';
  analysisParamFileName = analysisParamFileName';
  
  %% Process the runs
  [errorMsg, startTimes, endTimes, analysisOutFilename] = deal(cell(size(analysisParamFileList)));
  [errorMsg{:}] = deal('None');
  
  tmp  = load(analysisParamFileList{1}, 'outputVolume');
  outputVolume = tmp.outputVolume;
  
  if usePreprocessed
    for path_ind = 1:length(analysisParamFileList)
      analysisParamFileList{path_ind} = [fileparts(analysisParamFileList{path_ind}) filesep 'preprocessedData.mat'];
    end
  end
  
  if license('test','Distrib_Computing_Toolbox') && runParallel
    
    parfor run_ind = 1:length(analysisParamFileList)
      fprintf('run_ind reads %d... \n', run_ind);
      fprintf('Processing %s... \n', analysisParamFileList{run_ind});
      try
        if usePreprocessed
          [~, analysisOutFilename{run_ind}] = processRun('paramBuilder','buildAnalysisParamFileSocialVids','preprocessed',analysisParamFileList{run_ind});
        else
          [~, analysisOutFilename{run_ind}] = processRun('paramFile', analysisParamFileList{run_ind});
        end
      catch MyErr
        errorMsg{run_ind} = MyErr.message;
      end
      close all;
      fprintf('Done! \n');
    end
    
  else
    
    for run_ind = 1:length(analysisParamFileList)
      fprintf('Processing %s... \n', analysisParamFileList{run_ind});
      startTimes{run_ind} = datetime('now');
      if ~debugNoTry
        
        try
          if usePreprocessed
            [~, analysisOutFilename{run_ind}] = processRun('paramBuilder','buildAnalysisParamFileSocialVids','preprocessed',analysisParamFileList{run_ind});
          else
            [~, analysisOutFilename{run_ind}] = processRun('paramFile', analysisParamFileList{run_ind});
          end
        catch MyErr
          errorMsg{run_ind} = MyErr.message;
        end
        
      else
        
        % This loop is identical to the one above, without the try clause.
        % Should be used for testing problems in batch non-parallel run.
        if usePreprocessed
          [~, analysisOutFilename{run_ind}] = processRun('paramBuilder','buildAnalysisParamFileSocialVids','preprocessed',analysisParamFileList{run_ind});
        else
          [~, analysisOutFilename{run_ind}] = processRun('paramFile', analysisParamFileList{run_ind});
        end
      end
      
      endTimes{run_ind} = datetime('now');
      close all;
      fprintf('Done! \n');
    end
  end
  
  % Append the error messages generated to the 'analyzedData' files
  % produced.
  for run_ind = 1:length(analysisParamFileList)
    analysisFile = [fileparts(analysisParamFileList{run_i}) '/analyzedData.mat'];
    errorMsgRun = errorMsg{run_i};
    if exist(dataLog, 'file')
      save(analysisFile, 'errorMsgRun', '-append');
    else
      save(analysisFile, 'errorMsgRun');
    end
  end
  
end
  
if replaceAnalysisOut
  %% If generating excel...  
  if nargin >= 1
    %If there is only 1 file, it loads the analysisParamFile and composes a
    %list from all the data files in the ephysVolume.
    vararginInputs = varargin;
    analysisParamFile = varargin{1};
    load(feval(analysisParamFile), 'analysisLabel');
  end
  
  disp('Recompiling structures for excel output')
  
  % Find all the param files and the analyzedData files.
  analyzedDataList = dir(fullfile(outputVolumeBatch, '**', analysisLabel,'**','analyzedData.mat'));
  analysisParamList = dir(fullfile(outputVolumeBatch, '**', analysisLabel,'**','AnalysisParams.mat'));
  analysisParamFilename = {analysisParamList.folder}';
  analysisOutFilename = fullfile(analysisParamFilename, 'analyzedData.mat');

  analyzedDataListDir = {analyzedDataList.folder}';
  [~, B] = setdiff(analysisParamFilename, analyzedDataListDir);
  
  errorMsg = cell(length(analysisParamFilename),1);

  % Iterate through the files, sometimes error message is missing for some
  % unknown reason. To avoid crashes, loop below. 
  errorCount = 0;
  for analysis_i = 1:length(analysisOutFilename)
    tmp = load(analysisOutFilename{analysis_i}, 'errorMsgRun');
    if isfield(tmp, 'errorMsgRun')
      errorMsg{analysis_i} = tmp.errorMsgRun;
    else
      errorCount = errorCount + 1;
      errorMsg{analysis_i} = 'ERROR MISSING';
    end
  end
  fprintf('Missing Error Count %d \n', errorCount)
  
  % Geneates run name lists
  breakString = @(x) strsplit(x, filesep);
  joinStrings = @(x) strjoin([x(length(x)-2), x(length(x))],'');
  
  tmpAnalysisParamFileName = cellfun(breakString, analysisParamFilename, 'UniformOutput',0);
  analysisParamFileName = cellfun(joinStrings, tmpAnalysisParamFileName, 'UniformOutput',0);
  
  if ~isempty(B)
    analysisOutFilename(B) = deal(cell(1,1));
    errorMsg(B) = deal({'No analyzedData.mat'});
  end
end

%% Compose Excel sheet output
%analysisOutFilename is now a cell array of the filepaths to the outputs of
%runAnalyses. Cycle through them and extract desired information (# of
%units, significance), which you can add to the output file.
%[analysisParamFileName(ii), startTimes(ii), endTimes(ii), errorMsg(ii), channel, UnitCount, UnsortSESig, UnitSESig, MUASESig, sigUnits, sigUnsorted, sigMUA, sigStimLen, sigStimNames, ANOVASigString, ' ']

% Find the first non-empty Entry to retrieve some general information for
% the loops below

firstEntry = find(~cellfun('isempty', analysisOutFilename), 3);
entry2check = min(length(firstEntry), 3);
firstEntry = firstEntry(entry2check);
if ~isempty(firstEntry)
  tmp = load(analysisOutFilename{firstEntry},'sigStruct');
  assert(isfield(tmp, 'sigStruct'), 'sigStruct not in this analysis param, likely due to rerunning the single run and stopping prior to this structure being made. Rerun whole or pick different run')
  epochType = tmp.sigStruct.IndInfo{1};
  dataType = tmp.sigStruct.IndInfo{3};
  epochCount = length(epochType);
else
  error('No entries');
end

% Generate a cell array and the empty row to be used as a template.
tableVar = {'File_Analyzed', 'Error', 'Channel', 'Unit_Count', 'Pupil Dil', 'spikePupilCorr', 'subEvent Unsorted', 'subEvent Unit', 'subEvent MUA', ...
  'Sig Unit count', 'Sig Unsorted', 'Sig MUA', 'Stimuli_count', 'Stimuli','ANOVA','Other Info'};
emptyTable = cell2table(cell(0,length(tableVar)), 'VariableNames', tableVar);
emptyRow = cell([1, length(tableVar)]);

[emptyRow{strcmp(tableVar, 'Unit_Count')}, emptyRow{strcmp(tableVar, 'Sig Unsorted')}, ...
  emptyRow{strcmp(tableVar, 'Sig Unit count')}, emptyRow{strcmp(tableVar, 'Sig MUA')}, ...
  emptyRow{strcmp(tableVar, 'Stimuli_count')}]  = deal(NaN);

[emptyRow{strcmp(tableVar, 'Pupil Dil')}, emptyRow{strcmp(tableVar, 'spikePupilCorr')}] = deal("");

dataArray = cell(epochCount, 1);
[dataArray{:}] = deal(emptyTable);

% Generate a list of events to be counted across the generation of the
% excel table.
tmp2 = load(eventDataPath);
eventList = tmp2.eventData.Properties.VariableNames';
eventList = [eventList; 'blinks'; 'saccades'];
eventSigPerEpoch = cell(epochCount,1);
eventSigPerEpoch{:} = deal(zeros(length(eventList), 3));
eventSigPerEpochTable = cell(epochCount, length(analysisParamFileName));

%eventSigPerEpoch{epoch}{groupingType}
sigUnitMat = zeros(0, length(eventList));
sigUnitLabels = cell(0,1);
sigUnitSamples = zeros(3, length(eventList));

for ii = 1:length(analysisParamFileName)
  dataRow = emptyRow;
  dataRow(strcmp(tableVar, 'File_Analyzed')) = analysisParamFileName(ii);
  dataRow(strcmp(tableVar, 'Error')) = errorMsg(ii);
  
  % Check for analyzedData.mat to open
  if ~isempty((analysisOutFilename{ii}))
    tmp = load(analysisOutFilename{ii}, 'stimStatsTable','sigStruct','frEpochs', 'subEventSigStruct', 'eyeDataStruct', 'spikePupilCorrStruct'); %Relies on genStats function in runAnalyses.
    
    if isfield(tmp, 'sigStruct')
      channelCount = length(tmp.sigStruct.data);
    else
      channelCount = 1;
    end
    
    for epoch_i = 1:epochCount
      for channel_ind = 1:channelCount
        
        % Null model based significance testing
        if isfield(tmp, 'sigStruct')
          dataRow(strcmp(tableVar, 'Channel')) = tmp.sigStruct.channelNames(channel_ind);
          dataRow{strcmp(tableVar, 'Unit_Count')} = tmp.sigStruct.unitCount(channel_ind);
          dataRow{strcmp(tableVar, 'Sig Unsorted')} = tmp.sigStruct.sigInfo(channel_ind, 1);
          dataRow{strcmp(tableVar, 'Sig Unit count')} = tmp.sigStruct.sigInfo(channel_ind, 2);
          dataRow{strcmp(tableVar, 'Sig MUA')} = tmp.sigStruct.sigInfo(channel_ind, 3);
          
          sigStim = vertcat(tmp.sigStruct.data{channel_ind}{epoch_i,:,strcmp(dataType, 'Stimuli')});
          if ~isempty(sigStim)
            dataRow{strcmp(tableVar, 'Stimuli_count')} = length(sigStim);
            dataRow{strcmp(tableVar, 'Stimuli')} = strjoin(sigStim, '|');
          end
        end
        
        % ANOVA based significance
        if isfield(tmp, 'stimStatsTable')
          unitStructs = tmp.stimStatsTable{channel_ind};
          ANOVASigString = ' ';
          for unit_ind = 1:length(unitStructs)
            anovaSig = 0;
            if isfield(unitStructs{unit_ind}, 'tTest')
              anovaSig = unitStructs{unit_ind}.tTest.pVals{epoch_i} < 0.05;
            end
            ANOVASigString = [ANOVASigString ['[' num2str(unitStructs{unit_ind}.taskModulatedP < 0.05) ';' num2str(anovaSig(1)) ']']];
          end
          
          dataRow{strcmp(tableVar, 'ANOVA')} = ANOVASigString;
          
        end
        
        % Event based analysis significance
        if isfield(tmp, 'subEventSigStruct')
          subEventSigStruct = tmp.subEventSigStruct;
          if ~subEventSigStruct.noSubEvent
            % run the function which parses this structure
            [sigString, allEventSigCount, allEventShownCount, allEventSigExpanded] = returnEventNames(subEventSigStruct, channel_ind, eventList);
            
            eventSigPerEpochTable{epoch_i, ii}{channel_ind} = allEventSigCount;
            % Add to the total counts
            eventSigPerEpoch{epoch_i} = eventSigPerEpoch{epoch_i} + allEventSigCount;
            
            % Fill in value for this data row
            dataRow{strcmp(tableVar, 'subEvent Unsorted')} = sigString{1};
            dataRow{strcmp(tableVar, 'subEvent Unit')} = sigString{2};
            dataRow{strcmp(tableVar, 'subEvent MUA')} = sigString{3};
            
            % Second type of formating, for NDTB - generate a table with
            % values to be added to each unit.
            
            ULabels = convertUnitToName(1:length(subEventSigStruct.testResults{channel_ind}), length(subEventSigStruct.testResults{channel_ind}), 1);
            chLabel = ['Ch' num2str(channel_ind)];
            
            unitLabels = strcat(analysisParamFileName(ii), chLabel, ULabels);
            
            sigUnitLabels = [sigUnitLabels; unitLabels];
            sigUnitMat = [sigUnitMat; allEventSigExpanded'];
            sigUnitSamples = sigUnitSamples + allEventShownCount';
            
          end
        end
        
        % Analysis of pupil dilation
        if isfield(tmp, 'eyeDataStruct')
          catStrings = [];
          for cat_i = 1:length(tmp.eyeDataStruct.pupilStats.catComp)
            if ~isempty(tmp.eyeDataStruct.pupilStats.catStats{cat_i}) && tmp.eyeDataStruct.pupilStats.catStats{cat_i}{1} < 0.05
              statsStr = cellfun(@(x) num2str(x, 3), tmp.eyeDataStruct.pupilStats.catStats{cat_i}, 'UniformOutput', 0);
              statsStr = sprintf('(%s)', strjoin(statsStr, ';'));
              nameStats = strjoin([tmp.eyeDataStruct.pupilStats.catComp{cat_i} " " statsStr]);
              catStrings = [catStrings; nameStats];
            end
          end
          
          if ~isempty(catStrings)
            dataRow{strcmp(tableVar, 'Pupil Dil')} = strjoin(catStrings, '| ');
          end
          
        end
        
        % Analysis of spike Pupil correlation
        if isfield(tmp, 'spikePupilCorrStruct')
          allStimPerUnit = tmp.spikePupilCorrStruct.spikePupilUnitAllStim.pVals{channel_ind};
          stimPerUnit = sum(tmp.spikePupilCorrStruct.spikePupilUnitStim.pVals{channel_ind} < 0.05);
          funInd = 1:length(stimPerUnit);
          
          % Combine
          entryCells = arrayfun(@(x) sprintf("%d(%d)", allStimPerUnit(x), stimPerUnit(x)), funInd);
          dataRow{strcmp(tableVar, 'spikePupilCorr')} = strjoin(entryCells, ' | ');
        end
        
        % If these analyses have succeeded, erase the errors.
        
        % Package Outputs into structure for table
        dataArray{epoch_i} = [dataArray{epoch_i}; dataRow];
                                                                    
      end
    end
    
  else
    
    % Append empty row to each epoch's table.
    for epoch_i = 1:epochCount
      dataArray{epoch_i} = [dataArray{epoch_i}; dataRow]; 
    end
    
  end

end

subEventBatchStruct = struct();
subEventBatchStruct.dataTable = eventSigPerEpochTable;
subEventBatchStruct.sigUnitGrid.sigUnitLabels = sigUnitLabels;
subEventBatchStruct.sigUnitGrid.sigUnitMat = sigUnitMat;
subEventBatchStruct.sigUnitGrid.sigUnitSamples = sigUnitSamples;

subEventBatchStruct.indicies = '{epoch, analysisParamFileName}{channel_ind}(eventList, groupingType)';
subEventBatchStruct.epoch = epochType;
subEventBatchStruct.eventList = eventList;
subEventBatchStruct.groupingType = {'Unsorted', 'Units', 'MUA'};
subEventBatchStruct.analysisParamFileName = analysisParamFileName;

save(subEventBatchStructPath, 'subEventBatchStruct');

if isempty(firstEntry)
  
  fprintf('Done - No Excel sheet to produce')
else
  
  %Save Batch Run Results
  for table_ind = 1:length(dataArray)
    dataArray{table_ind}.("Other Info"){3} = sprintf('%d - %d ms', tmp.frEpochs(table_ind,1), tmp.frEpochs(table_ind,2));
    writetable(dataArray{table_ind},sprintf('%s/BatchRunResults.xlsx',outputVolumeBatch),'Sheet', sprintf('%s Epoch', epochType{table_ind}))
  end
end

% %PDF Summary
% %Remove empty cells from analysisOutFilename
% nonEmptyCellInd = ~(cellfun('isempty',analysisOutFilename));
% analysisOutFilename = analysisOutFilename(nonEmptyCellInd);
% 
% %Find the directory
% figDirPath = @(cell) (fileparts(cell));
% analysisOutFigDirs = cellfun(figDirPath, analysisOutFilename, 'UniformOutput', 0);
% 
% %run through the createSummaryDoc Function.
% for figDir_ind = 1:length(analysisOutFigDirs)
%   createSummaryDoc('buildLayoutParamFile', analysisOutFigDirs{figDir_ind})
% end

end

function [spike, LFP, names] = autoDetectChannels(parsedFolderName)
    channelFiles = dir([parsedFolderName '/*.NC5']);
    channelNames = {channelFiles.name};
    channelsTmp = [];
    channelNameTmp = {};
    for chan_ind = 1:length(channelNames)
      startInd = regexp(channelNames{chan_ind}, 'Ch') + 2;
      stopInd = regexp(channelNames{chan_ind},'.NC5') - 1;
      chanNum = channelNames{chan_ind}(startInd:stopInd);
      channelNameTmp{chan_ind} = ['Ch' chanNum];
      channelsTmp(chan_ind) = str2double(channelNames{chan_ind}(startInd:stopInd));
    end
    [spike, LFP] = deal(channelsTmp); %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
    names = channelNameTmp;
end

function [sigString, allEventSigCount, allEventPerGroupingType, allEventSigExpanded] = returnEventNames(subEventSigStruct, channel_ind, allEventString)
% a function which parses the 'subEventSigStruct' produced by individual
% calls to 'processRun'.
% Input
% - subEventSigStruct, with fields testResults and events.
% - channel_ind - the channel currently being checked.
% - allEventString - a string of all events in the entire data set.
% Output
% - sigString - a cell array with names of events significant in
% {Unsorted}, {Units}, and {MUA}.
% - allEventSigCount - an allEvent * Grouping Type sig count
% - allEventPerGroupingType - an allEvent * Grouping type of possible sig
% units. This is for coming up with the right denominator.
% - allEventSigExpanded - an allEvent * True unit count of significance for
% adding to a large table or for Neural decoding (per true unit).

% Find the appropriate indices for events to add.
chanInfo = subEventSigStruct.testResults{channel_ind};
chanInfo = cell2mat([chanInfo{:}]);
chanSig = chanInfo < 0.05;

% use chanSig as an index to add the appropriate titles for
% subEvents for adding to the row in the excel sheet.
UnsortSESig = strjoin(subEventSigStruct.events(chanSig(:,1)), ', ');
MUASESig = strjoin(subEventSigStruct.events(chanSig(:,end)), ', ');
if size(chanSig,2) > 2
  unitSig = chanSig(:, 2:end-1);
  UnitSESig = strjoin(subEventSigStruct.events(any(unitSig,2)), ', ');
else
  UnitSESig = '';
end

% Add counts to the overall counts across runs
sigString{1} = UnsortSESig;
sigString{2} = UnitSESig;
sigString{3} = MUASESig;

% Produce a count of the number of units with significant activation for a
% particular event, according to the list of events in allEventString
allEventRunSigCount = zeros(length(subEventSigStruct.events), 3);
[allEventSigCount] = deal(zeros(length(allEventString), 3));
[allEventShownCount, allEventSigExpanded] = deal(zeros(length(allEventString), size(chanInfo, 2)));

% trueUnitCount = size(chanSig,2);
% allEventRunSigCount = zeros(trueUnitCount, length(allEventString));
% 
% for event_i = 1:length(subEventSigStruct.events)
%   event_table_ind = strcmp(subEventSigStruct.events{event_i}, allEventString);
%   allEventRunSigCount(:, event_table_ind) = chanSig(event_i,:);
% end

% Change the matrix from localEvent * Units to LocalEvent * Grouptype
allEventRunSigCount(:,1) = chanSig(:,1);
allEventRunSigCount(:,end) = chanSig(:,end);
if size(chanSig, 2) > 2
  allEventRunSigCount(:,2) = sum(chanSig(:, 2:end-1), 2);
end

% Change from LocalEvent * GroupType to AllEvent * groupType
for event_i = 1:length(subEventSigStruct.events)
  event_table_ind = strcmp(allEventString, subEventSigStruct.events{event_i});
  allEventShownCount(event_table_ind, :) = 1;
  for unit_i = 1:3
    allEventSigCount(event_table_ind, unit_i) = allEventRunSigCount(event_i, unit_i);
  end
  
  %Create expanded presentation
  for unit_i = 1:size(chanSig, 2)
    allEventSigExpanded(event_table_ind, unit_i) = chanSig(event_i, unit_i);
  end
end

allEventPerGroupingType = zeros(length(allEventString), 3);
allEventPerGroupingType(:,1) = allEventShownCount(:,1);
allEventPerGroupingType(:,end) = allEventShownCount(:,end);
if size(allEventShownCount, 2) > 2
  allEventPerGroupingType(:,2:end-1) = sum(allEventShownCount(:,2:end-1), 2);
end

end
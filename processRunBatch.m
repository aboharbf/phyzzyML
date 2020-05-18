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

replaceAnalysisOut = 1;                                                       % This generates an excel file at the end based on previous analyses. Don't use when running a new.
outputVolume = 'H:\Analyzed';                                                 % Only used for the excel doc. change in analysisParamFile to change destination.
dataLog = 'C:\Onedrive\Lab\ESIN_Ephys_Files\Data\analysisParamTable.xlsx';    % Only used to find recording log, used to overwrite params.
usePreprocessed = 0;                                                          % uses preprocessed version of Phyzzy, only do when changing plotSwitch or calcSwitch and nothing else.
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
    paramTable = readtable(dataLog,'ReadRowNames', true, 'PreserveVariableNames', true);
  end
  
  for dateSubj_i = 1:length(runList)
    dateSubject = runList{dateSubj_i}{1};
    for run_i = 1:length(runList{dateSubj_i}{2})
      runNum = runList{dateSubj_i}{2}{run_i};
      
      % Look through paramTable (if available) and replace variables
      % Load variables from paramTable, if the row is present.
      if any(strcmp(paramTable.Properties.RowNames, {[dateSubject runNum]}))
        paramTableRow = paramTable([dateSubject runNum], :);
        paramTableVars = paramTable.Properties.VariableNames;
        for param_i = 1:width(paramTableRow)
          if isa(paramTableRow.(paramTableVars{param_i}), 'double')
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
      
      % Generate appropriate Paths
      analogInFilename = sprintf('%s/%s/%s%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);   %#ok
      [lfpFilename, photodiodeFilename, lineNoiseTriggerFilename] = deal(sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum));
      spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
      taskFilename = sprintf('%s/%s/%s%s.bhv2',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance
      [outDir, stimSyncParams.outDir, stimSyncParams.outDir]  = deal(sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum));
      analysisParamFilename = strcat(outDir,analysisParamFilenameStem);
      preprocessedDataFilename = strcat(outDir,preprocessedDataFilenameStem);                     %#ok
      ephysParams.outDir = sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum);
      
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
      startTimes{run_ind} = datetime('now');
      try
        if usePreprocessed
          [~, analysisOutFilename{run_ind}] = processRun('paramBuilder','buildAnalysisParamFileSocialVids','preprocessed',analysisParamFileList{run_ind});
        else
          [~, analysisOutFilename{run_ind}] = processRun('paramFile', analysisParamFileList{run_ind});
        end
        errorMsg{run_ind} = 'None';
      catch MyErr
        errorMsg{run_ind} = MyErr.message;
      end
      endTimes{run_ind} = datetime('now');
      close all;
      fprintf('Done! \n');
    end
    
    %Since save can't be called in a parfor loop.
    for run_ind = 1:length(analysisParamFileList)
      errorMsgRun = errorMsg{run_ind};
      save(analysisParamFileList{run_ind}, 'errorMsgRun', '-append');
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
          [errorMsg{run_ind}, errorMsgRun] = deal('None');
        catch MyErr
          [errorMsg{run_ind}, errorMsgRun] = deal(MyErr.message);
        end
        
      else
        % This loop is identical to the one above, without the try clause.
        if usePreprocessed
          [~, analysisOutFilename{run_ind}] = processRun('paramBuilder','buildAnalysisParamFileSocialVids','preprocessed',analysisParamFileList{run_ind});
        else
          [~, analysisOutFilename{run_ind}] = processRun('paramFile', analysisParamFileList{run_ind});
        end
        [errorMsg{run_ind}, errorMsgRun] = deal('None');
      end
      save(analysisParamFileList{run_ind}, 'errorMsgRun', '-append');
      endTimes{run_ind} = datetime('now');
      close all;
      fprintf('Done! \n');
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
  analyzedDataList = dir(fullfile(outputVolume, '**', analysisLabel,'**','analyzedData.mat'));
  analysisParamList = dir(fullfile(outputVolume, '**', analysisLabel,'**','AnalysisParams.mat'));
  
  analyzedDataListDir = {analyzedDataList.folder}';
  analysisParamFilename = {analysisParamList.folder}';
  [~, B] = setdiff(analysisParamFilename, analyzedDataListDir);
  
  
  errorMsg = cell(length(analysisParamFilename),1);
  % If preprocessing was used, the error message is in a different spot.
  if usePreprocessed
    analysisList = fullfile(analyzedDataListDir, 'preprocessedData.mat');
  else
    analysisList = fullfile(analyzedDataListDir, 'AnalysisParams.mat');
  end
  
  % Iterate through the files, sometimes error message is missing for some
  % unknown reason. To avoid crashes, loop below. 
  errorCount = 0;
  for analysis_i = 1:length(analysisList)
    tmp = load(analysisList{analysis_i}, 'errorMsgRun');
    if isfield(tmp, 'errorMsgRun')
      errorMsg{analysis_i} = tmp.errorMsgRun;
    else
      errorCount = errorCount + 1;
      errorMsg{analysis_i} = 'ERROR MISSING';
    end
  end
  fprintf('Missing Error Count %d \n', errorCount)
  
  [startTimes, endTimes] = deal(cell(length(analyzedDataList),1));
  [startTimes{:}, endTimes{:}] = deal(datetime(now,'ConvertFrom','datenum'));
  
  % Process the strings to generate the run labels and the paths to
  % analyzedData.mat
  breakString = @(x) strsplit(x, filesep);
  joinStrings = @(x) strjoin([x(length(x)-2), x(length(x))],'');
  
  tmpAnalysisParamFileName = cellfun(breakString, analysisParamFilename, 'UniformOutput',0);
  analysisParamFileName = cellfun(joinStrings, tmpAnalysisParamFileName, 'UniformOutput',0);
  analysisOutFilename = fullfile(analysisParamFilename, 'analyzedData.mat');
  analysisOutFilename(B) = deal(cell(1,1));
  errorMsg(B) = deal({'No analyzedData.mat'});
end

%% Compose Excel sheet output
%analysisOutFilename is now a cell array of the filepaths to the outputs of
%runAnalyses. Cycle through them and extract desired information (# of
%units, significance), which you can add to the output file.
%[analysisParamFileName(ii), startTimes(ii), endTimes(ii), errorMsg(ii), channel, UnitCount, UnsortSESig, UnitSESig, MUASESig, sigUnits, sigUnsorted, sigMUA, sigStimLen, sigStimNames, ANOVASigString, ' ']

% Find the first non-empty Entry to retrieve some general information for
% the loops below
%save('postBatchRun')
%load('postBatchRun')

firstEntry = find(~cellfun('isempty', analysisOutFilename), 1);
if ~isempty(firstEntry)
  tmp = load(analysisOutFilename{firstEntry},'sigStruct');
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
            chanInfo = subEventSigStruct.testResults{channel_ind};
            chanInfo = cell2mat([chanInfo{:}]);
            chanSig = chanInfo < 0.05;
            UnsortSESig = strjoin(subEventSigStruct.events(chanSig(:,1)), ', ');
            MUASESig = strjoin(subEventSigStruct.events(chanSig(:,end)), ', ');
            if size(chanSig,2) > 2
              unitSig = chanSig(:, 2:end-1);
              UnitSESig = strjoin(subEventSigStruct.events(any(unitSig,2)), ', ');
            end
            
            dataRow{strcmp(tableVar, 'subEvent Unit')} = UnitSESig;
            dataRow{strcmp(tableVar, 'subEvent Unsorted')} = UnsortSESig;
            dataRow{strcmp(tableVar, 'subEvent MUA')} = MUASESig;
            
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

if isempty(firstEntry)
  fprintf('Done - No Excel sheet to produce')
else
  %Save Batch Run Results
  for table_ind = 1:length(dataArray)
    dataArray{table_ind}.("Other Info"){3} = sprintf('%d - %d ms', tmp.frEpochs(table_ind,1), tmp.frEpochs(table_ind,2));
    writetable(dataArray{table_ind},sprintf('%s/BatchRunResults.xlsx',outputVolume),'Sheet', sprintf('%s Epoch', epochType{table_ind}))
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

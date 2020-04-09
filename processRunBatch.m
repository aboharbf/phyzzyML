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

replaceAnalysisOut = 0;                                                       % This generates an excel file at the end based on previous analyses. Don't use when running a new.
outputVolume = 'C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\Analyzed';          % Only used for the excel doc. change in analysisParamFile to change destination.
dataLog = 'C:\Onedrive\Lab\ESIN_Ephys_Files\Data\analysisParamTable.xlsx';    % Only used to find recording log, used to overwrite params.
usePreprocessed = 0;                                                          % uses preprocessed version of Phyzzy, only do when changing plotSwitch or calcSwitch and nothing else.
runParallel = 1;                                                              % Use parfor loop to go through processRun. Can't be debugged within the loop.
debugNoTry = 1;

%% Load Appropriate variables and paths
addpath(genpath('D:\Onedrive\Lab\ESIN_Ephys_Files\Analysis\phyzzy'))
if ~replaceAnalysisOut
  if nargin == 1 || 2
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
  [analysisParamFileList, analysisParamFileName] = deal(cell(0));
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
      analysisParamFileList{meta_ind} = analysisParamFilename;
      analysisParamFileName{meta_ind} = [dateSubject runNum];
      meta_ind = meta_ind + 1;
    end
  end
  
  analysisParamFileList = analysisParamFileList';
  analysisParamFileName = analysisParamFileName';
  
  %% Process the runs  
  [errorMsg, startTimes, endTimes, analysisOutFilename] = deal(cell(length(analysisParamFileList),1));
  
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
  
else
  disp('Recompiling structures for excel output')
  addEnd = @(x) fullfile(x, 'analyzedData.mat');
  breakString = @(x) strsplit(x, filesep);
  joinStrings = @(x) strjoin([x(length(x)-2), x(length(x))],'');
  
  analysisOutFilename = dir([outputVolume '\**\analyzedData.mat']);
  analysisParamList = dir([outputVolume filesep '**' filesep 'AnalysisParams.mat']);
  analysisOutFilename = {analysisOutFilename.folder}';
  analysisParamFilename = {analysisParamList.folder}';
  [~, B] = setdiff(analysisParamFilename, analysisOutFilename);
  analysisParamFilename(B) = [];
  if usePreprocessed
    tmp = cellfun(@(x) load(fullfile(x, 'preprocessedData.mat'), 'errorMsgRun'), analysisParamFilename, 'UniformOutput', false);
  else
    tmp = cellfun(@(x) load(fullfile(x, 'AnalysisParams.mat'), 'errorMsgRun'), analysisParamFilename, 'UniformOutput', false);
  end
  errorMsg = cellfun(@(x) x.errorMsgRun, tmp, 'UniformOutput', false);
  [startTimes, endTimes] = deal(cell(length(analysisOutFilename),1));
  
  tmpAnalysisParamFileName = cellfun(breakString, analysisOutFilename, 'UniformOutput',0);
  analysisParamFileName = cellfun(joinStrings, tmpAnalysisParamFileName, 'UniformOutput',0);
  analysisOutFilename = cellfun(addEnd, analysisOutFilename,'UniformOutput',0);
  [startTimes{:}, endTimes{:}] = deal(datetime(now,'ConvertFrom','datenum'));
end

%analysisOutFilename is now a cell array of the filepaths to the outputs of
%runAnalyses. Cycle through them and extract desired information (# of
%units, significance), which you can add to the output file.
%[analysisParamFileName(ii), startTimes(ii), endTimes(ii), errorMsg(ii), channel, UnitCount, UnsortSESig, UnitSESig, MUASESig, sigUnits, sigUnsorted, sigMUA, sigStimLen, sigStimNames, ANOVASigString, ' ']

firstEntry = find(~cellfun('isempty', analysisOutFilename), 1);
if ~isempty(firstEntry)
  tmp = load(analysisOutFilename{1},'sigStruct');
  epochType = tmp.sigStruct.IndInfo{1};
  groupingType = tmp.sigStruct.IndInfo{2};
  dataType = tmp.sigStruct.IndInfo{3};
end

titles = {'File_Analyzed', 'Start_Time', 'End_time', 'Error', 'Channel', 'Unit_Count', 'SubEvent Unsorted', 'SubEvent Unit', 'SubEvent MUA', 'Sig Unit count', 'Sig Unsorted', 'Sig MUA', 'Stimuli_count', 'Stimuli','ANOVA','Other Info'};
sigMUAInd = strcmp(titles,'Sig MUA');
sigUnsortedInd = strcmp(titles,'Sig Unsorted');

epochCount = size(tmp.sigStruct.data{1},1);
table = cell(epochCount,1);
[table{:}] = deal(cell(length(analysisOutFilename), length(titles)));
true_ind = 1;

for ii = 1:length(analysisParamFileName)
  if ~isempty((analysisOutFilename{ii}))
    tmp = load(analysisOutFilename{ii}, 'stimStatsTable','sigStruct','frEpochs', 'subEventSigStruct'); %Relies on genStats function in runAnalyses.
    true_ind_page = true_ind;
    for epoch_i = 1:epochCount
      for channel_ind = 1:length(tmp.sigStruct.data)
        
        % Null model based significance testing 
        channel = tmp.sigStruct.channelNames(channel_ind);
        UnitCount = tmp.sigStruct.unitCount(channel_ind);
        sigUnsorted = tmp.sigStruct.sigInfo(channel_ind, 1);
        sigUnits = tmp.sigStruct.sigInfo(channel_ind, 2);
        sigMUA = tmp.sigStruct.sigInfo(channel_ind, 3);
        
        sigStim = vertcat(tmp.sigStruct.data{channel_ind}{epoch_i,:,strcmp(dataType, 'Stimuli')});
        sigStimLen = length(sigStim);
        if ~isempty(sigStim)
          sigStimNames = strjoin(sigStim, ' ');
        else
          sigStimNames = ' ';
        end
        
        % ANOVA based significance
        unitStructs = tmp.stimStatsTable{channel_ind};
        ANOVASigString = ' ';
        for unit_ind = 1:length(unitStructs)
          anovaSig = 0;
          if isfield(unitStructs{unit_ind}, 'tTest')
            anovaSig = unitStructs{unit_ind}.tTest.pVals{epoch_i} < 0.05;
          end
          ANOVASigString = [ANOVASigString ['[' num2str(unitStructs{unit_ind}.taskModulatedP < 0.05) ';' num2str(anovaSig(1)) ']']];
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
            subECell = [{UnsortSESig}, {UnitSESig}, {MUASESig}];
            for jj = 1:length(subECell)
              if isempty(subECell{jj})
                subECell{jj} = ' ';
              end
            end
          else
            subECell = [{' '}, {' '}, {' '}];
          end
        else
          subECell = [{' '}, {' '}, {' '}];
          
        end
        % Package Outputs into structure for table
        table{epoch_i}(true_ind, :) = [analysisParamFileName(ii), startTimes(ii), endTimes(ii), errorMsg(ii), channel, UnitCount, subECell, sigUnits, sigUnsorted, sigMUA, sigStimLen, sigStimNames, ANOVASigString, {' '}];
        if ii == 1 && epoch_i == 1
          % table{epoch_i}{true_ind, 11} = sprintf('Comparison: %s Vs %s', strjoin(unitStructs{1}.ANOVA.stats.grpnames{1}), strjoin(unitStructs{1}.ANOVA.stats.grpnames{2}));
        end
        true_ind = true_ind + 1;     
      end
      
      if epoch_i ~= length(tmp.sigStruct.IndInfo{1}) % More pages to do = let count reset.
        true_ind = true_ind_page;
      end
      
    end
  else
    true_ind_page = true_ind;
    for epoch_i = 1:epochCount
      table{epoch_i}(true_ind, :) = [analysisParamFileName(ii), startTimes(ii), endTimes(ii), errorMsg(ii), {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {' '}];
      true_ind = true_ind + 1;
      if epoch_i ~= epochCount % More pages to do = let count reset.
        true_ind = true_ind_page;
      end
    end
  end
end

if isempty(firstEntry)
  fprintf('Done - No Excel sheet to produce')
else
  %Save Batch Run Results
  for table_ind = 1:length(table)
    table{table_ind}{3,end} = sprintf('%d - %d ms', tmp.frEpochs(table_ind,1), tmp.frEpochs(table_ind,2));
    T = cell2table(table{table_ind});
    T.Properties.VariableNames = titles;
    writetable(T,sprintf('%s/BatchRunResults.xlsx',outputVolume),'Sheet', sprintf('%s Epoch', epochType{table_ind}))
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

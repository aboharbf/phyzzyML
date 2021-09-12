function [errorStack, errorStackMsg, files2Check] = processRunBatch(varargin)

% applies the processRun function to a set of analysisParamFiles created
% with differente dateSubj and run variables
%   Inputs:
%   - varargin can have the following forms:
%     - one argument of the form paramBuilderFunctionName as a string i.e.
%       'buildAnalysisParamFileSocVid'. User will be prompted to select
%       directory with files.
%     - two arguments, the first as described above, the 2nd in the form of
%     a cell array with desired paradigms {'naturalSocial', 'headTurnCon',
%     'headTurnIso', 'familiarFace'}
%     - Currently commented out: two arguments, the first as described above, the second of the form 
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
  case 'vera.rockefeller.edu'
    outputVolumeBatch = '/Freiwald/faboharb/EphysAnalysis/Analyzed';                                            % The output folder for analyses performed.
    dataLog = '/Freiwald/faboharb/EphysAnalysis/EphysData/analysisParamTable.xlsx';                        % Only used to find recording log, used to overwrite params.
    eventDataPath = '/Freiwald/faboharb/EphysAnalysis/stimDir/Stimuli and Code/SocialCategories/eventData.mat';
    subEventBatchStructPath = sprintf('%s/subEventBatchStruct.mat',outputVolumeBatch);
  otherwise
    error('Matching machine not found')
end

replaceAnalysisOut = 0;                                                       % This generates an excel file at the end based on previous analyses. Don't use when running a new.
usePreprocessed = 1;                                                          % uses preprocessed version of Phyzzy, only do when changing plotSwitch or calcSwitch and nothing else.
runParallel = 1;                                                              % Use parfor loop to go through processRun. Can't be debugged within the loop.
debugMode = 0;                                                                % Allows for easier debugging of the non-parallel loop. 'runParallel' must be false.

%% Load Appropriate variables and paths
addpath('buildAnalysisParamFileLib');
addpath(genpath('dependencies'));

if ~replaceAnalysisOut
  if nargin == 1 || nargin == 2 || nargin == 3
    %If there is only 1 file, it loads the analysisParamFile and composes a
    %list from all the data files in the ephysVolume.
    vararginInputs = varargin;
    analysisParamFile = varargin{1};
    analysisParamPath = feval(analysisParamFile);
    load(analysisParamPath);
    
    varargin = vararginInputs;
    runListFolder = ephysVolume;
    
    if nargin == 2
      runList = buildRunList(runListFolder, 'paradigm', varargin{2});
    elseif nargin == 3
      runList = buildRunList(runListFolder, 'paradigm', varargin{2}, varargin{3});
    elseif nargin ~= 2
      try
        runList = buildRunList(runListFolder, 'nev');
      catch
        error('Nothing selected');
      end
    end
  elseif nargin >= 4 || nargin == 0
    disp('Must have 1  or 2 inputs.')
    return
  end
  
  %% Create all of the appropriate AnalysisParamFiles as a cell array.
  analysisCount = sum(cellfun(@(x) length(x{2}), runList));
  [analysisParamFileList, analysisParamFileName, outDirList] = deal(cell(analysisCount,1));
  meta_ind = 1;
  
  % Load in a page from an excel sheet in the data directory.
  if exist(dataLog, 'file')
    paramTable = readtable(dataLog,'ReadRowNames', true, 'PreserveVariableNames', true,'Format','auto');
  else
    paramTable = [];
  end
  
  for dateSubj_i = 1:length(runList)
    dateSubject = runList{dateSubj_i}{1};
    for run_i = 1:length(runList{dateSubj_i}{2})
      runNum = runList{dateSubj_i}{2}{run_i};
      dateSubjRun = [dateSubject runNum];
      
      % Look through paramTable (if available) and replace variables
      % Load variables from paramTable, if the row is present.
      if ~isempty(paramTable) && any(strcmp(paramTable.Properties.RowNames, {dateSubjRun}))
        paramTableRow = paramTable(dateSubjRun, :);
        paramTableVars = paramTable.Properties.VariableNames';
        for param_i = 1:width(paramTableRow)
          if isa(paramTableRow.(paramTableVars{param_i}), 'double') && ~isnan(paramTableRow.(paramTableVars{param_i}))
            eval(sprintf('%s = %d;', paramTableVars{param_i}, paramTableRow.(paramTableVars{param_i})))
          elseif isa(paramTableRow.(paramTableVars{param_i}), 'cell')
            eval(sprintf('%s = %s;', paramTableVars{param_i}, paramTableRow.(paramTableVars{param_i}){1}))
          end
        end
      end
      
      % Evaluates each row of 2nd argument input as a variable argument pair.
%       if nargin == 3
%         for arg_i = 1:size(varargin{2},1)
%           eval(sprintf('%s = %s;', varargin{2}{arg_i, 1}, num2str(varargin{2}{arg_i, 2})));
%         end
%       end
      
      % Update the psthPre w/ the ITI. 
      if exist('paramTableVars', 'var') && any(strcmp(paramTableVars, 'psthParams.psthPre'))
                        
        % Update timing related variables for ANOVA analysis
        preFix = [-500 0];
        Fix = [-800 0];
        stimEarly = [0 500];
        stimLate = [500 psthParams.psthImDur];
        reward = [psthParams.psthImDur psthParams.psthImDur + 350];
        epochTargParams.times = [preFix; Fix; stimEarly; stimLate; reward];
        
        epochSWparams.startTime = -psthParams.psthPre;
        
        % Update subevent Analysis stuff
        subEventAnalysisParams.psthParams.movingWin = psthParams.movingWin;
        subEventAnalysisParams.stimPlotParams.psthPre = psthParams.psthPre;
        subEventAnalysisParams.stimPlotParams.psthImDur = psthParams.psthImDur;
        subEventAnalysisParams.stimPlotParams.psthPost = psthParams.psthPost;
        
        % Update lfpalign stuff
        lfpAlignParams.msPreAlign = psthParams.psthPre+tfParams.movingWin(1)/2;
        lfpAlignParams.msPostAlign = psthParams.psthImDur+psthParams.psthPost+tfParams.movingWin(1)/2;
        
        spikeAlignParams.preAlign = psthParams.psthPre+3*psthParams.smoothingWidth;
        spikeAlignParams.postAlign = psthParams.psthImDur+psthParams.psthPost+3*psthParams.smoothingWidth;   %#ok
                
        % Note to self - this stuff should be turned into a function handle
        % and evaluated in runAnalysis.
      end
      
      % If Channels have been changed, redefine all the relevant channel
      % names, if not, check for what was parsed before.
      parsedFolderName = sprintf('%s/%s/%s%s_parsed',ephysVolume,dateSubject,dateSubject,runNum);
      if exist('paramTableVars', 'var') && any(strcmp(paramTableVars, 'channels2Read'))
        
        ephysParams.spikeChannels = channels2Read; %note: spikeChannels and lfpChannels must be the same length, in the same order, if analyzing both
        ephysParams.lfpChannels = channels2Read;
        ephysParams.channelNames = arrayfun(@(x) {sprintf('Ch%d', x)}, channels2Read);
        ephysParams.lfpChannelScaleBy = repmat(8191/32764, [length(channels2Read),1]); %converts raw values to microvolts
        ephysParams.unitsToUnsort = cell(length(channels2Read),1); %these units will be re-grouped with u0
        ephysParams.unitsToDiscard = cell(length(channels2Read),1); %these units will be considered noise and discarded
        
      elseif exist(parsedFolderName,'dir') == 7 && ~isempty(dir(parsedFolderName))
        % Autodetecting spike channels - If the file is parsed, retrieve channels present
        [ephysParams.spikeChannels, ephysParams.lfpChannels, ephysParams.channelNames] = autoDetectChannels(parsedFolderName);
      end
            
      % Generate appropriate Paths
      outputVolume = outputVolumeBatch;
      analogInFilename = sprintf('%s/%s/%s%s.ns2',ephysVolume,dateSubject,dateSubject,runNum);   %#ok
      [lfpFilename, photodiodeFilename, lineNoiseTriggerFilename] = deal(sprintf('%s/%s/%s%s.ns5',ephysVolume,dateSubject,dateSubject,runNum));
      spikeFilename = sprintf('%s/%s/%s%s.nev',ephysVolume,dateSubject,dateSubject,runNum); %note that this file also contains blackrock digital in events
      taskFilename = sprintf('%s/%s/%s%s.bhv2',stimulusLogVolume,dateSubject,dateSubject,runNum); %information on stimuli and performance
      outDir = deal(sprintf('%s/%s/%s/%s/',outputVolume,dateSubject,analysisLabel,runNum));
      
      % Update the many additional parts where outDir is required (function
      % specific plots, make sure to update).
      [eyeStatsParams.outDir, eyeStimOverlayParams.outDir, stimSyncParams.outDir, stimSyncParams.outDir, ephysParams.outDir, subEventAnalysisParams.outDir]  = deal(outDir);
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

      % Save files
      save(analysisParamFilename);
      outDirList{meta_ind} = outDir;
      analysisParamFileList{meta_ind} = analysisParamFilename;
      analysisParamFileName{meta_ind} = [dateSubject runNum];
      meta_ind = meta_ind + 1;
    end
  end
  
  %% Process the runs
  [errorMsg, startTimes, endTimes, analysisOutFilename] = deal(cell(size(analysisParamFileList)));
  [errorMsg{:}] = deal('notRun');
  
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
          [~, analysisOutFilename{run_ind}] = processRun('paramBuilder','buildAnalysisParamFileSocialVids', 'preprocessed',analysisParamFileList{run_ind});
        else
          [~, analysisOutFilename{run_ind}] = processRun('paramFile', analysisParamFileList{run_ind});
        end
        errorMsg{run_ind} = 'None';
        
      catch MyErr
        errorMsg{run_ind} = MyErr;
      end
      close all;
      fprintf('Done! \n');
    end
    
  else
    
    for run_ind = 1:length(analysisParamFileList)
      fprintf('Processing %s... \n', analysisParamFileList{run_ind});
      startTimes{run_ind} = datetime('now');
      if ~debugMode
        try
          if usePreprocessed
            [~, analysisOutFilename{run_ind}] = processRun('paramBuilder','buildAnalysisParamFileSocialVids', 'preprocessed', analysisParamFileList{run_ind});
          else
            [~, analysisOutFilename{run_ind}] = processRun('paramFile', analysisParamFileList{run_ind});
          end
          errorMsg{run_ind} = 'None';

        catch MyErr
          errorMsg{run_ind} = MyErr;
        end
        
      else
        
        % This loop is identical to the one above, without the try clause.
        % Should be used for testing problems in batch non-parallel run.
        if usePreprocessed
          [~, analysisOutFilename{run_ind}] = processRun('paramBuilder','buildAnalysisParamFileSocialVids', 'preprocessed',analysisParamFileList{run_ind});
        else
          [~, analysisOutFilename{run_ind}] = processRun('paramFile', analysisParamFileList{run_ind});
        end
      end
      
      endTimes{run_ind} = datetime('now');
      close all;
      fprintf('Done! \n');
    end
  end
  
end

%% Report on Errors - package outputs

errorInd = ~strcmp(errorMsg, 'None');
if any(errorInd)
  errorMsg = errorMsg(errorInd);
  errorStack = [errorMsg{:}];
  errorStackMsg = {errorStack.message}';
  files2Check = analysisParamFileList(errorInd);
else
  [errorStack, errorStackMsg, files2Check] = deal('None');
end
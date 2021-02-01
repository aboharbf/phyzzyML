function paramTableUpdate()
% Function which updates specific fields on the excel document, to account
% for different paradigms potentially having different parameters
% (specifically psthParams, here).

path2Excel = 'D:\EphysData\Data\analysisParamTable_2020.xlsx';
dataDir = 'D:\EphysData\Data';

% Fields to swap, and the values below for which they should be swapped.
% Must match in length.
fields2Update = ["psthParams.psthPre", "psthParams.psthImDur", "psthParams.psthPost"];

normalTimes = [800 2800 500];
shortTimes = [500; 1000; 500];

% Load in the Excel document
analysisParamTable = readtable(path2Excel,'ReadRowNames', true, 'PreserveVariableNames', true);


% Identify all the runs to find
runList = analysisParamTable.recording;

% Grab the corresponding monkeyLogic files
MonkeyLogicFiles = dir(fullfile(dataDir, '**', '*.bhv2'));
MLFileNamesExt = {MonkeyLogicFiles.name}';
MLFileNames = extractBefore(MLFileNamesExt, '.bhv2');
MLFullPath = fullfile({MonkeyLogicFiles.folder}, {MonkeyLogicFiles.name})';

for ii = 1:length(runList)
  % Load the corresponding monkeyLogic file.
  fileInd = strcmp(runList{ii}, MLFileNames);
  
  if any(fileInd)
    
    [~, ~, C] = mlread(MLFullPath{fileInd});
    
    % Identify the stimuli
    runStim = C.TaskInfo.Stimuli;
    
    % Use the stimuli list (or one or more entries) to ID a paradigm
    if any(contains(runStim, 'monkey'))
      % Soc V NonSoc Paradigm
      time2Use = normalTimes;
    elseif any(contains(runStim, 'headTurnCon'))
      % headTurnCon Paradigm
      time2Use = normalTimes;
    elseif any(contains(runStim, 'headTurnIso'))
      % headTurnIso
      time2Use = shortTimes;
    elseif any(contains(runStim, 'Otis'))
      % FamFace Paradigm
      time2Use = shortTimes;
    end
    
    % Use this paradigm to populate the column of the excel table correctly
    for jj = 1:length(fields2Update)
      analysisParamTable(ii, :).(fields2Update(jj)) = time2Use(jj);
    end
    
  end
end


% Overwrite the original excel file.
writetable(analysisParamTable, path2Excel);
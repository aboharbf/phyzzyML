function paramTableUpdate()
% Function which updates specific fields on the excel document, to account
% for different paradigms potentially having different parameters
% (specifically psthParams, here).

path2Excel = 'D:\EphysData\Data\analysisParamTable.xlsx';
dataDir = 'D:\EphysData\Data';

% Fields to swap, and the values below for which they should be swapped.
% Must match in length.
fields2Update = ["psthParams.psthPre", "psthParams.psthImDur", "psthParams.psthPost"];

normalTimes = [800 2800 500];
shortTimes = [500; 1000; 500];

% Load in the Excel document
analysisParamTable = readtable(path2Excel,'ReadRowNames', true, 'PreserveVariableNames', true,'Format','auto');
tableColumns = analysisParamTable.Properties.VariableNames;

% Identify all the runs to find
runList = analysisParamTable.recording;

% Grab the corresponding monkeyLogic files
MonkeyLogicFiles = dir(fullfile(dataDir, '**', '*.bhv2'));
MLFileNamesExt = {MonkeyLogicFiles.name}';
MLFileNames = extractBefore(MLFileNamesExt, '.bhv2');
MLFullPath = fullfile({MonkeyLogicFiles.folder}, {MonkeyLogicFiles.name})';
[~, MLrunList, ~] = fileparts(MLFullPath);

% % Do only for 2021.
% MLFullPath = MLFullPath(contains(MLrunList, '2021012'));
% MLrunList = MLrunList(contains(MLrunList, '2021012'));
% channel2AddByHand = cell(length(MLrunList),1);

for ii = 1:length(MLrunList)
  
  % See if the existing excel doc has a row for that.
  fileInd = strcmp(runList, MLrunList{ii});

  if any(fileInd) && isnan(analysisParamTable(fileInd, :).(fields2Update(1)))
    addTimes = true;
    rowInd = find(fileInd);
    
  elseif any(fileInd)
    % Do nothing if its on there.
    addTimes = false;
    rowInd = find(fileInd);
  else
    % The run isn't represented on the excel sheet, add it.
    newRow = analysisParamTable(ii,:); % Copy the first row
    newRow.Properties.RowNames =  MLrunList(ii);  % Add the name.
    addTimes = true;
    rowInd = size(analysisParamTable,1)+1;
    analysisParamTable = [analysisParamTable; newRow];  % Add the Row.    
  end
  
  % If you need to populate the row, do so
  if addTimes
  % Load the corresponding monkeyLogic file.
    [~, ~, C] = mlread(MLFullPath{ii});
    
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
      analysisParamTable(rowInd, :).(fields2Update(jj)) = time2Use(jj);
    end
  end
  
  if isnan(analysisParamTable.Channel(rowInd))
    % Find the spike file, fill this in correctly.
    [dirPath, fileN, ~] = fileparts(MLFullPath{ii});
    ns5File = dir(fullfile(dirPath, [fileN '.ns5']));
    if ~isempty(ns5File)
      ns5path = fullfile(ns5File.folder, ns5File.name);
      
      % find channel numbers
      ns5Data = openNSx(ns5path);
      ns5Chan = double([ns5Data.ElectrodesInfo.ElectrodeID]);
      ns5Chan = ns5Chan(ns5Chan < 129);
      
      % Add to the row
      channel2AddByHand{ii} = ns5Chan;
    end
  end
  
end

analysisParamTable = [analysisParamTable.recording, analysisParamTable];
analysisParamTable.Properties.VariableNames{1} = 'recording';

% Overwrite the original excel file.
writetable(analysisParamTable, path2Excel);
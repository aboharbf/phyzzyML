function [trueCellInd, runInfo, unitCountStruct, resultTable, nullSigCell] = trueCellCount(batchRunxls, recordingLogxls)
%% Function
% uses information from batchRunxls and recordingLogxls to determine which 
% units within recordings are new (as opposed to phase 2, or slight changes in depth.)

% batchRunxls = 'D:\DataAnalysis\Jan2020\BatchRunResults.xlsx';
% recordingLogxls = 'D:\Onedrive\Lab\ESIN_Ephys_Files\Data\RecordingsMoUpdated.xlsx';

% Input
% batchRunxls - the path to an xls generated by the processRunBatch, with columns of
% 'File_Analyzed' (in the typical sessionNameRun format), 
% recordingLogxls - the path to a personally annotated log of recording information.
% has to have properly formatted information on grid hole, depth, and
% session name in the correct format. 
performSigCounts = 1;
%% Analysis
% Read both excel tables.
logTable = readtable(recordingLogxls);                % Grab the record with task information from Recordings excel.
analysisTable = readtable(batchRunxls, 'Sheet', 1);   % Grab the first sheet of the analysis table (1 sheet per Epoch).
unitCountStruct = struct();

% Process the task information to find a vector of cells with run title and
% whether this run was a follow up or not.
% Specifically, if a run has the same date, depth, and channel as a previous
% recording, do not include it.
uniqueChannels = unique(logTable.Channel);
recordingDepth = cellfun(@str2double, logTable.putativeDistance);
taskDateRunNum = logTable.recording;  % Repeats represent multiple channels.
taskDate = extractBetween(taskDateRunNum,1,8);
taskDateNum = str2num(cell2mat(taskDate));
% the date changed (different date = 1)
dateChangeInd = [1; logical(diff(taskDateNum))]; %Still works for 2 channels.

%For every recording...
chanDepthChangeInd = zeros(size(dateChangeInd));
for ch_ind = 1:length(uniqueChannels)
  % find the recordings to compare
  matchInd = (logTable.Channel == uniqueChannels(ch_ind));
  % the depth changed < 100 microns from the last record (>100 = 1)
  chanRecordingDepth = recordingDepth(matchInd);
  depthChange = diff(chanRecordingDepth);
  sigDepthChangeInd = [1; abs(depthChange) > 100];
  chanDepthChangeInd(matchInd) = sigDepthChangeInd;
end

%If the run took place on a new day or over 100 microns deeper, count it.
validUnitInd = logical(dateChangeInd+chanDepthChangeInd);

%Use this validUnitInd to generate an index of valid runs for the analysis
%performed, details of which are in the analysisTable.
analyzedTaskDates = analysisTable.File_Analyzed;
unitCounts = analysisTable.Unit_Count;
[~, B] = ismember(analyzedTaskDates,taskDateRunNum); % Find matching entries

% Code Below for finding unique grid holes and their counts.
gridHolesTotal = [logTable.ML logTable.AP];
firstAugRunInd = find(strncmp(taskDateRunNum,'201808',6),1);
gridHolesG1 = [logTable.ML(1:firstAugRunInd) logTable.AP(1:firstAugRunInd)];          % June/July
gridHolesG2 = [logTable.ML(firstAugRunInd+1:end) logTable.AP(firstAugRunInd+1:end)];  % August
unitCountStruct.uniqueGridHolesG1 = unique(gridHolesG1,'rows');
unitCountStruct.uniqueGridHolesG2 = unique(gridHolesG2,'rows');

% Use the indicies of the matches to generate an index for units to count
% in a particular analysis, and to modify the gridHole counts generated.
trueCellInd = validUnitInd(B);
recordingDepth = recordingDepth(B);
gridHoles = gridHolesTotal(B,:);
unitCount = sum(unitCounts(trueCellInd));
validTaskRuns = analyzedTaskDates(trueCellInd);
validGridHoles = gridHoles(trueCellInd,:);

[gridHole, recDepth] = deal(cell(length(gridHoles),1));

for row_i = 1:length(gridHoles)
  gridHole{row_i} = gridHoles(row_i,:);
  recDepth{row_i} = recordingDepth(row_i);
end

runInfo = [analyzedTaskDates, gridHole, recDepth];

analyzedTaskRuns = unique(analyzedTaskDates);

unitCountStruct.nonExclude.runList = analyzedTaskRuns;
unitCountStruct.nonExclude.AugustInd = find(strcmp(cellfun(@(x) extractBetween(x, 1, 6), analyzedTaskRuns), '201808'), 1);
unitCountStruct.nonExclude.Y2019Ind = find(strcmp(cellfun(@(x) extractBetween(x, 1, 6), analyzedTaskRuns), '201909'), 1);
unitCountStruct.nonExclude.unitCount = sum(unitCounts);

validAnalyzedTaskRuns = unique(analyzedTaskDates(trueCellInd));

unitCountStruct.Exclude.runList = validAnalyzedTaskRuns;
unitCountStruct.Exclude.AugustInd = find(strcmp(cellfun(@(x) extractBetween(x, 1, 6), validAnalyzedTaskRuns), '201808'), 1);
unitCountStruct.Exclude.Y2019Ind = find(strcmp(cellfun(@(x) extractBetween(x, 1, 6), validAnalyzedTaskRuns), '201909'), 1);
unitCountStruct.Exclude.unitCount = sum(unitCounts(trueCellInd));


%% Signifiance Counts (Null Model and ANOVA)
if performSigCounts
  % the processRunBatch currently compiles results of Null model testing and ANOVAs performed within
  % each run. ANOVAs are represented as a string of brackets, 1 per unit, including Unsorted and MUA (per phyzzy convention).
  sheetNames = {'Presentation','Fixation','Reward'};
  unitCountStruct.epochs = sheetNames;
  resultTable = cell(3,1);
  lengthPerBlock = 0;
  [unitCountStruct.sigUnitCountNull , unitCountStruct.sigUnitCountANOVA] = deal(zeros(length(sheetNames),1));
  
  for table_ind = 1:3
    % Load the sheet corresponding to the correct epoch.
    analysisTable = readtable(batchRunxls, 'Sheet', table_ind);
    
    % Null Model counts
    sigUnitCounts = analysisTable.SigUnitCount;
    validSigUnitCounts = sigUnitCounts(trueCellInd);
    sigUnitCounts = validSigUnitCounts(logical(validSigUnitCounts));                % significant unit count
    sigTaskRuns = validTaskRuns(logical(validSigUnitCounts));
    sigGridHoles = validGridHoles(logical(validSigUnitCounts),:);
    
    fractionSig = sum(sigUnitCounts)/unitCount;                            % Percent Significant
    
    if table_ind == 1
      nullSigCell = sum(sigUnitCounts);
    end
    % ANOVA counts
    
    % Parse ANOVA strings.
    ANOVAStrings = analysisTable.ANOVA;
    ANOVAStrings = strrep(ANOVAStrings,' ','');
    extractB = @(x) extractBetween(x, 2, length(x)-1);
    ANOVAStrings = cellfun(extractB, ANOVAStrings);
    breakCells = @(x) strsplit(x, '][');
    ANOVACells = cellfun(breakCells, ANOVAStrings,'UniformOutput',0);
    
    % Use trueCellInd to exclude repeated units.
    ANOVACells = ANOVACells(trueCellInd);
    
    % Initialize vectors for counts
    groupType = {'Unsorted','Units','MUA'};
    dataType = {'Total Tally', 'Task Modulated', 'SocialSig', 'sigRun', 'sigHole'};
    countMat = cell(length(groupType),length(dataType));
    [TaskModUnsortedCount, TaskModMUACount, TaskModUnitCount] = deal(0);
    [SocSigUnsortedCount,  SocSigMUACount,  SocSigUnitCount] = deal(0);
    [sigUnitRunIndA, sigUnitRunIndTM, sigUnitGridHole] = deal([]);
    
    % Iteratate accross runs and count.
    [sigANOVAUnitCounts, sigANOVA, sigANOVAGrid] = deal([]);
    for run_ind = 1:length(ANOVACells)
      for cell_ind = 1:length(ANOVACells{run_ind})
        ANOVAResult = ANOVACells{run_ind}{cell_ind};
        % first index identifies task modulation.
        if cell_ind == 1
          groupInd = 1;
        elseif cell_ind == length(ANOVACells{run_ind})
          groupInd = 3;
        else
          groupInd = 2;
        end
        % {'Total Tally', 'Task Modulated', 'SocialSig', 'sigRun', 'sigHole'};
        countMat{groupInd, strcmp(dataType, 'Total Tally')} = [countMat{groupInd, strcmp(dataType, 'Total Tally')}; 1];
        if str2double(ANOVAResult(1))
          countMat{groupInd, strcmp(dataType, 'Task Modulated')} = [countMat{groupInd, strcmp(dataType, 'Task Modulated')}; 1];
        end
        
        % 2nd number (index 3) identifies social selectivity on ANOVA.
        if str2double(ANOVAResult(3))
          countMat{groupInd, strcmp(dataType, 'SocialSig')} = [countMat{groupInd, strcmp(dataType, 'SocialSig')}; 1];
          countMat{groupInd, strcmp(dataType, 'sigRun')} = [countMat{groupInd, strcmp(dataType, 'sigRun')}; validTaskRuns(run_ind)];
          countMat{groupInd, strcmp(dataType, 'sigHole')} = [countMat{groupInd, strcmp(dataType, 'sigHole')}; validGridHoles(run_ind,:)];
        end
      end
%       if addRunAndHoleToList    
%         sigANOVAUnitCounts = [sigANOVAUnitCounts; SocSigUnitCountTmp];
%         sigANOVA = [sigANOVA; validTaskRuns(run_ind)];
%         sigANOVAGrid = [sigANOVAGrid; validGridHoles(run_ind,:)];
%       end
    end
    
    % Save the total counts into a cell for the table.
    resultTable{table_ind} = countMat;
%     resultTable(table_ind).taskMod = unique(validTaskRuns(sigUnitRunIndTM));
%     resultTable(table_ind).sigNullUnitCounts = sigUnitCounts;
%     resultTable(table_ind).sigNullTaskRuns = sigTaskRuns;
%     resultTable(table_ind).sigNullGrid = sigGridHoles;
%     resultTable(table_ind).sigANOVAUnitCounts = sigANOVAUnitCounts;
%     resultTable(table_ind).sigANOVATaskRuns = sigANOVA;
%     resultTable(table_ind).sigANOVAGrid = sigANOVAGrid;

%     % Report
%     if table_ind == 1
%       fprintf('Total Non-repeat units: %d \n',unitCount)
%     end
%     fprintf('Total Null Model %s Sig units: %d \n',sheetNames{table_ind}, sum(sigUnitCounts))
%     fprintf('\n');
%     fprintf('Total task modulated units: %d \n',TaskModUnitCount)
%     fprintf('Total %s units found via ANOVA : %d \n', sheetNames{table_ind},SocSigUnitCount)
%     
%     % output will be an excel table. This number below will be used for proper formating.
%     lengthPerBlock = max([lengthPerBlock; length(sigUnitRunIndA); length(sigTaskRuns)]);
  end
%   
  % Significant unit counts per grid hole for Presentation epoch (table_ind = 1)
%   AugustInd = find(strncmp(resultTable(1).sigNullTaskRuns,'201808',6),1);
%   unitCountStruct.nullSigNullTaskRunsPresJJ = resultTable(1).sigNullTaskRuns(1:AugustInd-1,:);
%   unitCountStruct.nullSigNullTaskRunsPresAug = resultTable(1).sigNullTaskRuns(AugustInd:end,:);
%   unitCountStruct.nullSigGridHolesPresJJ = resultTable(1).sigNullGrid(1:AugustInd-1,:);
%   unitCountStruct.nullSigGridHolesPresAug = resultTable(1).sigNullGrid(AugustInd:end,:);
%   unitCountStruct.nullSigCountPerHolePresJJ = resultTable(1).sigNullUnitCounts(1:AugustInd-1,:);
%   unitCountStruct.nullSigCountPerHolePresAug = resultTable(1).sigNullUnitCounts(AugustInd:end,:);
%   
%   AugustInd = find(strncmp(resultTable(1).sigANOVATaskRuns,'201808',6),1);
%   unitCountStruct.ANOVASigNullTaskRunsPresJJ = resultTable(1).sigANOVATaskRuns(1:AugustInd-1,:);
%   unitCountStruct.ANOVASigNullTaskRunsPresAug = resultTable(1).sigANOVATaskRuns(AugustInd:end,:);
%   unitCountStruct.ANOVASigGridHolesPresJJ = resultTable(1).sigANOVAGrid(1:AugustInd-1,:);
%   unitCountStruct.ANOVASigGridHolesPresAug = resultTable(1).sigANOVAGrid(AugustInd:end,:);
%   unitCountStruct.ANOVASigCountPerHolePresJJ = resultTable(1).sigANOVAUnitCounts(1:AugustInd-1,:);
%   unitCountStruct.ANOVASigCountPerHolePresAug = resultTable(1).sigANOVAUnitCounts(AugustInd:end,:);
%   
%   tableCell = cell((lengthPerBlock + 1) * 2,3);
%   fileParts = split(fileparts(batchRunxls),filesep);
%   resultName = fileParts{end};
%   [tableCell(1,1:3), tableCell(2 + lengthPerBlock ,1:3)] = deal(sheetNames);
%   
%   for result_ind = 1:length(resultTable)
%     columnLength = length(resultTable(result_ind).sigNullTaskRuns);
%     columnLength2 = length(resultTable(result_ind).sigANOVATaskRuns);
%     tableCell(2:columnLength+1,result_ind) = resultTable(result_ind).sigNullTaskRuns;
%     tableCell(lengthPerBlock+3:lengthPerBlock+2+columnLength2,result_ind) = resultTable(result_ind).sigANOVATaskRuns;
%   end
%   
%   T = cell2table(tableCell);
%   writetable(T,fullfile(fileparts(batchRunxls),sprintf('%s_sigRuns.xlsx',resultName)))
end


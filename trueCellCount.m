%% True Cell Count
batchRunxls = 'D:\DataAnalysis\ANOVA_FullTime\BatchRunResults.xlsx';

recordingLogxls = 'D:\Onedrive\Lab\ESIN_Ephys_Files\Data 2018\RecordingsMoUpdated.xlsx';

% Grab the record with task information from Recordings excel.
logTable = readtable(recordingLogxls);

% How many grid holes were examined?
gridHolesTotal = [logTable.ML logTable.AP];
gridHolesG1 = [logTable.ML(1:44) logTable.AP(1:44)];
gridHolesG2 = [logTable.ML(45:end) logTable.AP(45:end)];
uniqueGridHolesG1 = unique(gridHolesG1,'rows');
uniqueGridHolesG2 = unique(gridHolesG2,'rows');

% Process the task information to find a vector of cells with run title and
% whether this run was a follow up or not.
%Specifically, if a run has the same date, depth, and channel as a previous
%recording, do not include it.
uniqueChannels = unique(logTable.Channel);
recordingDepth = logTable.putativeDistance;
taskDateRunNum = logTable.recording;
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

%% ANOVA time
sheetNames = {'Presentation','Fixation','Reward'};
resultTable = struct();
lengthPerBlock = 0;

for table_ind = 1:3
  % Grab the file with the cell counts
  analysisTable = readtable(batchRunxls, 'Sheet', table_ind);
  if table_ind == 1
    resultTable(table_ind).file = batchRunxls;
    analyzedTaskDates = analysisTable.File_Analyzed;
    unitCounts = analysisTable.Unit_Count;
    % Find matching entries
    [~, B] = ismember(analyzedTaskDates,taskDateRunNum);
    % Use the indicies of the matches to generate an index for units to count.
    validUnitInd = validUnitInd(B);
    AnalyzedTaskRuns = taskDateRunNum(B);
    gridHoles = gridHolesTotal(B,:);
    unitCount =  sum(unitCounts(validUnitInd));
    validTaskRuns = AnalyzedTaskRuns(validUnitInd);
    validGridHoles = gridHoles(validUnitInd,:);
  end
  sigUnitCounts = analysisTable.SigUnitCount;
  validSigUnitCounts = sigUnitCounts(validUnitInd);
  sigTaskRuns = validTaskRuns(logical(validSigUnitCounts));
  sigGridHoles = validGridHoles(logical(validSigUnitCounts),:);
  
  sigUnitCount = sum(sigUnitCounts(validUnitInd));
  fractionSig = sigUnitCount/unitCount;   %Percent Sig
    
  %Use Valid Ind to find relevant recording sites here.
  ANOVAStrings = analysisTable.ANOVA;
  ANOVAStrings = strrep(ANOVAStrings,' ','');
  extractB = @(x) extractBetween(x, 2, length(x)-1);
  ANOVAStrings = cellfun(extractB, ANOVAStrings);
  breakCells = @(x) strsplit(x, '][');
  ANOVACells = cellfun(breakCells, ANOVAStrings,'UniformOutput',0);
  
  ANOVACells = ANOVACells(validUnitInd);
  
  TaskModUnsortedCount = 0;
  TaskModMUACount = 0;
  TaskModUnitCount = 0;
  SocSigUnsortedCount = 0;
  SocSigMUACount = 0;
  SocSigUnitCount = 0;
  sigUnitRunIndA = [];
  sigUnitRunIndTM = [];
  sigUnitGridHole = [];

  
  for run_ind = 1:length(ANOVACells)
    for cell_ind = 1:length(ANOVACells{run_ind})
      ANOVAResult = ANOVACells{run_ind}{cell_ind};
      if str2double(ANOVAResult(1))
        if cell_ind == 1
          TaskModUnsortedCount = TaskModUnsortedCount + 1;
        elseif cell_ind == length(ANOVACells{run_ind})
          TaskModMUACount = TaskModMUACount + 1;
        else
          TaskModUnitCount = TaskModUnitCount + 1;
          sigUnitRunIndTM = [sigUnitRunIndTM; run_ind];
        end
      end
      if str2double(ANOVAResult(3))
        if cell_ind == 1
          SocSigUnsortedCount = SocSigUnsortedCount + 1;
        elseif cell_ind == length(ANOVACells{run_ind})
          SocSigMUACount = SocSigMUACount + 1;
        else
          SocSigUnitCount = SocSigUnitCount + 1;
          sigUnitRunIndA = [sigUnitRunIndA; run_ind];
          sigUnitGridHole = [sigUnitGridHole; validGridHoles(run_ind,:)];
        end
      end
    end
  end
  
  % Save
  resultTable(table_ind).taskMod = unique(validTaskRuns(sigUnitRunIndTM));
  resultTable(table_ind).SigANOVA = validTaskRuns(sigUnitRunIndA);
  resultTable(table_ind).SigNull = sigTaskRuns;
  resultTable(table_ind).sigANOVAGrid = sigUnitGridHole;
  resultTable(table_ind).sigNullGrid = sigGridHoles;

  
  lengthPerBlock = max([lengthPerBlock; length(sigUnitRunIndA); length(sigTaskRuns)]);
  % Report
  if table_ind == 1
    fprintf('Total Non-repeat units: %d \n',unitCount)
  end
  fprintf('Total Null Model %s Sig units: %d \n',sheetNames{table_ind}, sigUnitCount)
  fprintf('\n');
  fprintf('Total task modulated units: %d \n',TaskModUnitCount)
  fprintf('Total %s units found via ANOVA : %d \n', sheetNames{table_ind},SocSigUnitCount)
end

[sigHolesA, ~, sigHoleAC] = unique(resultTable(1).sigNullGrid(5:end,:),'row'); %5 starts August, first 4 dates are dates in J/J.
[sigHolesJJ, ~, sigHoleJJC] = unique(resultTable(1).sigNullGrid(1:4,:),'row'); %5 starts August, first 4 dates are dates in J/J.

A_counts = accumarray(sigHoleAC,1);
JJ_counts = accumarray(sigHoleJJC,1);

sigHolesA_Counts = [A_counts, sigHolesA];
sigHolesJJ_Counts = [JJ_counts, sigHolesJJ];

[sigHolesA, ~, sigHoleAC] = unique(resultTable(1).sigANOVAGrid(6:end,:),'row'); %5 starts August, first 4 dates are dates in J/J.
[sigHolesJJ, ~, sigHoleJJC] = unique(resultTable(1).sigANOVAGrid(1:5,:),'row'); %5 starts August, first 4 dates are dates in J/J.

A_counts = accumarray(sigHoleAC,1);
JJ_counts = accumarray(sigHoleJJC,1);

sigHolesA_Counts = [A_counts, sigHolesA];
sigHolesJJ_Counts = [JJ_counts, sigHolesJJ];


distanceFromNew = zeros(size(validTaskRuns));
reset = [1 2 3];
for ii = 1:length(validTaskRuns)
  
end

tableCell = cell((lengthPerBlock + 1) * 2,3);
fileParts = split(resultTable(1).file,filesep);
resultName = fileParts{end-1};
[tableCell(1,1:3), tableCell(2 + lengthPerBlock ,1:3)] = deal(sheetNames);

for result_ind = 1:length(resultTable)
  columnLength = length(resultTable(result_ind).SigNull);
  columnLength2 = length(resultTable(result_ind).SigANOVA);
  tableCell(2:columnLength+1,result_ind) = resultTable(result_ind).SigNull;
  tableCell(lengthPerBlock+3:lengthPerBlock+2+columnLength2,result_ind) = resultTable(result_ind).SigANOVA;
end

T = cell2table(tableCell);
writetable(T,sprintf('%s.xlsx',resultName))

%% Comparing Lists
taskMod = resultTable(1).taskMod;
taskModIntersect = setdiff(resultTable(1).SigNull, taskMod);
DoubleSigPres = intersect(resultTable(1).SigANOVA, resultTable(1).SigNull);
DoubleSigFix = intersect(resultTable(2).SigANOVA, resultTable(2).SigNull);
DoubleSigReward = intersect(resultTable(3).SigANOVA, resultTable(3).SigNull);


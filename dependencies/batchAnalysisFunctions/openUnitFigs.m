function openUnitFigs(unitChannelString, figTag, addFilt, plotDir)
% Function which operates at the end of the report functions.
% Feed in a cell array, where each element has the dateSubj, runNum,
% channel, and unit to open.
% ex. "20201117Mo001_Ch32U1"

if isa(unitChannelString, 'table')
  unitChannelString = strcat(unitChannelString.dateSubj, unitChannelString.runNum, '_', unitChannelString.channel, unitChannelString.unitType);
end

% plotDir = 'D:\DataAnalysis_Old';

allFigs = dir(fullfile(plotDir, '**', figTag));
allFigs = fullfile({allFigs.folder}, {allFigs.name})';

% Extract channels present
ChIndPerTitle = extractBetween(allFigs, 'Ch', '.fig');
[~, b] = regexp(ChIndPerTitle, '\d+(\.\d+)?'); 
numStrLen = cellfun(@(x) x(1), b);

ChPerFig = nan(size(ChIndPerTitle));
for fig_i = 1:length(ChPerFig)
  ChPerFig(fig_i) = str2double(ChIndPerTitle{fig_i}(1:numStrLen(fig_i)));
end

% Channel and unit vec
chanVec = str2double(extractBetween(unitChannelString, 'Ch', 'U'));
unitVec = extractAfter(unitChannelString, 'U');

tmp = extractBefore(unitChannelString, '_');
dateSubjVec = cellfun(@(x) x(1:8), tmp, 'UniformOutput', false);
runNumVec = cellfun(@(x) x(end-2:end), tmp, 'UniformOutput', false);

for string_i = 1:length(unitChannelString)
  
  % For each string, search the plotDir for the psth.
  chInd = ChPerFig == chanVec(string_i);
  unitInd = contains(allFigs, strcat('U', unitVec(string_i)));
  dateSubjInd = contains(allFigs, dateSubjVec(string_i));
  runNumInd = contains(allFigs, runNumVec(string_i));
  if ~isempty(addFilt)
    addFiltInd = contains(allFigs, addFilt);
    plotInd = chInd & unitInd & dateSubjInd & runNumInd & addFiltInd;
  else
    plotInd = chInd & unitInd & dateSubjInd & runNumInd;
  end
  
  unitPlot = allFigs(plotInd);
  for fig_i = 1:length(unitPlot)
    open(unitPlot{fig_i});
  end
  
end
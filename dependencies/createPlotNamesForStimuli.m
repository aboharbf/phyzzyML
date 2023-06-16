function [stimuliPlotList, newSortInd] = createPlotNamesForStimuli(stimuliList)

% If headTurnCon, translate
stimuliList = eventIDMerge(stimuliList);

% first resort
socInd = find(contains(stimuliList, {'Fight', 'Groom', 'Chas', 'Mount'}));
gdInd = find(contains(stimuliList, {'Goal'}));
idleInd = find(contains(stimuliList, {'Idle'}));
objInd = find(contains(stimuliList, {'obj'}));
landInd = find(contains(stimuliList, {'land'}));
sceneInd = find(contains(stimuliList, {'scene'}));

newSortInd = [socInd; gdInd; idleInd; objInd; landInd; sceneInd];

stimuliPlotListTmp = stimuliList(newSortInd);
if contains(stimuliList{1}, '.avi')
  stimuliPlotListTmp = extractBefore(stimuliPlotListTmp, '.avi');
else
  stimuliPlotList = stimuliPlotListTmp;
  return
end

intNum = cellfun(@(x) x(end), stimuliPlotListTmp);
coreTag = extractBefore(stimuliPlotListTmp, '_');

stimuliPlotList = strcat(coreTag, intNum);

end
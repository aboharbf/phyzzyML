function reportSubEventCounts(subEventStructPath, trueCellStruct)
% A function which loads a struct containing selectivity information about
% subEvents across the population, then applies the trueCellInd, and
% generates a table/reports the summed values. 
disp('calculating subEvent selectivity totals')

tmp = load(subEventStructPath);
subEventStruct = tmp.subEventBatchStruct;
subEventData = subEventStruct.dataTable;
subEventData_channelCounts = cellfun(@(x) length(x), subEventData);
subEventStruct_runs = subEventStruct.analysisParamFileName;
tmp2 = arrayfun(@(x) repmat(subEventStruct_runs(x), [subEventData_channelCounts(x), 1]), 1:length(subEventData_channelCounts), 'UniformOutput', false);
subEvent_runs = vertcat(tmp2{:});

trueCell_runs = trueCellStruct.taskDates;
trueCell_ind = trueCellStruct.ind;

validData = subEventData;
for ii = 1:length(subEventStruct_runs)
  trueCell_i = strcmp(trueCell_runs, subEventStruct_runs(ii));
  keep_ind = trueCell_ind(trueCell_i);
  validData{ii} = validData{ii}(keep_ind);
end

% Combine across cells to tally counts for selecitivity
validDataAll = horzcat(validData{:});
validDataAllMat = zeros(size(validDataAll{1},1), size(validDataAll{1},2), length(validDataAll));
for ii = 1:length(validDataAll)
    validDataAllMat(:,:, ii) = validDataAll{ii};
end
eventUnitCounts = sum(validDataAllMat, 3);
validUnitCounts = trueCellStruct.unitCounts(trueCell_ind);
totalCounts = [length(validUnitCounts), sum(validUnitCounts), length(validUnitCounts)];

eventSelectivityPercent = eventUnitCounts./totalCounts * 100;
eventSelTable = array2table(eventSelectivityPercent);
eventSelTable.Properties.VariableNames = subEventStruct.groupingType;
eventSelTable.Properties.RowNames = subEventStruct.eventList;

end

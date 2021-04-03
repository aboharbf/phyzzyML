function selTable = replaceChanNum_1123(selTable)
% Run 20201123Mo001 - 004 has 32 channels, labeled Ch33 - Ch64. This causes
% problems. 

% used in selCount, spikePathBank_to_rasterData, meanPSTH.

runList = selTable.dateSubj;

if any(strcmp(runList, "20201123Mo"))
  runInd2Change = strcmp(runList, "20201123Mo");
  table2Change = selTable(runInd2Change, :);
  % Ch33 - Ch64 --> Ch1 - Ch32
  table2Change.channel = strcat("Ch", string(double(extractAfter(table2Change.channel, "Ch"))-32));
  selTable(runInd2Change, :) = table2Change;
end

end
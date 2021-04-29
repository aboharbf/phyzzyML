function crossParadigmCheck(spikePathBank, batchAnalysisParams)

[selTablePerRun] = spikePathLoad(spikePathBank, {'selTable'}, batchAnalysisParams.spikePathLoadParams);

paradigmList = unique(spikePathBank.paradigmName);
paradigm2Compare = paradigmList(contains(paradigmList, 'headTurn'));
unitType = {'US', 'MUA', digitsPattern};
unitTypePlot = {'US', 'MUA', 'Units'};

for par_i = 1:length(paradigm2Compare)
  
  NS_parInd = strcmp(spikePathBank.paradigmName, 'NaturalSocial');
  HTC_parInd = strcmp(spikePathBank.paradigmName, paradigm2Compare{par_i});
  
  selTable_NS = selTablePerRun(NS_parInd);
  selTable_NS = vertcat(selTable_NS{:});
  selTable_HTC = selTablePerRun(HTC_parInd);
  selTable_HTC = vertcat(selTable_HTC{:});
  
  selTable_NS = expandSelTableComboEvents(selTable_NS, batchAnalysisParams.selParam);
  selTable_HTC = expandSelTableComboEvents(selTable_HTC, batchAnalysisParams.selParam);
  
  selTable_NS = selTable_NS(selTable_NS.headTurnSel_any ~= 0, :);
  selTable_HTC = selTable_HTC(selTable_HTC.headTurnSel_any ~= 0, :);
  
  headTurnRun_NS_Units = strcat(selTable_NS.dateSubj, '_', selTable_NS.channel, '_', selTable_NS.unitType);
  headTurnRun_HTC_Units = strcat(selTable_HTC.dateSubj, '_', selTable_HTC.channel, '_', selTable_HTC.unitType);
    
  for unit_i = 1:length(unitType)
    
    unitIndNS = contains(selTable_NS.unitType, unitType{unit_i});
    unitIndHTC = contains(selTable_HTC.unitType, unitType{unit_i});
    
    NS_unit = headTurnRun_NS_Units(unitIndNS);
    HTC_unit = headTurnRun_HTC_Units(unitIndHTC);
    
    A = length(NS_unit);
    B = length(HTC_unit);
    AB = length(intersect(HTC_unit, NS_unit));
    A = A - AB;
    B = B - AB;
    
    A_fraction = sum(unitIndNS)/length(unitIndNS);
    B_fraction = sum(unitIndHTC)/length(unitIndHTC);
    
    vennX([A, B, AB], 0.01)
    title(sprintf('headTurning %s activity in Natural Social vs %s', unitTypePlot{unit_i}, paradigm2Compare{par_i}))
    
  end
    
end

end
function selCount(spikePathBank, batchAnalysisParams)
% A Function to tally across processed runs and the 'selTable' produced.

% Collect unit selectivity
outputDir = batchAnalysisParams.selParam.outputDir;

% Now, selectivity across paradigms...
% selCountCrossParadigm(spikePathBank, selTablePerRun, batchAnalysisParams);

paradigmList = unique(spikePathBank.paradigmName);

for par_i = 1:length(paradigmList)
  
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  selTableParadigmPerRun = spikePathBank.selTable(pInd);
  
  % Commented code below is for checking specific stimulus sets within
  % naturalSocial paradigm.
    spikePathBankParadigm = spikePathBank(pInd,:);
%   SS1Ind = cellfun(@(x) any(contains(x, '4003')), spikePathBankParadigm.stimuli);
%   SS2Ind = ~SS1Ind;
%   selTableParadigmPerRun = selTableParadigmPerRun(SS2Ind); for checking
    
  % Combine across tables
  selTableParadigm = vertcat(selTableParadigmPerRun{:});
   
  % Simple Counts
%   UnitTypes = batchAnalysisParams.selParam.UnitTypes;
%   UnitTypePlot = batchAnalysisParams.selParam.UnitTypePlot;
%   fprintf(' Paradigm %s has %d runs, which contain ... \n', paradigmList{par_i}, length(unique(strcat(selTableParadigm.dateSubj, selTableParadigm.runNum))))
%   fprintf(' Runs, collected over %d days...\n', length(unique(selTableParadigm.dateSubj)))
%   for unit_i = 1:length(UnitTypes)
%     unitCount = sum(contains(selTableParadigm.unitType, UnitTypes{unit_i}));
%     fprintf('%d %s traces \n', unitCount, UnitTypePlot{unit_i})
%   end
  
  % Find cells with predefined selectivities for the sake of illustrating
  % effects.
  exampleCellFinder(selTableParadigm);
  exampleCellFinderObject(selTableParadigm)
  
  % Make the directory specific to paradigm.
  batchAnalysisParams.selParam.outputDir = fullfile(outputDir, paradigmList{par_i});
  
  % Make a bar plot and distribution describing epoch preferences (Activity
  % above baseline during an epoch).
  EpochPreferenceBarGraphsAndVennDiagram(selTableParadigm, paradigmList{par_i}, batchAnalysisParams.selParam)
  
  % Make bar plots per Epoch Selectivity - comparing an rates within an
  % epoch across categories. Barplots or Venn Diagram.
  selectivityPerEpochBarGraphsAndVennDiagram(selTableParadigm, paradigmList{par_i}, batchAnalysisParams.selParam)
  
  % Make bar plots for fixed events (Reward delivery, Fixation)
  selectivityPerEventBarGraphs(selTableParadigm, paradigmList{par_i}, batchAnalysisParams.selParam)
  
  % Make bar plots counting eye related activity selectivity
  selectivityPerEyeEvent(selTableParadigm, paradigmList{par_i}, batchAnalysisParams.selParam)
  
  % Make bar plots for subevents (Head turn, body turn)
  selectivityPerSubEventBarGraphs(selTableParadigm, paradigmList{par_i}, batchAnalysisParams.selParam)
  
  % cross section of different categories and 
  selectivityPerAnyAndHeadTurn(selTableParadigm, paradigmList{par_i}, batchAnalysisParams.selParam)

  % Now, for each unitType and selectivity, map out gridHole contents.
%   batchAnalysisParams.selParam.paradigm = paradigmList{par_i};
%   selectivityPerGridHole(spikePathBankParadigm, selTableParadigm, batchAnalysisParams)
%   
%   plotCountsOnBrain(selTableParadigm, batchAnalysisParams);

end

end
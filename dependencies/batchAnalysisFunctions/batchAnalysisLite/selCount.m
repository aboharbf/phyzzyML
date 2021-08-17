function selCount(spikePathBank, batchAnalysisParams)
% A Function to tally across processed runs and the 'selTable' produced.

% Now, selectivity across paradigms...
% selCountCrossParadigm(spikePathBank, selTablePerRun, batchAnalysisParams);

paradigmList = unique(spikePathBank.paradigmName);
monkeyList = {'Sam', 'Mo', 'Combo'};
unitList = {'MU', 'U'};
unitTags = {'MUA', digitsPattern};
selParams = batchAnalysisParams.selParam;
mainOutDir = selParams.outputDir;
selParams.unitList = unitList;
selParams.unitTags = unitTags;

for m_i = 3%1:length(monkeyList)
  for par_i = 1:length(paradigmList)
    
    selParams.paradigmTag = paradigmList{par_i};
    selParams.monkeyTag = monkeyList{m_i};
    
    % Filter at the stage of spikePathBank - monkey and paradigm
    pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
    if ~strcmp(monkeyList{m_i}, 'Combo')
      mInd = contains(spikePathBank.Row, monkeyList{m_i});
    else
      mInd = true(size(pInd));
    end
    
    selTableParadigmPerRun = spikePathBank.selTable(pInd & mInd);
    
    % Combine across tables
    selTableParadigm = vertcat(selTableParadigmPerRun{:});
    
    % Set output directory
    selParams.outputDir = fullfile(mainOutDir, paradigmList{par_i}, monkeyList{m_i});
    
    % Now, for each unitType and selectivity, map out gridHole contents.
    % plotCountsOnBrain(selTableParadigm, selParams);
    
    % Cycle through cell types, create the plots
    for unit_i = 1:length(unitList)
      
      % Set output directory
      selParams.outputDir = fullfile(mainOutDir, paradigmList{par_i}, monkeyList{m_i}, unitList{unit_i});
      
      unitIndex = contains(selTableParadigm.unitType, unitTags{unit_i});
      selTableParadigmUnit = selTableParadigm(unitIndex, :);
      selParams.unitTag = unitList{unit_i};
      
      % Find cells with predefined selectivities for the sake of illustrating
      % effects.
      exampleCellFinderSubEvent(selTableParadigmUnit, selParams);
      exampleCellFinderObject(selTableParadigmUnit, selParams)
      
      % jointTuningPlot(selTableParadigm, paradigmList{par_i}, selParams)
      
      % Make a bar plot and distribution describing epoch preferences (Activity
      % above baseline during an epoch).
      epochPreferenceBarGraphsAndVennDiagram(selTableParadigmUnit, selParams)
      
      % Make bar plots per Epoch Selectivity - comparing an rates within an
      % epoch across categories. Barplots or Venn Diagram.
      % Figures likely not needed anymore
      % selectivityPerEpochBarGraphsAndVennDiagram(selTableParadigmUnit, paradigmList{par_i}, selParams)
      
      % Make bar plots for fixed events (Reward delivery, Fixation)
      selectivityPerEventBarGraphs(selTableParadigmUnit, selParams)
      
      % Make bar plots counting eye related activity selectivity
      selectivityPerEyeEvent(selTableParadigmUnit, selParams)
      
      % Make bar plots for subevents (Head turn, body turn)
      selectivityPerSubEventBarGraphs(selTableParadigmUnit, selParams)
      
      % Similar to previous function, redundant in most ways, does binary
      % plots the previous function doesn't do
      % selectivityPerAnyAndHeadTurn(selTableParadigmUnit, selParams)
      
    end
  end
end

end
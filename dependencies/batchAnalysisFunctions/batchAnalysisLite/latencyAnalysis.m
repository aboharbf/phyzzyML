function latencyAnalysis(spikePathBank, batchAnalysisParams)
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
    eventIDsPerRun = spikePathBank.stimuli(pInd & mInd);
    eventIDPerRun = horzcat(eventIDsPerRun{:});
    
    uniqueStim = unique(eventIDPerRun);     % Identify unique ones
    
    % Combine across tables
    selTableParadigm = vertcat(selTableParadigmPerRun{:});
    unitsPerRun = cellfun(@(x) size(x,1), selTableParadigmPerRun);
    runIndPerUnit = [];
    for run_i = 1:length(unitsPerRun)
      runIndPerUnit = [runIndPerUnit; repmat(run_i, [unitsPerRun(run_i) 1])];
    end
    
    % Set output directory
    selParams.outputDir = fullfile(mainOutDir, paradigmList{par_i}, monkeyList{m_i});
    
    % Cycle through cell types, create the plots
    for unit_i = 2%1:length(unitList)
      
      % Set output directory
      selParams.outputDir = fullfile(mainOutDir, paradigmList{par_i}, monkeyList{m_i}, unitList{unit_i});
      
      % Per Stim histogram information
      unitIndex = contains(selTableParadigm.unitType, unitTags{unit_i});
      selTableParadigmUnit = selTableParadigm(unitIndex, :);
      runIndPerUnitUnit = runIndPerUnit(unitIndex);
      selParams.unitTag = unitList{unit_i};
      
      minLatencyPerUnit = selTableParadigmUnit.latency_min;
      latencyPerStimPerUnit = horzcat(selTableParadigmUnit.latency_perStim{:});
      
      % Scatter of R^2 adjusted and slopeacross epochs
      figH = figure('units', 'normalized', 'position', [0.1 0.1 0.43 0.63]);
      histogram(minLatencyPerUnit, 100)
      hold on
      ylim(ylim());
      medData = median(minLatencyPerUnit, 'omitnan');
      plot([medData, medData], ylim(), 'linewidth', 3, 'color', 'k');
      title(sprintf('Minimum Latency Per Unit - %s ms', num2str(medData,3)))
      xlim([1 1000]);
      xlabel('Time (ms)');
      ylabel('Count');
      figH.Children.FontSize = 14;
      
      figH = figure('units', 'normalized', 'position', [0.1 0.1 0.44 0.82]);
      sgtitle('Histogram of Latency (ms) Per Stim + Median values')
      aviInd = strfind(uniqueStim, '.avi');
      stimNumInd = [aviInd{:}]' - 1;
      stimNum = arrayfun(@(x) str2num(uniqueStim{x}(stimNumInd(x))), 1:length(stimNumInd));
      
      uniqueStimPlot =  extractBefore(uniqueStim, '_');
      uniqueStimPlot =  strrep(uniqueStimPlot, 'monkey', '');
      uniqueStimPlot =  lower(uniqueStimPlot');
      uniqueStimPlot = strcat(uniqueStimPlot, string(stimNum));
      
      for stim_i = 1:length(uniqueStim)
        subplot(8,4,stim_i)
        
        % Identify which belong
        stimInd = strcmp(eventIDPerRun, uniqueStim{stim_i});
        run2Take = find(any(stimInd));
        stimIndInRun = any(stimInd,2);
        latency2Sample = ismember(runIndPerUnitUnit, run2Take);
        latencyData = latencyPerStimPerUnit(stimIndInRun,latency2Sample);
        
        histogram(latencyData, 100);
        distMedian = median(latencyData, 'omitnan');
        hold on
        ylim(ylim());
        plot([distMedian, distMedian], ylim(), 'linewidth', 3, 'color', 'k')
        title(sprintf('%s - %s ms', uniqueStimPlot{stim_i}, num2str(mean(distMedian), 3)));
        xlim([1 1000]);
        ylabel('Count');

      end
      
    end
  end
end

end
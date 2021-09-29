function rampingAnalysis(spikePathBank, batchAnalysisParams)
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
    
    % Cycle through cell types, create the plots
    for unit_i = 2%1:length(unitList)
      
      % Set output directory
      selParams.outputDir = fullfile(mainOutDir, paradigmList{par_i}, monkeyList{m_i}, unitList{unit_i});
      
      unitIndex = contains(selTableParadigm.unitType, unitTags{unit_i});
      selTableParadigmUnit = selTableParadigm(unitIndex, :);
      selParams.unitTag = unitList{unit_i};
      
      rampingUnits = selTableParadigmUnit(~cellfun('isempty', selTableParadigmUnit.rampStats), :);
      totalUnitCount = size(selTableParadigmUnit,1);
      unitsCountWithRamps = size(rampingUnits,1);
      
      % Ramp across epochs count.
      rampingData = rampingUnits.rampStats;
      rampingDataPool = vertcat(rampingData{:});
      rampPerUnit = cellfun(@(x) size(x,1), rampingData);
      unitLabels = strcat(rampingUnits.dateSubj, rampingUnits.runNum, '_', rampingUnits.channel, rampingUnits.unitType);

      unitLabelVec = [];
      for ii = 1:length(rampPerUnit)
        unitLabelVec = [unitLabelVec; repmat(unitLabels(ii), [rampPerUnit(ii), 1])];
      end
      
      % identify epoch correctly.
      epochTimes = [-799 1 501 2801 3201];
      epochNames = {'Fix', 'stimEarly', 'stimLate', 'Reward'};

      % Epoch labeling was done incorrectly - fix it.
      realEpochLabel = zeros(size(rampingDataPool,1),1);
      for epoch_i = 1:length(epochTimes)-1
        epochInd = rampingDataPool(:,2) >= epochTimes(epoch_i) & rampingDataPool(:,2) < epochTimes(epoch_i+1);
        realEpochLabel(epochInd) = deal(epoch_i);
      end
      rampingDataPool(:,1) = realEpochLabel;
      
      % Set up the image
      % a function which stacks rasters, and create a line PSTH at the top. for
      % every raster submitted, this raster adds a raster to the bottom of the
      % stack, and another line to the line plot at the top.
      
      % Targets to load
      epochNames = ["Fixation", "Stimulus Early", "Stimulus Late", "Reward"];
      epochLabelXVals = [-800, 1, 500, 2800, 3200];
               
      % Generate Period shading
      ITI = 0;%psthParams.ITI - 100;
      pre = 800 - ITI;
      ImDur = 2800;
      post = 400;
      
      totalTime =  ITI + pre + ImDur + post;
      
      periodImage = zeros(3000, totalTime, 3);
      periodWindows = [ITI+1 ITI+pre; ITI+pre+1 ITI+pre+500; ITI+pre+500 + 1 ITI+pre+ImDur; ITI+pre+ImDur + 1 ITI+pre+ImDur+post];
      periodColor = [0.8 0.6 0.6;...
        0.6 0.8 0.6;...
        0.4 0.4 0.8;...
        0.4 0.4 0.4;];
      lineColorPerEpoch = periodColor - 0.2;

      for epoch_i = 1:length(periodWindows)
        for layer_i = 1:3
          periodImage(:, periodWindows(epoch_i,1):periodWindows(epoch_i,2), layer_i) = deal(periodColor(epoch_i,layer_i));
        end
      end
      
      % Plotting
      % figTitle = sprintf('Ramping activity Per Epoch');
      plotTitle = sprintf('Ramps per Epoch - %s - %s - %s', paradigmList{par_i}, monkeyList{m_i}, unitList{unit_i});
      figH = figure('Name', plotTitle, 'NumberTitle','off', 'Units', 'normalized', 'Position', [0.3443 0.1176 0.4937 0.7185]);
      title(plotTitle);
      
      pImageHand = image(periodImage);
      pImageHand.YData = [-1500 1500];
      pImageHand.XData = [-pre ImDur + post];
      xlim([-800 3200])
      ylim([-1500 1500]);
      hold on
      plot([-800 3200], [0 0], 'LineWidth', 2, 'color', 'k');
      set(gca,'Ydir','normal')
      
      % Create vector to translate b/t times defined in the ramping data
      % and the actual bins.
      binStartVec = min(rampingDataPool(:,2)):25:max(rampingDataPool(:,3));
      binEndVec = -775:25:max(rampingDataPool(:,3));
      
      % convert columns to indices
      [~, rampStart] = ismember(rampingDataPool(:,2), binStartVec);
      [~, rampEnd] = ismember(rampingDataPool(:,3), binEndVec);
                  
      % Plot lines
      lineEpochs = rampingDataPool(:,1);
      lineStart = rampingDataPool(:,2);
      lineEnd = rampingDataPool(:,3);
      
      % Line widths according to ranking priority (R^2 atm, will be SSE) 
      lineWidth = rampingDataPool(:,6);
      lineWidth = lineWidth/min(lineWidth);
      lineWidth = lineWidth/(max(lineWidth)/3);
      
      % Compile stuff for plot
      lineX = [lineStart, lineEnd];
      lineY = [zeros(length(lineStart),1), rampingDataPool(:,4).* rampEnd - rampStart];
      
      hold on
      for line_i = 1:size(rampingDataPool, 1)
        plot(lineX(line_i, :)', lineY(line_i, :)', 'color', lineColorPerEpoch(lineEpochs(line_i),:)', 'LineWidth', lineWidth(line_i))
      end
      
      % Distribution of slopes across epochs
      figH = figure('units', 'normalized', 'position', [0.1 0.1 0.93 0.35]);
      sgtitle(sprintf('Slope Adjusted Per Epoch - %d/%d Units', unitsCountWithRamps, totalUnitCount));
      for epoch_i = 1:length(epochNames)
        subplot(1,4,epoch_i)
        histogram(rampingDataPool(rampingDataPool(:,1) == epoch_i, 4), 40)
        title(epochNames{epoch_i});
        ylabel('Count'); xlabel('Slope')
      end
      
      % Distribution of R^2 adjusted across epochs
      figH = figure('units', 'normalized', 'position', [0.1 0.1 0.93 0.35]);
      sgtitle(sprintf('R^2 Adjusted Per Epoch - %d/%d Units', unitsCountWithRamps, totalUnitCount));
      for epoch_i = 1:length(epochNames)
        subplot(1,4,epoch_i)
        histogram(rampingDataPool(rampingDataPool(:,1) == epoch_i, 6), 40)
        title(epochNames{epoch_i});
        ylabel('Count'); xlabel('R^2 Adjusted')
      end
      
      % Scatter of R^2 adjusted and slopeacross epochs
      figH = figure('units', 'normalized', 'position', [0.1 0.1 0.99 0.47]);
      sgtitle(sprintf('R^2 Adjusted vs Slope Per Epoch - %d/%d Units', unitsCountWithRamps, totalUnitCount));
      for epoch_i = 1:length(epochNames)
        subplot(1,4,epoch_i)
        epochInd = rampingDataPool(:,1) == epoch_i;
        scatter(rampingDataPool(epochInd, 6), rampingDataPool(epochInd, 4), 5, 'filled')
        title(epochNames{epoch_i});
        ylabel('Slope'); xlabel('R^2 Adjusted')
        ylabel('Slope'); xlabel('SSE')
        xlim([0 1])
        ylim([ylim()])
        hold on
        plot([0 0], [-100 100], 'color', 'k', 'linewidth', 2);
        plot([-100 100], [0 0], 'color', 'k', 'linewidth', 2);
      end
      
      % Report topX per epoch
      units2Report = 10;
      for epoch_i = 1:length(epochNames)
        % Title
        epochInd = find(rampingDataPool(:,1) == epoch_i);
        epochRamps = rampingDataPool(epochInd, 4);  % 4 is Slope, 6 is R^2
        [~, epochInd2] = sort(epochRamps, 'descend');
        
        % Create sorted list
        epochInd = epochInd(epochInd2);
        sortedUnitNames = unitLabelVec(epochInd);
        sortedUnitRamps = epochRamps(epochInd2);
        
        % Report on top X
        fprintf('=== Top Ramps in %s === \n', epochNames{epoch_i})
        for ramp_i = 1:units2Report
          fprintf('%s - R^2 %s \n', sortedUnitNames(ramp_i), num2str(sortedUnitRamps(ramp_i), 3));
        end
        
        % Report bottom X
        fprintf('=== Bottom Ramps in %s === \n', epochNames{epoch_i})
        for ramp_i = 0:units2Report-1
          fprintf('%s - R^2 %s \n', sortedUnitNames(end-ramp_i), num2str(sortedUnitRamps(end-ramp_i), 3));
        end
        
        % open the PSTHes.
        % openUnitFigs(sortedUnitNames(1:units2Report), batchAnalysisParams.analysisDirectory)
        % openUnitFigs(sortedUnitNames(end:-1:end-units2Report+1), batchAnalysisParams.analysisDirectory)
        
        % Extra Section - How about event aligned activity?
        switch epochNames{epoch_i}
          case {"Reward"}
            % for rewards, ramping activity may end after or before reward
            % delivery.
            num2Check = 9;
            openUnitFigs(sortedUnitNames(num2Check), 'subEventPSTH_*_reward*',  batchAnalysisParams.analysisDirectory)
            openUnitFigs(sortedUnitNames(num2Check), 'imPSTH_*',  batchAnalysisParams.analysisDirectory)
            fprintf('%s \n', sortedUnitNames(num2Check))

            num2Check = 0;
            openUnitFigs(sortedUnitNames(end-num2Check), 'subEventPSTH_*_reward*',  batchAnalysisParams.analysisDirectory)
            openUnitFigs(sortedUnitNames(end-num2Check), 'imPSTH_*',  batchAnalysisParams.analysisDirectory)
            fprintf('%s \n', sortedUnitNames(end-num2Check))
            
            openUnitFigs("20201202Mo003_Ch29U1", 'subEventPSTH_*_reward*',  batchAnalysisParams.analysisDirectory)
            openUnitFigs("20201202Mo003_Ch29U1", 'imPSTH_*',  batchAnalysisParams.analysisDirectory)
            fprintf('%s \n', "20201202Mo003_Ch29U1")

%           case 'Fixation'
%             openUnitFigs(sortedUnitNames(1:10), 'subEventPSTH_*_fix*',  batchAnalysisParams.analysisDirectory)
%             openUnitFigs(sortedUnitNames(1:10), 'imPSTH_*',  batchAnalysisParams.analysisDirectory)
        end
        
      end

    end
  end
end

end
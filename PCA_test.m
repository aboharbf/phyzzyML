% PCA Notepad
function PCA_prep()

monk = {'Sam', 'Mo'};
plotIndParams.stimParamsFilename = 'C:\Onedrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\stimParamFileLib\StimParamFileSocialVids_Full.mat';

for m_i = 1:length(monk)
  % Paths to datafiles.
  dataFiles = {sprintf('C:\\batchAnalysis%s\\NeuralDecodingTB\\headTurnCon\\rasterData\\rasterData_binned_100ms_bins_100ms_sampled.mat', monk{m_i}),...
    sprintf('C:\\batchAnalysis%s\\NeuralDecodingTB\\NaturalSocial\\rasterData\\rasterData_binned_100ms_bins_100ms_sampled.mat', monk{m_i})};
  paradigm = {'HTC', 'NS'};
  outputDir = sprintf('D:\\DataAnalysis%s\\batchAnalysis\\PCA\\preProc', monk{m_i});
  
  if ~exist(outputDir, 'dir')
    mkdir(outputDir)
  end
  
  % binned_data - N Site cell array, where each entry is trials*bins
  % binned_labels - a struct with some labels, each of which is a N site cell
  % array, where each entry is trials*1.
  % binned_site_info -  a struct with some labels, each of which is a N * 1
  % structure with different values each can take on (many are logical).
  % paradigm = {'NS'};
  unitIndName = {'MUA', 'SocMUA', 'US&U', 'SocUS&U'};
  % unitIndName = {'US&U'};
  labels2Include = {{'chasing', 'fighting', 'grooming', 'mounting', 'goalDirected', 'idle', 'objects', 'scenes'}}; % Label sorting
  % labels2Include = {{'socialInteraction', 'goalDirected', 'idle', 'objects', 'scene'}}; % Label sorting
  
  figStruct.saveFig = 1;      % save the figure in its output directory.
  figStruct.closeFig = 0;     % close the figure once it is saved
  figStruct.exportFig = 1;    % export figure using export_fig.
  figStruct.saveFigData = 0;  % save data with the figure.
  figStruct.noOverWrite = 1;  % If a figure is already there, don't make it again.
  verbosity = 'INFO';         %other options, 'DEBUG', 'VERBOSE';
  
  storeData = cell(length(paradigm), length(unitIndName));
  
  for paradigm_i = 1:length(paradigm)
    
    % load the binned data
    binnedData = load(dataFiles{paradigm_i});
    
    % Smooth it
    binned_data = binnedData.binned_data;
    parfor ii = 1:length(binnedData.binned_data)
      binData = binned_data{ii};
      for jj = 1:size(binData,1)
        binData(jj,:) = conv(binData(jj,:), gausswin(3), 'same');
      end
      binned_data{ii} = binData;
    end
    
    unitType = binnedData.binned_site_info.UnitType';
    
    unitInd1 = contains(unitType, 'M');
    unitInd2 = contains(unitType, 'M') & binnedData.binned_site_info.sVns_any_selInd;
    unitInd3 = ~contains(unitType, 'M');
    unitInd4 = ~contains(unitType, 'M') & binnedData.binned_site_info.sVns_any_selInd;
    
    unitIndStack = {unitInd1; unitInd2; unitInd3; unitInd4};
    
    % Scramble labels
    % unitInd = unitInd(randperm(length(unitInd)));
    
    binStarts = binnedData.binned_site_info.binning_parameters.the_bin_start_times;
    %   starts = binnedData.binned_site_info.binning_parameters.alignment_event_time;
    starts = 1300;
    points_to_label = [-300, 0, 500, 1000, 1500, 2000, 2500, 3000];
    points_for_lines = [0, 2800];
    
    the_bin_start_times_shift = binStarts - starts;
    bins_to_label = interp1(the_bin_start_times_shift, 1:length(the_bin_start_times_shift), points_to_label);
    x_for_lines = interp1(the_bin_start_times_shift, 1:length(the_bin_start_times_shift), points_for_lines);
    
    binLabelInfo.the_bin_start_times_shift = the_bin_start_times_shift;
    binLabelInfo.bins_to_label = bins_to_label;
    binLabelInfo.points_to_label = points_to_label;
    binLabelInfo.x_for_lines = x_for_lines;
    binLabelInfo.points_for_lines = points_for_lines;
    
    trialsPerSite = cellfun('length', binnedData.binned_labels.socialCat)';
    binCount = size(binnedData.binned_data{1}, 2);
    binnedDataPerGroup = nan(sum(trialsPerSite), binCount);
    
    % Grab Related labels
    stimPerSite = binnedData.binned_labels.stimuli;
    
    for unit_Type_i = 1:length(unitIndStack)
      unitInd = unitIndStack{unit_Type_i};
      for group_i = 1:length(labels2Include)
        
        % Look for the indices to include
        plotIndParams.plotLabels = labels2Include(group_i);
        indices2Grab = cell(sum(unitInd), length(labels2Include{group_i}));
        
        % Initialize array for storing data.
        binnedDataArray = cell(length(labels2Include{group_i}),1);
        [binnedDataArray{:}] = deal(binnedDataPerGroup);
        
        % Collect the indices which will be used to collect data.
        MUA_ind = find(unitInd)';
        for site_i = 1:length(MUA_ind)
          
          % Identify trials to be collected from the binned_data.
          stimuliIndPerLabel = plotIndex(stimPerSite{MUA_ind(site_i)}, plotIndParams);
          [A, ~, C] = unique(stimuliIndPerLabel);
          a_counts = accumarray(C, 1); % 1 is 0, and shouldn't be counted.
          counts2Keep = find(A ~= 0)';
          a_counts = a_counts(counts2Keep);
          
          % Grab the number of labels which are
          trialsPerLabel = min(a_counts);
          
          for label_i = counts2Keep
            % Identify trials to collect
            trialInd = C == label_i;
            
            % Collect their indicies
            label_inds = find(trialInd);
            
            % If there are too many trials, grab a random set of trials to reach the minimum
            if length(label_inds) > trialsPerLabel
              label_inds = label_inds(randperm(length(label_inds), trialsPerLabel));
            end
            
            % Store in the larger matrix
            indices2Grab{site_i, A(label_i)} = label_inds;
            
          end
          
        end
        
        % Find the sites
        trialsPerSite2Gather = cellfun('length', indices2Grab(:,1));
        
        % Cycle through each site, and label
        labelCountTotal = length(labels2Include{group_i});
        for site_i = 1:size(indices2Grab,1)
          
          % Initialize vectors for saving
          tmpArray = nan(labelCountTotal, binCount);
          
          for label_i = 1:labelCountTotal
            
            % Find the indices to index into the data
            indices2Use = indices2Grab{site_i, label_i};
            
            % Collect the data using the indicies
            data2Store = mean(binnedData.binned_data{MUA_ind(site_i)}(indices2Use,:));
            
            % Store it in this spot
            %       binnedDataArray{label_i}(nanInd:nanInd+size(data2Store, 1)-1, :) = data2Store;
            tmpArray(label_i,:) = data2Store;
            
          end
          
          % Normalize all the activity from this specific site, across
          % conditions.
          tmpArray = normalize(tmpArray);
          
          for label_i = 1:labelCountTotal
            % Find the next slot to put data into.
            nanInd = find(isnan(binnedDataArray{label_i}(:,1)),1);
            binnedDataArray{label_i}(nanInd:nanInd+size(data2Store, 1)-1, :) = tmpArray(label_i,:);
          end
          
        end
        
        % Remove the NaNs
        for grp_i = 1:length(binnedDataArray)
          nanInd = find(isnan(binnedDataArray{grp_i}(:,1)),1);
          binnedDataArray{grp_i} = binnedDataArray{grp_i}(1:nanInd-1, :);
        end
        
        % Combine across the arrays
        labelCount = size(binnedDataArray{1}, 1);
        
        % Dummy check
        %   A = ones(1, 72);
        %   B = ones(1, 72) * 2;
        %   binnedDataArrayD = {A; B; A; B};
        %
        %   imagesc(vertcat(binnedDataArrayD{:}))
        %   imagesc(horzcat(binnedDataArrayD{:}))
        
        binnedDataStack = horzcat(binnedDataArray{:});
        binnedLabelStack = repmat(1:length(binnedDataArray), [labelCount, 1]);
        binnedLabelStack = binnedLabelStack(:);
        binnedDataStack(isnan(binnedDataStack)) = 0;
        
        storeData{paradigm_i, unit_Type_i} = binnedDataStack';
        
        % pca time
        [coeff, score, latent, tsquared, explained, ~] = pca(binnedDataStack');
        
        figTitle = sprintf('Data input to PCA, blocked by type - %s - %s', paradigm{paradigm_i}, unitIndName{unit_Type_i});
        figure('name', figTitle)
        imagesc(binnedDataStack)
        hold on
        ylabel('Unit')
        xlabel('Bin')
        start = 1:binCount:size(score, 1);
        yLims = ylim();
        ylim(yLims);
        for ii = 2:length(start)
          plot([start(ii), start(ii)], yLims, 'color', 'k', 'linewidth', 2)
        end
        title(figTitle);
        saveFigure(outputDir, ['1. ' figTitle], [], figStruct, [])
        
        figTitle = sprintf('Score per Unit, blocked by type - %s - %s', paradigm{paradigm_i}, unitIndName{unit_Type_i});
        figure('name', figTitle)
        imagesc(score)
        hold on
        ylabel('Unit')
        xlabel('Bin')
        start = 1:binCount:size(score, 1);
        yLims = ylim();
        ylim(yLims);
        for ii = 2:length(start)
          plot([start(ii), start(ii)], yLims, 'color', 'k', 'linewidth', 2)
        end
        title(figTitle);
        saveFigure(outputDir, ['2. ' figTitle], [], figStruct, [])
        
        figTitle = sprintf('PC Coefficients for each Unit - %s - %s', paradigm{paradigm_i}, unitIndName{unit_Type_i});
        figure('name', figTitle)
        imagesc(coeff)
        hold on
        ylabel('Eigenvector')
        xlabel('Bin start time')
        xticks(bins_to_label);
        xticklabels(points_to_label)
        xlim([1, length(the_bin_start_times_shift)]);
        for ii = 1:length(x_for_lines)
          plot([x_for_lines(ii), x_for_lines(ii)], ylim(), 'linewidth', 4, 'color', [0.2, 0.2, 0.2])
        end
        title(figTitle);
        saveFigure(outputDir, ['3. ' figTitle], [], figStruct, [])
        
        figTitle = sprintf('Latent variables cumulative sum - %s - %s', paradigm{paradigm_i}, unitIndName{unit_Type_i});
        figure('name', figTitle)
        plot(cumsum(explained))
        hold on
        ylim([0 100])
        xlim([1 length(explained)])
        xlabel('Principal Component')
        ylabel('Cumulative Explained Variance')
        title(figTitle)
        for ii = 2:length(start)
          plot([start(ii), start(ii)], yLims, 'color', 'k', 'linewidth', 2)
        end
        saveFigure(outputDir, ['4. ' figTitle], [], figStruct, [])
        
        bhvLabels = labels2Include{1};
        paradigmLabel = paradigm{paradigm_i};
        save(fullfile(outputDir, sprintf('pcaData_%s_%s', paradigm{paradigm_i}, unitIndName{unit_Type_i})), 'coeff', 'score', 'latent', 'binCount', 'paradigmLabel', 'bhvLabels', 'binLabelInfo')
        
      end
    end
  end
end
  error('Done');
  
  % Combo PCA Stuff - needs to be improved to actually match up units with
  % each other across the paradigms and linked up to batchAnalysisLite.
  save('PCA_data_tmp', 'storeData')
  
  for unit_Type_i = [1, 3]
    
    % Combine Animated and non-Animated
    comboData = storeData(:, unit_Type_i);
    
    % Find whichever is shorter, use that
    comboLengths = cellfun('length', comboData);
    units2Take = min(comboLengths);
    
    comboDataTrimmed = cellfun(@(x) x(:, 1:units2Take), comboData, 'UniformOutput', false);
    
    comboDataStack = vertcat(comboDataTrimmed{:});
    
    % pca time
    [coeff, score, latent, tsquared, explained, ~] = pca(comboDataStack);
    
    figTitle = sprintf('Data input to PCA, blocked by type - Combo Paradigm - %s', unitIndName{unit_Type_i});
    figure('name', figTitle)
    imagesc(comboDataStack)
    hold on
    ylabel('Unit')
    xlabel('Bin')
    start = 1:binCount:size(score, 1);
    yLims = ylim();
    ylim(yLims);
    for ii = 2:length(start)
      plot([start(ii), start(ii)], yLims, 'color', 'k', 'linewidth', 2)
    end
    title(figTitle);
    saveFigure(outputDir, ['1. ' figTitle], [], figStruct, [])
    
    figTitle = sprintf('Score per Unit, blocked by type - Combo Paradigm - %s', unitIndName{unit_Type_i});
    figure('name', figTitle)
    imagesc(score)
    hold on
    ylabel('Unit')
    xlabel('Bin')
    start = 1:binCount:size(score, 1);
    yLims = ylim();
    ylim(yLims);
    for ii = 2:length(start)
      plot([start(ii), start(ii)], yLims, 'color', 'k', 'linewidth', 2)
    end
    title(figTitle);
    saveFigure(outputDir, ['2. ' figTitle], [], figStruct, [])
    
    figTitle = sprintf('PC Coefficients for each Unit - Combo Paradigm - %s', unitIndName{unit_Type_i});
    figure('name', figTitle)
    imagesc(coeff)
    hold on
    ylabel('Eigenvector')
    xlabel('Bin start time')
    xticks(bins_to_label);
    xticklabels(points_to_label)
    xlim([1, length(the_bin_start_times_shift)]);
    for ii = 1:length(x_for_lines)
      plot([x_for_lines(ii), x_for_lines(ii)], ylim(), 'linewidth', 4, 'color', [0.2, 0.2, 0.2])
    end
    title(figTitle);
    saveFigure(outputDir, ['3. ' figTitle], [], figStruct, [])
    
    figTitle = sprintf('Latent variables cumulative sum - Combo Paradigm - %s', unitIndName{unit_Type_i});
    figure('name', figTitle)
    plot(cumsum(explained))
    hold on
    ylim([0 100])
    xlim([1 length(explained)])
    xlabel('Principal Component')
    ylabel('Cumulative Explained Variance')
    title(figTitle)
    for ii = 2:length(start)
      plot([start(ii), start(ii)], yLims, 'color', 'k', 'linewidth', 2)
    end
    saveFigure(outputDir, ['4. ' figTitle], [], figStruct, [])
    
    bhvLabels = [labels2Include{1} labels2Include{1}];
    paradigmLabel = 'Combo';
    save(sprintf('pcaData_%s_%s', paradigmLabel, unitIndName{unit_Type_i}), 'coeff', 'score', 'latent', 'binCount', 'paradigmLabel', 'bhvLabels', 'binLabelInfo')
    
  end
  
end
function pcaPrep(spikePathBank, params)

% Switches
normalizeTraces = true;     % A switch which normalizes a trace, following the mean being taken.

% Identify monkeys in set
dateSubjClipped = extractBefore(spikePathBank.Row, '00');
nameVec = cellfun(@(x) x(10:end), dateSubjClipped, 'UniformOutput', false);
outputDir = params.pcaParams.outputDir;
preprocDir = params.pcaParams.preprocDir;

if ~exist(preprocDir, 'dir')
  mkdir(preprocDir)
end

% Identify present monkeys in dataset
nameArray = unique(nameVec);
paradigmVec = unique(spikePathBank.paradigmName);
plotIndParams.stimParamsFilename = params.stimParamsFilename;
dataFilesFound = params.NDTParams.rasterDir;
rasterDataDir = fullfile(outputDir, 'rasterData');

% Files to make
% Units 
unitTypes = {'M', 'U'};
unitTypesLabel = {'MUA', 'Units'};

% Labels to collect - the code grabs data belonging to each of these labels
% and takes a mean to create a block per unit. 
labels2Include = {{'chasing', 'fighting', 'grooming', 'mounting', 'goalDirected', 'idle', 'objects', 'scene'}}; % Label sorting
% labels2Include = {{'socialInteraction', 'goalDirected', 'idle', 'objects', 'scene'}}; % Label sorting
labelSet = {'categories'};

% storeData = cell(length(paradigmVec), length(unitTypes), length(nameArray));

for m_i = 1:length(nameArray)
  for paradigm_i = 1:length(paradigmVec)
    
    paradigmLabel = paradigmVec{paradigm_i};    % Stated like this due to need to save later.
    
    % Define path to data file and load it.
    rasterDataDirMP = fullfile(rasterDataDir, paradigmVec{paradigm_i}, nameArray{m_i});
    spikePathBankMP = spikePathBank(contains(spikePathBank.Row, nameArray{m_i}) & contains(spikePathBank.paradigmName, paradigmVec{paradigm_i}), :);
    
    % Check if Raster data is present, if not, generate.
    if ~exist(rasterDataDirMP, 'dir') | length(dir(rasterDataDirMP)) == 2
      paradigmParams = params.NDTParams.spikeToRasterParams;
      paradigmParams.spikePathLoadParams = params.spikePathLoadParams;
      paradigmParams.rasterParams = params.NDTParams.spikeToRasterParams.(paradigmLabel);
      spikePathBank_to_rasterData(spikePathBankMP, rasterDataDirMP, paradigmParams);
    end
    
    % Bin the data
    binnedDataPath = pcaBin(rasterDataDirMP, 100, 100, false);
    binnedData = load(binnedDataPath);
        
    % If desired, smooth the data
    smoothKernal = gausswin(5);
    smoothKernal = smoothKernal/sum(smoothKernal);
    parfor ii = 1:length(binnedData.binned_data)
      binData = binnedData.binned_data{ii};
      for jj = 1:size(binData,1)
        binData(jj,:) = conv(binData(jj,:), smoothKernal, 'same');
      end
      binned_data{ii} = binData;
    end
    
    % If desired, normalize the data
    if normalizeTraces
      parfor ii = 1:length(binnedData.binned_data)
        % Normalize the data for the site
        normData = binned_data{ii};
        dataMean = mean(normData(:), 'omitnan');
        dataSD = std(normData(:), 'omitnan');
        if dataSD ~= 0
          binned_data{ii} = (normData - dataMean)/dataSD;
        else
          continue
        end
      end
    end
    
    % Scramble labels - Done by hand.
    % unitInd = unitInd(randperm(length(unitInd)));
    binStarts = binnedData.binned_site_info.binning_parameters.the_bin_start_times;
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
    binCount = size(binned_data{1}, 2);
    
    % Grab Related labels
    stimPerSite = binnedData.binned_labels.stimuli;
    unitType = binnedData.binned_site_info.UnitType;
    
    for unit_Type_i = 1:length(unitTypes)
            
      % Identify data to extract, initialize structure
      unitInd = strcmp(unitType, unitTypes{unit_Type_i});
      binnedDataPerGroup = nan(sum(unitInd), binCount);

      for group_i = 1:length(labels2Include)
        
        % Final destination file name - if it exists, skip this step
        fileOutName = sprintf('pcaData_%s_%s_%s_%s.mat', labelSet{group_i}, nameArray{m_i}, paradigmLabel, unitTypesLabel{unit_Type_i});
        dataOutFile = fullfile(preprocDir, fileOutName);
        
        if exist(dataOutFile, 'file')
          fprintf('Found %s, skipping... \n', fileOutName)
          continue
        else
          fprintf('Generating %s... \n', fileOutName)
        end
        
        % Look for the indices to include in each group.
        plotIndParams.plotLabels = labels2Include(group_i);
        labelCountTotal = length(labels2Include{group_i});

        % Initialize array for storing data.
        binnedDataArray = cell(length(labels2Include{group_i}),1);
        [binnedDataArray{:}] = deal(binnedDataPerGroup);
        
        % Collect the indices which will be used to collect data.
        unit_index = find(unitInd)';

        for site_i = 1:length(unit_index)
          
          % Identify trials to be collected from the binned_data.
          stimuliIndPerLabel = plotIndex(stimPerSite{unit_index(site_i)}, plotIndParams);
          [A, ~, C] = unique(stimuliIndPerLabel);
          a_counts = accumarray(C, 1); % 1 is 0, and shouldn't be counted.
          counts2Keep = find(A ~= 0)';
          a_counts = a_counts(counts2Keep);
                              
          for label_i = counts2Keep

            % Identify trials, collect their indicies
            label_inds = find(C == label_i);
            
            % If there are too many trials, grab a random set of trials to reach the minimum
            if length(label_inds) > min(a_counts)
              label_inds = label_inds(randperm(length(label_inds), min(a_counts)));
            end
            
            % Grab the data
            binnedDataArray{A(label_i)}(site_i, :) = mean(binned_data{unit_index(site_i)}(label_inds,:));
                        
          end
                              
        end

        binnedDataStack = horzcat(binnedDataArray{:});
        binnedDataStack(isnan(binnedDataStack)) = 0;
        
        % storeData{paradigm_i, unit_Type_i, m_i} = binnedDataStack';
        
        % pca time
        [coeff, score, latent, tsquared, explained, ~] = pca(binnedDataStack');
        
        figTitle = sprintf('Data input to PCA, blocked by type - %s - %s', paradigmVec{paradigm_i}, unitTypesLabel{unit_Type_i});
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
        saveFigure(outputDir, ['1. ' figTitle], [], params.figStruct, [])
        
        figTitle = sprintf('Score per Unit, blocked by type - %s - %s', paradigmVec{paradigm_i}, unitTypesLabel{unit_Type_i});
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
        saveFigure(outputDir, ['2. ' figTitle], [], params.figStruct, [])
        
        figTitle = sprintf('PC Coefficients for each Unit - %s - %s', paradigmVec{paradigm_i}, unitTypesLabel{unit_Type_i});
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
        saveFigure(outputDir, ['3. ' figTitle], [], params.figStruct, [])
        
        figTitle = sprintf('Latent variables cumulative sum - %s - %s', paradigmVec{paradigm_i}, unitTypesLabel{unit_Type_i});
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
        saveFigure(outputDir, ['4. ' figTitle], [], params.figStruct, [])
        
        bhvLabels = labels2Include{1};
        save(dataOutFile, 'coeff', 'score', 'latent', 'binCount', 'paradigmLabel', 'bhvLabels', 'binLabelInfo')
        
      end
    end
  end
end
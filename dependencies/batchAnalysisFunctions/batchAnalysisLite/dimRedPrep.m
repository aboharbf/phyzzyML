function dimRedPrep(spikePathBank, params)

% Switches
normalizeTraces = true;     % A switch which normalizes a trace, following the mean being taken.
smoothTraces = true;        % A switch that smooths a trace.
HzThreshold = true;         % A switch which removes units with low firing rate
HzThresNum = 1;             % The threshold used.
binSize = 100;
binStep = 100;

% Identify monkeys in set
dateSubjClipped = extractBefore(spikePathBank.Row, '00');
nameVec = cellfun(@(x) x(10:end), dateSubjClipped, 'UniformOutput', false);
outputDir = params.dimRedParams.outputDir;
preprocDir = params.dimRedParams.preprocDir;

if ~exist(preprocDir, 'dir')
  mkdir(preprocDir)
end

% Identify present monkeys in dataset
nameArray = [unique(nameVec); 'Combo'];
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
labels2Include = {{'chasing', 'fighting', 'grooming', 'mounting', 'goalDirected', 'idle', 'objects', 'scene'}, {'socialInteraction', 'goalDirected', 'idle', 'objects', 'scene'}}; % Label sorting
% labels2Include = {{'socialInteraction', 'goalDirected', 'idle', 'objects', 'scene'}}; % Label sorting
labelSet = {'categories', 'condensedCategories'};

% storeData = cell(length(paradigmVec), length(unitTypes), length(nameArray));

for m_i = 1:length(nameArray)
  for paradigm_i = 1:length(paradigmVec)
    
    paradigmLabel = paradigmVec{paradigm_i};    % Stated like this due to need to save later.
    
    % Define path to data file and load it.
    rasterDataDirMP = fullfile(rasterDataDir, paradigmVec{paradigm_i}, nameArray{m_i});
    if ~strcmp(nameArray{m_i}, 'Combo')
      spikePathBankMP = spikePathBank(contains(spikePathBank.Row, nameArray{m_i}) & contains(spikePathBank.paradigmName, paradigmVec{paradigm_i}), :);
    else
      spikePathBankMP = spikePathBank(contains(spikePathBank.paradigmName, paradigmVec{paradigm_i}), :);
    end
    
    % Check if Raster data is present, if not, generate.
    if ~exist(rasterDataDirMP, 'dir') || length(dir(rasterDataDirMP)) == 2
      paradigmParams = params.NDTParams.spikeToRasterParams;
      paradigmParams.spikePathLoadParams = params.spikePathLoadParams;
      paradigmParams.rasterParams = params.NDTParams.spikeToRasterParams.(paradigmLabel);
      paradigmParams.fixShorten = 0;
      spikePathBank_to_rasterData(spikePathBankMP, rasterDataDirMP, paradigmParams);
    end
    
    % Bin the data
    binnedDataPath = pcaBin(rasterDataDirMP, binSize, binStep, false);
    binnedData = load(binnedDataPath);
    minRate = round((binnedData.binned_site_info.binning_parameters.end_time/1000) * HzThresNum, 1);
    
    % Process switches below on smoothing, normalizing, or thresholding.
    if smoothTraces || normalizeTraces || HzThreshold
      preProcTagAll = 'SNT';
      smoothKernal = gausswin(5);
      smoothKernal = smoothKernal/sum(smoothKernal);
      excludeInd = false(length(binnedData.binned_data),1);
      binned_data = cell(size(binnedData.binned_data));
      for ii = 1:length(binnedData.binned_data)
        binData = binnedData.binned_data{ii};
        
        if HzThreshold
          % Check that the unit meets the desired threshold
          if mean(sum(binData,2)) < minRate
            excludeInd(ii) = true;
            continue
          end
        end
        
        % Normalize
        if normalizeTraces
          dataMean = mean(binData(:), 'omitnan');
          dataSD = std(binData(:), 'omitnan');
          binData = (binData - dataMean)/dataSD;
        end
        
        if smoothTraces
          for jj = 1:size(binData,1)
            % Smoooth
            binData(jj,:) = conv(binData(jj,:), smoothKernal, 'same');
          end
          
          binned_data{ii} = binData;
        end
        
        
      end
      
      binned_data = binned_data(~excludeInd);
      stimPerSite = binnedData.binned_labels.stimuli(~excludeInd);
      unitTypeVec = binnedData.binned_site_info.UnitType(~excludeInd);
      preProcTag = preProcTagAll([smoothTraces normalizeTraces HzThreshold]);

    else
      
      % Grab Related labels, process accordingly
      preProcTag = '';
      binned_data = binnedData.binned_data;
      stimPerSite = binnedData.binned_labels.stimuli;
      unitTypeVec = binnedData.binned_site_info.UnitType;
      
    end

    
    % Scramble labels - Done by hand.
    % unitInd = unitInd(randperm(length(unitInd)));
    
    binStarts = binnedData.binned_site_info.binning_parameters.the_bin_start_times - 1;
    preStim = unique([binnedData.binned_site_info.alignment_event_time{:}]);

    points_to_label = [-300, 0, 500, 1000, 1500, 2000, 2500, 3000];
    points_for_lines = [0, 2800];
    
    the_bin_start_times_shift = binStarts - preStim;
    bins_to_label = interp1(the_bin_start_times_shift, 1:length(the_bin_start_times_shift), points_to_label);
    x_for_lines = interp1(the_bin_start_times_shift, 1:length(the_bin_start_times_shift), points_for_lines);
        
    binLabelInfo.the_bin_start_times_shift = the_bin_start_times_shift;
    binLabelInfo.bins_to_label = bins_to_label;
    binLabelInfo.points_to_label = points_to_label;
    binLabelInfo.x_for_lines = x_for_lines;
    binLabelInfo.points_for_lines = points_for_lines;
    binCount = size(binned_data{1}, 2);
        
    for unit_Type_i = 1:length(unitTypes)
      
      unitType = unitTypes{unit_Type_i};
            
      % Identify data to extract, initialize structure
      unitInd = strcmp(unitTypeVec, unitTypes{unit_Type_i});
      binnedDataPerGroup = nan(sum(unitInd), binCount);

      for group_i = 1:length(labels2Include)
        
        % Final destination file name - if it exists, skip this step
        fileOutName = sprintf('dimRedData_%dW_%dS_%s_%s_%s_%s_%s.mat', binSize, binStep, labelSet{group_i}, nameArray{m_i}, paradigmLabel, unitTypesLabel{unit_Type_i}, preProcTag);
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
                
        labels = labels2Include{group_i};
        save(dataOutFile, 'binnedDataArray', 'binCount', 'paradigmLabel', 'labels', 'unitType', 'binLabelInfo')
                
      end
    end
  end
end
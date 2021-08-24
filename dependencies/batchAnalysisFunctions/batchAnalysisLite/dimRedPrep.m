function dimRedPrep(spikePathBank, params)

% Switches
switchStruct.normalizeTraces = true;     % A switch which normalizes a trace, following the mean being taken.
switchStruct.smoothTraces = true;        % A switch that smooths a trace.
switchStruct.smoothWindow = 100;        % A switch that smooths a trace.
switchStruct.HzThreshold = true;         % A switch which removes units with low firing rate
switchStruct.HzThresNum = 1;             % The threshold used.
binSize = 400;
binWidth = binSize;

points_to_label = [-200, 0, 500, 1000, 1500, 2000, 2500, 3000];
points_for_lines = [0, 2800];
binLabelInfo.points_to_label = points_to_label;
binLabelInfo.points_for_lines = points_for_lines;

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
unitTypes = {'U'};
unitTypesLabel = {'Units'};

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
    
    % Process and Bin the data
    [binnedDataPath, preProcTag] = preprocessAndBin(rasterDataDirMP, binWidth, binWidth, switchStruct);
    binnedData = load(binnedDataPath);
    
    % Scramble labels - Done by hand.
    % unitInd = unitInd(randperm(length(unitInd)));
    
    binStarts = binnedData.binned_site_info.binning_parameters.the_bin_start_times - 1;
    preStim = unique([binnedData.binned_site_info.alignment_event_time{:}]);
    
    the_bin_start_times_shift = binStarts - preStim;
    bins_to_label = interp1(the_bin_start_times_shift, 1:length(the_bin_start_times_shift), points_to_label);
    x_for_lines = interp1(the_bin_start_times_shift, 1:length(the_bin_start_times_shift), points_for_lines);
        
    binLabelInfo.the_bin_start_times_shift = the_bin_start_times_shift;
    binLabelInfo.bins_to_label = bins_to_label;

    binLabelInfo.x_for_lines = x_for_lines;
    binCount = size(binnedData.binned_data{1}, 2);
    unitTypeVec = binnedData.binned_site_info.UnitType;
    stimPerSite = binnedData.binned_labels.stimuli;
    binned_data = binnedData.binned_data;
        
    for unit_Type_i = 1:length(unitTypes)
      
      unitType = unitTypes{unit_Type_i};
            
      % Identify data to extract, initialize structure
      unitInd = strcmp(unitTypeVec, unitTypes{unit_Type_i});
      binnedDataPerGroup = nan(sum(unitInd), binCount);

      for group_i = 1:length(labels2Include)
        
        % Final destination file name - if it exists, skip this step
        fileOutName = sprintf('dimRedData_%dW_%dS_%s_%s_%s_%s_%s.mat', binSize, binWidth, labelSet{group_i}, nameArray{m_i}, paradigmLabel, unitTypesLabel{unit_Type_i}, preProcTag);
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
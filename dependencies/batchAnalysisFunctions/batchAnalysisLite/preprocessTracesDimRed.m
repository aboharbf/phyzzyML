function preprocessTracesDimRed(rasterDir, binSize, binStep, switchStruct)
% A function which works within dimRedPrep to applying smoothing,
% normalization, and thresholding, followed by binning.

% Find all the traces within the previous collected raster folder
% Identify output file name
binnedDataPath = fullfile(rasterDir, sprintf('binned_data_%dW_%dS.mat', binSize, binStep));

if 0%exist(binnedDataPath, 'file')
  fprintf('Data in %s already binned. Returning. \n', rasterDir)
  return
else
  fprintf('Binning data in %s... \n', rasterDir)
end

smoothTraces = switchStruct.smoothTraces;
smoothWindow = switchStruct.smoothWindow;
normalizeTraces = switchStruct.normalizeTraces;
HzThreshold = switchStruct.HzThreshold;

% Find the files in the directory
rasterFiles = dir(fullfile(rasterDir, 'S*.mat')); % Individual run data begins with an S.
rasterFiles = fullfile({rasterFiles.folder}, {rasterFiles.name})';

if smoothTraces
  smoothKernal = gausswin(smoothWindow);
  smoothKernal = smoothKernal/sum(smoothKernal);
end

% initialize structures
excludeInd = false(size(rasterFiles));
binned_data = cell(size(rasterFiles));
binned_labels = struct();
binned_site_info = struct();

for file_i = 1:length(rasterFiles)
  
  % print a message the the data is being binned
  curr_bin_string = [' Binning the data: ' num2str(file_i) ' of ' num2str(length(rasterFiles))];
  if file_i == 1
    disp(curr_bin_string);
  else
    fprintf([repmat(8,1,bin_str_len) curr_bin_string]);
  end
  bin_str_len = length(curr_bin_string);
  
  % Load file
  rasterFile = load(rasterFiles{file_i});
  
  % Grab the relevant labels
  if file_i == 1
    labels_to_get = fields(rasterFile.raster_labels);
    site_info_to_get = fields(rasterFile.raster_site_info);
    
    % Initialize arrays
    for label_i = 1:length(labels_to_get)
      binned_labels.(labels_to_get{label_i}) = cell(1, length(rasterFiles));
    end
    
    for site_info_i = 1:length(site_info_to_get)
      binned_site_info.(site_info_to_get{site_info_i}) = cell(1, length(rasterFiles));
    end
    
  end
  
  % Pull data
  binData = rasterFile.raster_data;
  
  % Check that the unit meets the desired threshold
  if HzThreshold
    if mean(sum(binData,2)) < minRate
      excludeInd(file_i) = true;
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
  end
  
  binned_data{file_i} = binData;
  
end

binned_data = binned_data(~excludeInd);
stimPerSite = binnedData.binned_labels.stimuli(~excludeInd);
unitTypeVec = binnedData.binned_site_info.UnitType(~excludeInd);

% Grab the raster data
[binned_data{file_i}, binStarts, binEnds] = binStuff(rasterFile.raster_data, binSize, binStep);

% Grab the labels and site info for the site
for label_i = 1:length(labels_to_get)
  binned_labels.(labels_to_get{label_i}){file_i} = rasterFile.raster_labels.(labels_to_get{label_i});
end

for site_info_i = 1:length(site_info_to_get)
  binned_site_info.(site_info_to_get{site_info_i}){file_i} = rasterFile.raster_site_info.(site_info_to_get{site_info_i});
end

% save extra information about the bin_width, sampling_interval, etc.
binned_site_info.binning_parameters.bin_width = binSize;
binned_site_info.binning_parameters.sampling_interval = binStep;
binned_site_info.binning_parameters.start_time = 1;
binned_site_info.binning_parameters.end_time = binEnds(end);
binned_site_info.binning_parameters.the_bin_start_times = binStarts;

disp(['Saving the binned data to the file: ' binnedDataPath])
save(binnedDataPath, 'binned_data', 'binned_labels', 'binned_site_info');

end
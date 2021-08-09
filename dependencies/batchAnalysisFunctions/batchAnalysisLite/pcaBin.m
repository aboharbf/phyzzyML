function binnedDataPath = pcaBin(rasterDir, binWidth, binStep, overwrite)
% a function which bins data in a directory to a single file.
% Inputs:
% - rasterDir: a path to a directory with binned files.
% - binWidth - the size of the window to bin across
% - binStep - the size of the step to take for every subsequent window.

% Outputs:
% binnedDataPath: a path to the binned file.
% The file has 3 fields
% binned_data - a 1*n cell array, each cell is the binned data per site
% binned_labels - a struct with many fields, each is 1*n cell array for the
% appropriate label of each trial.
% binned_site_info - a struct with fields session_ID and unit type, each a
% cell array of single entries with the label for the site.

% Identify output file name
binnedDataPath = fullfile(rasterDir, sprintf('binned_data_%dW_%dS.mat', binWidth, binStep));

if exist(binnedDataPath, 'file') && ~overwrite
  fprintf('Data in %s already binned. Returning. \n', rasterDir)
  return
else
  fprintf('Binning data in %s... \n', rasterDir)
end

% Find the files in the directory
rasterFiles = dir(fullfile(rasterDir, 'S*.mat')); % Individual run data begins with an S.
rasterFiles = fullfile({rasterFiles.folder}, {rasterFiles.name})';

% initialize structures
binned_data = cell(1, length(rasterFiles));
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
  
  % Grab the raster data 
  [binned_data{file_i}, binStarts, binEnds] = binStuff(rasterFile.raster_data, binWidth, binStep);
  
  % Grab the labels and site info for the site
  for label_i = 1:length(labels_to_get)
    binned_labels.(labels_to_get{label_i}){file_i} = rasterFile.raster_labels.(labels_to_get{label_i});
  end
  
  for site_info_i = 1:length(site_info_to_get)
    binned_site_info.(site_info_to_get{site_info_i}){file_i} = rasterFile.raster_site_info.(site_info_to_get{site_info_i});
  end
  
end

% save extra information about the bin_width, sampling_interval, etc.
binned_site_info.binning_parameters.bin_width = binWidth;
binned_site_info.binning_parameters.sampling_interval = binStep;
binned_site_info.binning_parameters.start_time = 1;
binned_site_info.binning_parameters.end_time = binEnds(end);
binned_site_info.binning_parameters.the_bin_start_times = binStarts;

disp(['Saving the binned data to the file: ' binnedDataPath])
save(binnedDataPath, 'binned_data', 'binned_labels', 'binned_site_info');



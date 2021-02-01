%% File clean up
% Function which cycles through all the NS5 files in a directory and:
% - Determines which channels have spikes using a generous threshold.
% - Saves a copy of the original file to another directory.
% - uses blackrock scripts to rearrange, and reassign channel values
% (probably should target solely at records in 2020,2021 to be in order.

%% Variables
dataDir = 'D:\EphysData\Data\**\2020*.ns5';    % directory w/ files to clean.
blackRockDir = 'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\dependencies\NPMK';    % directory w/ blackrock functions.
waveClusDir = 'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\dependencies\wave_clus';    % directory w/ blackrock functions.

outputDrive = 'E';                                  % Slower speed SSD, will hold all original data.
threshold = 2;                                      % STDs above
trace2Sample = 5;                                   % Time in minutes to sample for spikes.

% waveclus parameters
stdThreshold = 3;   % std over mean for something to be a spike. Be generous.




%% Script

% Add blackrock scripts.
addpath(genpath(blackRockDir));
addpath(genpath(waveClusDir));

% generate waveClus params
par = set_parameters();
par.stdmin = stdThreshold;
par.ref = double(par.ref_ms)*par.sr/1000; % Refractory time in samples.

% Calculate samples to check
samples2Check = par.sr * 60 * trace2Sample;

% Identify files
ns5Files2Clean = dir(dataDir);
ns5FilePaths = fullfile({ns5Files2Clean.folder}', {ns5Files2Clean.name}');
originalChannelNum = cell(size(ns5FilePaths));


for file_i = 1:length(ns5FilePaths)
  % Recover data for file
  tmp = openNSx(ns5FilePaths{file_i},'report','read');
  
  % Move the file to the other directory
  newPath = ns5FilePaths{file_i};
  newPath(1) = outputDrive;
  newDir = fileparts(newPath);
  
  if ~exist(newDir, 'dir')
    mkdir(newDir);
  end
  
  % copy the file w/ the new path
  copyfile(ns5FilePaths{file_i}, newPath);
  
  traceEnd = min(size(tmp.Data,2), par.sr * 60 * trace2Sample);
  
  % Identify channels to keep
  channels2Keep = false(size(tmp.Data, 1), 1);
  for chan_i = 1:length(channels2Keep)
    
    % Read the file and all its channels
    rawDataTrace = double(tmp.Data(chan_i,1:traceEnd));
    [spikes,~,~] = amp_detect(rawDataTrace, par);
    
    % If no spikes are found, label the channel to be kept.
    channels2Keep(chan_i) = ~isempty(spikes); 
    
  end
  
  % Save original channel numbers
  originalChannelNum{file_i} = find(channels2Keep);
  
  % Once the channels are all checked, see if any empty channels exist
  if ~all(channels2Keep)
    
    % Create a temporary file name
    [A, B, C] = fileparts(ns5FilePaths{file_i});
    tmpFile = [A, B, '_tmp', C];
    
    % Copy file to temporary new file
    saveChNSx_Simple(ns5FilePaths{file_i}, tmpFile, find(channels2Keep), 'reset')
    
    % Delete the original, and copy the temporary one with the new name,
    % then delete the temporary one.
    delete(ns5FilePaths{file_i});
    copyfile(tmpFile, ns5FilePaths{file_i});
    delete(tmpFile);
    
  end
  
end

disp('done');

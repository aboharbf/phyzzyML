function psidAndPlot(params)
% Implements PSID

psidPath = 'C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\dimRed\PSID';

addpath(genpath(psidPath))

% Identify files of interest
dimRedParams = dir(fullfile(params.dimRedParams.preprocDir, '*ST.mat'));
dimRedParams = fullfile({dimRedParams.folder}, {dimRedParams.name})';

% Define variables
outputDir = fullfile(params.dimRedParams.outputDir, 'PSID');
fileNames = extractBetween(dimRedParams, 'dimRedData_', '.mat'); % Specifically not normalized files.
testName = extractBefore(fileNames, '_');
titles = strrep(fileNames, '_', ' ');
colors = 'rbgkcmyk';
markerSize = 30;

for data_i = 1:length(dimRedParams)
  
  % Extract data
  tmp = load(dimRedParams{data_i});
  
  % Format correctly
  
  % Run function
  % y
  idSys = PSID(y, z, nx, n1, i);
  
  
  
  
end
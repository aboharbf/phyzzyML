%% Neural Decoding Toolbox Testing

% Add Phyzzy
addpath(genpath('C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML'))

params.NDTPath = 'C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\NeuralDecodingToolbox';
params.rasterDir = 'H:\Analyzed\batchAnalysis\NeuralDecodingTB\rasterData';
params.outputDir = 'H:\Analyzed\batchAnalysis\NeuralDecodingTB';
spikeDataDir = 'H:\Analyzed\batchAnalysis';


saveSpikeDataBank([], [], 'load', spikeDataDir);


% Add the path to the NB
addpath(genpath(params.NDTPath));

if ~exist(params.outputDir, 'dir')
  mkdir(params.outputDir)
end

if ~exist(params.rasterDir, 'dir')
  spikeDataBank_to_rasterData(spikeDataBank, params.rasterDir);
end

% Step 1 - generate the appropriate binned data
binWidth = 150;
stepSize = 50;
binnedFileIn = fullfile(params.rasterDir, 'rasterData_binned');
binnedFileOut = fullfile(params.rasterDir, sprintf('rasterData_binned_%dms_bins_%dms_sampled.mat', binWidth, stepSize));
rasterDataPath = fullfile(params.rasterDir, 'S20*');

if ~exist(binnedFileOut, 'file')
  binnedFileOut = create_binned_data_from_raster_data(rasterDataPath , binnedFileIn, binWidth, stepSize);
end

% Step 2 - generate the data source object
specific_binned_labels_names = 'stimuli'; % Which label in the binned data would you like to classify?
num_cv_splits = 5; % use 20 cross-validation splits (which means that 19 examples of each object are used for training and 1 example of each object is used for testing)

% create the basic datasource object
ds = basic_DS(binnedFileOut, specific_binned_labels_names,  num_cv_splits);
% Add in steps to ensure only sites with adequate repeats are used
[ds.sites_to_use, B, C, D] = find_sites_with_k_label_repetitions(ds.the_labels, num_cv_splits);

if isempty(ds.sites_to_use)
  error('No sites meet criteria above')
end

% Step 3 - Generate a classifier and data preprocessor objects compatible with that data source.
the_classifier = max_correlation_coefficient_CL;
the_feature_preprocessors{1} = zscore_normalize_FP; 


% Step 4 - Generate a cross validator
the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);


% Step 4.1 - define some parameters on newly generated object
the_cross_validator.num_resample_runs = 5;


% Step 5 - Run decoding analyses
DECODING_RESULTS = the_cross_validator.run_cv_decoding; 
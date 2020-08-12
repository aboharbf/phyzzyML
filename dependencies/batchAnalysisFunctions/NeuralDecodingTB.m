function NeuralDecodingTB(spikeDataBank, params)
% Function which implements the Neural Decoding Toolbox (1.0.4). Takes in
% the spikeDataBank, generates the appropriate structure for the toolbox if
% not already generated, and runs the analysis specified in the NBTDParams.

% Add the path to the NB
addpath(genpath(params.NDTPath));

if ~exist(params.outputDir, 'dir')
  mkdir(params.outputDir)
end

if ~exist(params.rasterDir, 'dir')
  spikeDataBank_to_rasterData(spikeDataBank, params.rasterDir, params.spikeToRasterParams);
end

% Generate the appropriate binned data
binnedFileName = fullfile(params.rasterDir, 'rasterData_binned');
rasterDataPath = fullfile(params.rasterDir, 'S20*');
binnedFileNameOut = sprintf('%s_%dms_bins_%dms_sampled.mat', binnedFileName, params.binWidth, params.stepSize);

if ~exist(binnedFileNameOut, 'file')
  binnedFileName = create_binned_data_from_raster_data(rasterDataPath, binnedFileName, params.binWidth, params.stepSize);
else
  binnedFileName = binnedFileNameOut;
end

% Grab a variable which will be useful later
tmp = dir(rasterDataPath);
tmp = tmp(1);
tmp = load(fullfile(tmp.folder, tmp.name));
params.end_time = size(tmp.raster_data, 2);

% Analyses specified by k_aid app in NDTB, stored in phyzzy folder.
analysesFiles = dir([params.NDTAnalysesPath, filesep, '*.mat']);

for ii = 1:length(analysesFiles)
  % Step 1 - Load the specified parameters
  tmp = load(fullfile(analysesFiles(ii).folder, analysesFiles(ii).name));
  tmp = tmp.analysisStruct;
  fieldsToReplace = fieldnames(tmp);
  
  % Load
  templateAnalysisParams = params.AnalysesDefault;
  for field_i = 1:length(fieldsToReplace)
    templateAnalysisParams.(fieldsToReplace{field_i}) = tmp.(fieldsToReplace{field_i});    % Load in fields
  end
  
  if ~isfield(templateAnalysisParams, 'num_features_to_use')
    % Adding k to all for the sake of seeing p values for features.
    templateAnalysisParams.preProc = [templateAnalysisParams.preProc, {'select_or_exclude_top_k_features_FP'}];
    templateAnalysisParams.num_features_to_exclude = 0;
    templateAnalysisParams.num_features_to_use = length(templateAnalysisParams.sites);
  end
  [~, fileName, ~] = fileparts(analysesFiles(ii).name);
  templateAnalysisParams.plotTitle = tmp.plotTitle;
  templateAnalysisParams.load_data_as_spike_counts = strcmp(templateAnalysisParams.classifier, 'poisson_naive_bayes_CL');
  params.Analyses.(fileName) = templateAnalysisParams;
  
end

analysesToRun = fields(params.Analyses);
analysesStructs = params.Analyses;

for analysis_i = 4:length(analysesToRun)
  % Step 2 - generate the data source object
  analysisStruct = analysesStructs.(analysesToRun{analysis_i}); % Which label in the binned data would you like to classify?
  
  % create the basic datasource object
  ds = basic_DS(binnedFileName, analysisStruct.label,  analysisStruct.num_cv_splits, analysisStruct.load_data_as_spike_counts);
  ds.num_times_to_repeat_each_label_per_cv_split = 1;
  
  % Add in steps to ensure only sites with adequate repeats are used
  ds.label_names_to_use = analysisStruct.label_names_to_use;
  ds.sites_to_use = analysisStruct.sites;

  % Step 3 - Generate a classifier and data preprocessor objects compatible with that data source.
  the_classifier = eval(analysisStruct.classifier);

  the_feature_preprocessors = cell(length(analysisStruct.preProc), 1);
  for pp_i = 1:length(analysisStruct.preProc)
    the_feature_preprocessors{pp_i} = eval(analysisStruct.preProc{pp_i});
    
    switch analysisStruct.preProc{pp_i}
      case 'select_or_exclude_top_k_features_FP'
        the_feature_preprocessors{pp_i}.num_features_to_use = analysisStruct.num_features_to_use;
        the_feature_preprocessors{pp_i}.num_features_to_exclude = analysisStruct.num_features_to_exclude;
        the_feature_preprocessors{pp_i}.save_extra_info = analysisStruct.save_extra_preprocessing_info;
      case 'select_pvalue_significant_features_FP'
        the_feature_preprocessors{pp_i}.pvalue_threshold = analysisStruct.pvalue_threshold;
        the_feature_preprocessors{pp_i}.save_extra_info = analysisStruct.save_extra_preprocessing_info;
    end
    
  end
    
  % Step 4 - Generate a cross validator, define some parameters on newly generated object
  the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
  the_cross_validator.num_resample_runs = analysisStruct.cross_validator_num_resample;

  % Step 5 - Run decoding analyses
  decoding_results = the_cross_validator.run_cv_decoding;
  
  % Step 6 - save the results
  [result_names{1}, save_file_name] = deal(fullfile(params.outputDir, sprintf('results_%s', analysesToRun{analysis_i})));
  save(save_file_name, 'decoding_results');
    
  % Step 7 - Plotting - 
  % Per Label Accuracy Trace, Figure 1
  figh = plot_per_label_accuracy(decoding_results, ds, analysisStruct);
  saveFigure(params.outputDir, ['1. ' figTitle], analysisStruct, params.figStruct, [])
  
  % create an object to plot the results, put line for stim start
  plot_obj = plot_standard_results_object(result_names);
  plot_obj.significant_event_times = 0;
  
  % optional argument, can plot different types of results
  %plot_obj.result_type_to_plot = 2;  % for example, setting this to 2 plots the normalized rank results
    
  % TCT Matrix, Figure 2
  plot_obj = plot_standard_results_TCT_object(save_file_name);
  plot_obj.significant_event_times = 0;   % the time when the stimulus was shown
  plot_obj.display_TCT_movie = 0;
  plot_obj.TCT_figure_number = length(findobj('type','figure')) + 1; % Open a new figure.
  plot_obj.plot_results;  % plot the TCT matrix and a movie showing if information is coded by a dynamic population code
  
  % Clean up the figure before saving
  j = gcf;
  j.Units = 'Normalized';
  j.Position = [0 0.05 0.5 0.6];
  figTitle = sprintf('Cross Classification Matrix of %s', analysisStruct.plotTitle);
  title(figTitle);
  saveFigure(params.outputDir, ['2. ' figTitle], analysisStruct, params.figStruct, [])
end

save(fullfile(params.outputDir, 'ndtParams.mat'), 'params')

end
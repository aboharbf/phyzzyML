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

for analysis_i = 1:length(analysesToRun)
  % Step 2 - generate the data source object
  analysisStruct = analysesStructs.(analysesToRun{analysis_i}); % Which label in the binned data would you like to classify?
  
  % create the basic datasource object
  % analysisStruct.load_data_as_spike_counts = 1; % FOR TESTING PNB CLASSIFIER
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
  
  % Create Paths
  save_file_dir = deal(fullfile(params.outputDir, sprintf('results_%s', analysesToRun{analysis_i})));
  if ~exist(save_file_dir, 'dir')
    mkdir(save_file_dir);
  end
  
  shuff_file_dir = fullfile(save_file_dir, 'shuff_dir');
  if ~exist(shuff_file_dir, 'dir')
    mkdir(shuff_file_dir);
    skip_decoding_load_instead = 0;
  else
    skip_decoding_load_instead = 1;
  end
    
  if ~skip_decoding_load_instead
    
  % Step 4 - Generate a cross validator, define some parameters on newly generated object
  cross_val_resample = analysisStruct.cross_validator_num_resample;
  
  parfor shuff_ind = 1:analysisStruct.null_shuffle_count+1
    fprintf('Running Shuffle ind %d, @ %s \n', shuff_ind, datetime())
    tic
    
    tmpds = ds;
    
    if shuff_ind ~= 1
      tmpds.randomly_shuffle_labels_before_running = 1;  % randomly shuffled the labels before running
      fileName = sprintf('decoding_results_shuffle%d', shuff_ind - 1);
      save_file_dir_itr = shuff_file_dir;
      tmpds.initialized = 0;       % Without this, the shuffle code isn't seen and shuffling doesn't take place until you create a new ds from scratch.
    else
      fileName = sprintf('decoding_results');
      save_file_dir_itr = save_file_dir;
    end
    
    the_cross_validator = standard_resample_CV(tmpds, the_classifier, the_feature_preprocessors);
    the_cross_validator.num_resample_runs = cross_val_resample;
    the_cross_validator.display_progress.zero_one_loss = 0;
    the_cross_validator.display_progress.resample_run_time = 0;
    
    % Step 5 - Run decoding analyses
    decoding_results = the_cross_validator.run_cv_decoding;
    
    % Step 6 - save the results
    parforsave(fullfile(save_file_dir_itr, fileName), decoding_results)
    
    fprintf('Run complete after %d minutes \n', toc/60);
  end
        
  end
  
  % Step 7 - Significance testing
  % Significance testing
  saved_results_struct_name = 'decoding_results';
  true_decoder_name = fullfile(save_file_dir, saved_results_struct_name);
  null_dist_dir = strcat(shuff_file_dir, filesep);
  pval_obj = pvalue_object(true_decoder_name, null_dist_dir);
  pval_obj.collapse_all_times_when_estimating_pvals = 1;       
  pval_obj.saved_results_structure_name = saved_results_struct_name;
  
  % Plot Significance testing
  [p_values, null_dist, ~] = pval_obj.create_pvalues_from_nulldist_files;
  params.sig_bins = p_values < params.p_val_threshold;
  params.decoding_threshold = prctile(null_dist, 100 - (params.p_val_threshold * 100), 'all');
  
  % Load the normal decoding results
  save_file_name = fullfile(save_file_dir, saved_results_struct_name);
  tmp_decoding = load(save_file_name);
  decoding_results = tmp_decoding.decoding_results;
  
  % Step 8 - Plotting
  % Per Label Accuracy Trace, Figure 1
  figTitle = sprintf('Per Label accuracy trace for %s', analysisStruct.plotTitle);
  plot_per_label_accuracy(decoding_results, ds, analysisStruct, params);
  saveFigure(params.outputDir, ['1. ' figTitle], analysisStruct, params.figStruct, [])
    
  % TCT Matrix, Figure 2
  params.figTitle = sprintf('Cross Temporal Decoding of %s', analysisStruct.plotTitle);
  if params.addTCTSigShading
    sigStr = sprintf(', %s%% threshold', num2str(100 - (params.p_val_threshold * 100)));
    params.figTitle = horzcat([params.figTitle, sigStr]);
  end
  generate_TCT_plot(analysisStruct, save_file_name, saved_results_struct_name, params)
  saveFigure(params.outputDir, ['2. ' params.figTitle], analysisStruct, params.figStruct, [])
  
end
save(fullfile(params.outputDir, 'ndtParams.mat'), 'params')

end
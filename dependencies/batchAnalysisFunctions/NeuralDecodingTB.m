function NeuralDecodingTB(spikePathBank, params)
% Function which implements the Neural Decoding Toolbox (1.0.4). Takes in
% the spikeDataBank, generates the appropriate structure for the toolbox if
% not already generated, and runs the analysis specified in the NBTDParams.

% Server version - segments everything in such a way that parfor can be
% spread out across all the runs - pre-decoding analyses and post-analysis
% plotting specifically are distinct.

% IMPORTANT SWITCHES
analysisChangeParams.expandLabelPerSplit = 1;  % Divides the number of labels per split by 3.
analysisChangeParams.swap2libsvm = 0;          % Swaps whatever classifier is defined in the file to libsvm.
analysisChangeParams.dontZScoreFeatures = 0;   % Turns off the Z scoring
analysisChangeParams.reportFeaturepVal = 0;   % Turns off the Z scoring

% Add the path to the NB
addpath(genpath(params.NDTPath));

if ~exist(params.outputDir, 'dir')
  mkdir(params.outputDir)
end

params.spikeToRasterParams.spikePathLoadParams = params.spikePathLoadParams;
paradigmList = unique(spikePathBank.paradigmName);

% Generate the appropriate binned data
for paradigm_i = 1:length(paradigmList)
  pName = paradigmList{paradigm_i};
  spikePathBankParadigm = spikePathBank(strcmp(spikePathBank.paradigmName, pName), :);
  pFolder = fullfile(params.outputDir, pName);
  pRasterDir = fullfile(pFolder, 'rasterData');
  
  % Check if Raster data is present, if not, generate.
  if ~exist(pRasterDir, 'dir') || length(dir(pRasterDir)) == 2
    paradigmParams = params.spikeToRasterParams;
    paradigmParams.rasterParams = params.spikeToRasterParams.(pName);
    spikePathBank_to_rasterData(spikePathBankParadigm, pRasterDir, paradigmParams);
  end
  
  % Check if the raster data is present, if not, make it
  binnedDirName = fullfile(pRasterDir, 'rasterData_binned');
  rasterDataPath = fullfile(pRasterDir, 'S20*');
  binnedFileNameOut = sprintf('%s_%dms_bins_%dms_sampled.mat', binnedDirName, params.binWidth, params.stepSize);
  
  if ~exist(binnedFileNameOut, 'file')
    binnedFileName = create_binned_data_from_raster_data(rasterDataPath, binnedDirName, params.binWidth, params.stepSize);
  else
    binnedFileName = binnedFileNameOut;
  end

  % Once these raster files have been generated, use this script to
  % generate analyses.
  analysisDir = fullfile(params.NDTAnalysesPath, pName);
  k_aid_generate_analyses(binnedFileName, pName, analysisDir, pFolder, analysisChangeParams, params.AnalysesDefault)
  
  % collect the generated analyses and turn them into the struct format.
%   analysesStructs = NDTB_prepareAnalysisStruct(analysisDir, pFolder, analysisChangeParams);
  analysesStructPaths = dir(fullfile(analysisDir, '*.mat'));
  analysesStructPaths = fullfile({analysesStructPaths.folder}, {analysesStructPaths.name})';
  
  logCellArray = cell(size(analysesStructPaths));
  parfor analysis_i = 1:length(analysesStructPaths)
    % Get Analysis structure
    tic
    tmp = load(analysesStructPaths{analysis_i});
    analysisStruct = tmp.analysisStruct;
    
    % Check if the analysis was already performed.
    if exist(analysisStruct.decoding_results_file, 'file')
      logCellArray{analysis_i} = sprintf('Done Ind %d', analysis_i);
      continue
    end
    
    % Make sure other directories exist.
    if ~exist(analysisStruct.shuff_file_dir, 'dir')
      mkdir(analysisStruct.shuff_file_dir);
    end
    
    if ~exist(analysisStruct.save_file_dir, 'dir')
      mkdir(analysisStruct.save_file_dir);
    end
    
    fprintf('Running Analysis ind %d, @ %s \n', analysis_i, datetime())
    
    % Step 2 - generate the data source object
    ds = basic_DS(analysisStruct.binnedFile, analysisStruct.label,  analysisStruct.num_cv_splits, analysisStruct.load_data_as_spike_counts, analysisStruct.pseudoGenString);
    ds.label_names_to_use = analysisStruct.label_names_to_use;
    ds.num_times_to_repeat_each_label_per_cv_split = analysisStruct.num_times_to_repeat_each_label_per_cv_split;
    
    % Add in steps to ensure only sites with adequate repeats are used
    ds.sites_to_use = analysisStruct.sites;
    
    % Step 3 - Generate a classifier and data preprocessor objects compatible with that data source.
    the_classifier = analysisStruct.the_classifier;
    
    the_feature_preprocessors = analysisStruct.the_feature_preprocessors;
    for pp_i = 1:length(analysisStruct.preProc)
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
    % Cycle through the runs
    ds.randomly_shuffle_labels_before_running = analysisStruct.shuffle_ds;  % randomly shuffled the labels before running
    
    the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
    the_cross_validator.num_resample_runs = analysisStruct.cross_val_resample;
    the_cross_validator.test_only_at_training_times = ~analysisStruct.crossTempDecode; % For speeding up testing.
    
    % Only stop once results have converged.
    if analysisStruct.shuffle_ds
      the_cross_validator.stop_resample_runs_only_when_specfic_results_have_converged.zero_one_loss_results = [];
    else
      the_cross_validator.stop_resample_runs_only_when_specfic_results_have_converged.zero_one_loss_results = 0.1;
    end
    
    the_cross_validator.display_progress.zero_one_loss = 0;
    the_cross_validator.display_progress.resample_run_time = 0;
    
    % Step 5 - Run decoding analyses
    decoding_results = the_cross_validator.run_cv_decoding;
    
    % Step 6 - save the results
    parforsave(analysisStruct.decoding_results_file, decoding_results)
%     save(analysisStruct.decoding_results_file, 'decoding_results')
  
    fprintf('Run complete after %s minutes \n', num2str(round(toc/60, 1)));
    
  end

  % Step 7 - Significance testing
  % Significance testing
%   saved_results_struct_name = 'decoding_results';
%   true_decoder_name = fullfile(save_file_dir, 'decoding_results_run1');
%   
%   pval_obj = pvalue_object(true_decoder_name, strcat(analysisStruct.shuff_file_dir, filesep));
%   pval_obj.collapse_all_times_when_estimating_pvals = 1;
%   pval_obj.saved_results_structure_name = saved_results_struct_name;
%   
%   % Plot Significance testing
%   [p_values, null_dist, ~] = pval_obj.create_pvalues_from_nulldist_files;
%   params.sig_bins = p_values < params.p_val_threshold;
%   params.decoding_threshold = prctile(null_dist, 100 - (params.p_val_threshold * 100), 'all');
%   
%   % Load the normal decoding results
%   save_file_name = dir(fullfile(save_file_dir, 'decoding_results*'));
%   save_file_name = fullfile({save_file_name.folder}, {save_file_name.name})';
%   tmp_decoding = load(save_file_name{1});
%   decoding_results = tmp_decoding.decoding_results;
%   
%   % Step 8 - Plotting
%   % Assign params for correct plotting.
%   params.plotParams = params.(pName).plotParams;
%   
%   % Per Label Accuracy Trace, Figure 1
%   plot_per_label_accuracy(decoding_results, analysisStruct, params);
%   
%   % TCT Matrix, Figure 2
%   if analysisStruct.crossTempDecode
%     generate_TCT_plot(analysisStruct, save_file_name{1}, saved_results_struct_name, params)
%   end
  
end

save(fullfile(params.outputDir, 'ndtParams.mat'), 'params')

end
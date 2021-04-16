function NeuralDecodingTBLite(spikePathBank, params)
% Function which implements the Neural Decoding Toolbox (1.0.4). Takes in
% the spikeDataBank, generates the appropriate structure for the toolbox if
% not already generated, and runs the analysis specified in the NBTDParams.

% IMPORTANT SWITCHES 
expandLabelPerSplit = 1;  % Divides the number of labels per split by 3.
swap2libsvm = 0;          % Swaps whatever classifier is defined in the file to libsvm.
dontZScoreFeatures = 0;   % Turns off the Z scoring

% Add the path to the NB
addpath(genpath(params.NDTPath));

if ~exist(params.outputDir, 'dir')
  mkdir(params.outputDir)
end

params.spikeToRasterParams.spikePathLoadParams = params.spikePathLoadParams;
paradigmList = unique(spikePathBank.paradigmName);

% Generate the appropriate binned data
for paradigm_i = 2:length(paradigmList)
  pName = paradigmList{paradigm_i};
  spikePathBankParadigm = spikePathBank(strcmp(spikePathBank.paradigmName, pName), :);
  pFolder = fullfile(params.outputDir, pName);
  pRasterDir = fullfile(pFolder, 'rasterData');
  
  if ~exist(pFolder, 'dir')
    mkdir(pFolder);
  end
  
  % Check if Raster data is present, if not, generate.
  if ~exist(pRasterDir, 'dir')
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
  
  % Grab a variable which will be useful later
  tmp = dir(rasterDataPath);
  tmp = tmp(1);
  tmp = load(fullfile(tmp.folder, tmp.name));
  params.end_time = size(tmp.raster_data, 2);
  
  % Analyses specified by k_aid app in NDTB, stored in phyzzy folder.
  analysesFiles = dir(fullfile(params.NDTAnalysesPath, pName, '*.mat'));
  % clear anything there before
  params.Analyses = struct();
  
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
    
    if swap2libsvm
      templateAnalysisParams.classifier = 'libsvm_CL';
    end
    
    if dontZScoreFeatures
      templateAnalysisParams.preProc = [];
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
    
    % For each analysis, make the folder to which output images will be
    % saved.
    if ~isempty(templateAnalysisParams.preProc)
      analysisSubType = [templateAnalysisParams.classifier(1:3) '_' templateAnalysisParams.preProc{1}(1:3)];
    else
      analysisSubType = [[templateAnalysisParams.classifier(1:3)]];
    end
    
    templateAnalysisParams.plotOutDir = fullfile(pFolder, analysisSubType);
    
    params.Analyses.(fileName) = templateAnalysisParams;
    
  end
  
  analysesToRun = fields(params.Analyses);
  analysesStructs = params.Analyses;
  
  for analysis_i = 1:length(analysesToRun)
    % Step 2 - generate the data source object
    analysisStruct = analysesStructs.(analysesToRun{analysis_i}); % Which label in the binned data would you like to classify?
    
    % TESTING 4/3/21
    if expandLabelPerSplit
      analysisStruct.num_cv_splits = floor(analysisStruct.num_cv_splits/3);
    end
    
    % create the basic datasource object
    % analysisStruct.load_data_as_spike_counts = 1; % FOR TESTING PNB CLASSIFIER
    ds = basic_DS(binnedFileName, analysisStruct.label,  analysisStruct.num_cv_splits, analysisStruct.load_data_as_spike_counts);
    if expandLabelPerSplit
      ds.num_times_to_repeat_each_label_per_cv_split = 3;
    else
      ds.num_times_to_repeat_each_label_per_cv_split = 1;
    end
    
    % Add in steps to ensure only sites with adequate repeats are used
    ds.label_names_to_use = analysisStruct.label_names_to_use;
    ds.sites_to_use = analysisStruct.sites;
    
    % Step 3 - Generate a classifier and data preprocessor objects compatible with that data source.
    the_classifier = eval(analysisStruct.classifier);
    
    if swap2libsvm
      the_classifier.kernel = 'rbf';
      the_classifier.gaussian_gamma = 1;
    end
    
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
    save_file_dir = deal(fullfile(analysisStruct.plotOutDir, sprintf('results_%s', analysesToRun{analysis_i})));
    if ~exist(save_file_dir, 'dir')
      mkdir(save_file_dir);
    end
    
    shuff_file_dir = fullfile(save_file_dir, 'shuff_dir');
    if ~exist(shuff_file_dir, 'dir') || length(dir(shuff_file_dir)) == 2
      mkdir(shuff_file_dir);
      skip_decoding_load_instead = 0;
    else
      skip_decoding_load_instead = 1;
    end
    
    if ~skip_decoding_load_instead
      
      % Step 4 - Generate a cross validator, define some parameters on newly generated object
      cross_val_resample = analysisStruct.cross_validator_num_resample;
      runShuff = [false(1, 1); true(7,1)]';
      
      % Cycle through the runs
      parfor shuff_ind = 1:length(runShuff)
        fprintf('Running Shuffle ind %d, @ %s \n', shuff_ind, datetime())
        tic
        
        tmpds = ds;
        
        if runShuff(shuff_ind)
          tmpds.randomly_shuffle_labels_before_running = 1;  % randomly shuffled the labels before running
          fileName = sprintf('decoding_results_shuffle%d', shuff_ind);
          save_file_dir_itr = shuff_file_dir;
          tmpds.initialized = 0;       % Without this, the shuffle code isn't seen and shuffling doesn't take place until you create a new ds from scratch.
        else
          fileName = sprintf('decoding_results_run%d', shuff_ind);
          save_file_dir_itr = save_file_dir;
        end
        
        the_cross_validator = standard_resample_CV(tmpds, the_classifier, the_feature_preprocessors);
        the_cross_validator.num_resample_runs = cross_val_resample;
        the_cross_validator.test_only_at_training_times = 1; % For speeding up testing.
        
        % Only stop once results have converged.
        if runShuff(shuff_ind)
          the_cross_validator.stop_resample_runs_only_when_specfic_results_have_converged.zero_one_loss_results = [];
        else
          the_cross_validator.stop_resample_runs_only_when_specfic_results_have_converged.zero_one_loss_results = 0.1;
        end
        
        the_cross_validator.display_progress.zero_one_loss = 0;
        the_cross_validator.display_progress.resample_run_time = 0;
        
        % Step 5 - Run decoding analyses
        decoding_results = the_cross_validator.run_cv_decoding;
        
        % Step 6 - save the results
        parforsave(fullfile(save_file_dir_itr, fileName), decoding_results)
        
        fprintf('Run complete after %s minutes \n', num2str(round(toc/60, 1)));
      end
      
    end
    
    % Step 7 - Significance testing
    % Significance testing
    saved_results_file_name = 'decoding_results_run1';
    saved_results_struct_name = 'decoding_results';
    true_decoder_name = fullfile(save_file_dir, saved_results_file_name);
    null_dist_dir = strcat(shuff_file_dir, filesep);
    pval_obj = pvalue_object(true_decoder_name, null_dist_dir);
    pval_obj.collapse_all_times_when_estimating_pvals = 1;
    pval_obj.saved_results_structure_name = saved_results_struct_name;
    
    % Plot Significance testing
    [p_values, null_dist, ~] = pval_obj.create_pvalues_from_nulldist_files;
    params.sig_bins = p_values < params.p_val_threshold;
    params.decoding_threshold = prctile(null_dist, 100 - (params.p_val_threshold * 100), 'all');
    
    % Load the normal decoding results
    save_file_name = dir(fullfile(save_file_dir, 'decoding_results*'));
    save_file_name = fullfile({save_file_name.folder}, {save_file_name.name})';
    decoding_results = cell(size(save_file_name));
    for result_i = 1:length(save_file_name)
      tmp_decoding = load(save_file_name{result_i});
      decoding_results{result_i} = tmp_decoding.decoding_results;
    end
    
    % Step 8 - Plotting
    % Assign params for correct plotting.
    params.plotParams = params.(pName).plotParams;
    
    % Per Label Accuracy Trace, Figure 1
    plot_per_label_accuracy(decoding_results, ds, analysisStruct, params);
    
    % TCT Matrix, Figure 2
%     params.figTitle = sprintf('Cross Temporal Decoding of %s', analysisStruct.plotTitle);
%     if params.addTCTSigShading
%       sigStr = sprintf(', %s%% threshold', num2str(100 - (params.p_val_threshold * 100)));
%       params.figTitle = horzcat([params.figTitle, sigStr]);
%     end
%     generate_TCT_plot(analysisStruct, save_file_name{1}, saved_results_struct_name, params)
%     saveFigure(pFolder, ['2. ' params.figTitle], analysisStruct, params.figStruct, [])
    
  end
end
save(fullfile(params.outputDir, 'ndtParams.mat'), 'params')

end
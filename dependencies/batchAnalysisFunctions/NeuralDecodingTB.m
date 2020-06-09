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

% Step 1 - generate the appropriate binned data
binnedFileName = fullfile(params.rasterDir, 'rasterData_binned');
rasterDataPath = fullfile(params.rasterDir, 'S20*');
binnedFileNameOut = sprintf('%s_%dms_bins_%dms_sampled.mat', binnedFileName, params.binWidth, params.stepSize);

if ~exist(binnedFileNameOut, 'file')
  binnedFileName = create_binned_data_from_raster_data(rasterDataPath, binnedFileName, params.binWidth, params.stepSize);
else
  binnedFileName = binnedFileNameOut;
end

groupingType = {{'U', 'US'}, {'M'}};
groupingLabel = {'Unsorted and Units', 'Multiunit Activity'};

% Additional analyses specified by App.
analysesFiles = dir([params.NDTAnalysesPath, filesep, '*.mat']);

for ii = 1:length(analysesFiles)
  % Load the specified parameters
  [~, fileName, ~] = fileparts(analysesFiles(ii).name);
  tmp = load(fullfile(analysesFiles(ii).folder, analysesFiles(ii).name));
  tmp = tmp.analysisStruct;
  fieldsToReplace = fieldnames(tmp);
  
  % Load
  templateAnalysisParams = params.AnalysesDefault;
  for field_i = 1:length(fieldsToReplace)
    templateAnalysisParams.(fieldsToReplace{field_i}) = tmp.(fieldsToReplace{field_i});    % Load in fields
  end
  templateAnalysisParams.plotTitle = fileName;
  params.Analyses.(fileName) = templateAnalysisParams;
end

analysesToRun = fields(params.Analyses);
analysesStructs = params.Analyses;

for analysis_i = 4:length(analysesToRun)
  % Step 2 - generate the data source object
  analysisStruct = analysesStructs.(analysesToRun{analysis_i}); % Which label in the binned data would you like to classify?
  
  % create the basic datasource object
  ds = basic_DS(binnedFileName, analysisStruct.label,  analysisStruct.num_cv_splits);
  ds.num_times_to_repeat_each_label_per_cv_split = analysisStruct.repeat_each_label_per_split;
  
  min_repeats = analysisStruct.num_cv_splits * ds.num_times_to_repeat_each_label_per_cv_split;
  
  % Add in steps to ensure only sites with adequate repeats are used
  if isfield(analysisStruct, 'sites')
    available_sites = analysisStruct.sites;
    ds.label_names_to_use = analysisStruct.label_names_to_use;
  else
    [available_sites, ~, ~, ds.label_names_to_use] = find_sites_with_k_label_repetitions(ds.the_labels, min_repeats, analysisStruct.label_names_to_use);
  end
  
  for group_i = 1:length(groupingType)
    if ~isfield(analysisStruct, 'sites')  % Analyses not specified w/ the app need to be run per group_i.
      ds.sites_to_use = available_sites;    % Step helps with iterating across loops.
      ds.sites_to_use = find_sites_meeting_label(ds, {'UnitType', groupingType{group_i}});
    else
      ds.sites_to_use = available_sites;
      group_i = length(groupingType);
    end
    % Step 3 - Generate a classifier and data preprocessor objects compatible with that data source.
    the_classifier = {max_correlation_coefficient_CL, poisson_naive_bayes_CL, libsvm_CL};
    the_classifier = the_classifier{logical(analysisStruct.classifier)};
    
    the_feature_preprocessors = {zscore_normalize_FP, select_pvalue_significant_features_FP, select_or_exclude_top_k_features_FP};
    the_feature_preprocessors = the_feature_preprocessors(logical(analysisStruct.featPre));
    
    % Step 4 - Generate a cross validator, define some parameters on newly generated object
    the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
    the_cross_validator.num_resample_runs = analysisStruct.cross_validator_num_resample;
    
    % Step 5 - Run decoding analyses
    DECODING_RESULTS = the_cross_validator.run_cv_decoding;
    
    % Plotting
    % save the results
    save_file_name = fullfile(params.outputDir, sprintf('results_%s', analysesToRun{analysis_i}));
    save(save_file_name, 'DECODING_RESULTS');
    
    % which results should be plotted (only have one result to plot here)
    result_names{1} = save_file_name;
    
    % create an object to plot the results
    plot_obj = plot_standard_results_object(result_names);
    
    % put a line at the time when the stimulus was shown
    plot_obj.significant_event_times = 0;
    
    % optional argument, can plot different types of results
    %plot_obj.result_type_to_plot = 2;  % for example, setting this to 2 plots the normalized rank results
    
    figure()
    plot_obj.plot_results;   % actually plot the results
    title(sprintf('Classification of %s, based on %s', analysisStruct.plotTitle, groupingLabel{group_i}));
    
    plot_obj = plot_standard_results_TCT_object(save_file_name);
    plot_obj.significant_event_times = 0;   % the time when the stimulus was shown
    plot_obj.display_TCT_movie = 0;
    plot_obj.TCT_figure_number = length(findobj('type','figure')) + 1; % Open a new figure.
    
    plot_obj.plot_results;  % plot the TCT matrix and a movie showing if information is coded by a dynamic population code
    title(sprintf('Classification of %s, based on %s', analysisStruct.plotTitle, groupingLabel{group_i}));
  end
  
end

save(fullfile(params.outputDir, 'ndtParams.mat'), 'params')

end
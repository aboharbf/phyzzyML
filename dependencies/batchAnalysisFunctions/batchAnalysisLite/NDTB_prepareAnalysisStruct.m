function analysisStructAll = NDTB_prepareAnalysisStruct(analysisDir, pFolder, params)

  % Analyses specified by k_aid_generate_analyses in NDTB, stored in phyzzy folder.
  analysesFiles = dir(fullfile(analysisDir, '*.mat'));
  analysisFileNames =   strrep({analysesFiles.name}', '.mat', '');
  analysesFiles = fullfile({analysesFiles.folder}, {analysesFiles.name});
  
  % Initialize output
  analysisStructAll = struct();
  
  % Cycle through the files found
  for ii = 1:length(analysesFiles)
    
    % Step 1 - Load the specified parameters
    tmp = load(analysesFiles{ii});
    analysisStruct = tmp.analysisStruct;
    
    % Depedning on switches set (at the top of NeuralDecodingTBLite, change
    % the analysis.
    if params.swap2libsvm
      analysisStruct.classifier = 'libsvm_CL';
    end
    
    if params.dontZScoreFeatures
      analysisStruct.preProc = [];
    end
    
    if params.reportFeaturepVal
      % Adding k to all for the sake of seeing p values for features.
      analysisStruct.preProc = [analysisStruct.preProc, {'select_or_exclude_top_k_features_FP'}];
      analysisStruct.num_features_to_exclude = 0;
      analysisStruct.num_features_to_use = length(analysisStruct.sites);
    end
    
    if params.expandLabelPerSplit
      analysisStruct.num_cv_splits = floor(analysisStruct.num_cv_splits/3);
    end
    
    % Define the plot name    
    analysisStruct.load_data_as_spike_counts = strcmp(analysisStruct.classifier, 'poisson_naive_bayes_CL');
    
    % Define the folder where the outputs for this analysis will be saved.
    if ~isempty(analysisStruct.preProc)
      analysisSubType = [analysisStruct.classifier(1:3) '_' analysisStruct.preProc{1}(1:3)];
    else
      analysisSubType = analysisStruct.classifier(1:3);
    end
    
    analysisStruct.plotOutDir = fullfile(pFolder, analysisSubType);
    
    analysisStructAll.(analysisFileNames{ii}) = analysisStruct;
  end

end
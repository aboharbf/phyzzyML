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
    
    if params.expandLabelPerSplit
      analysisStruct.num_cv_splits = floor(analysisStruct.num_cv_splits/4);
    end
    
    % Define the folder where the outputs for this analysis will be saved.
    if ~isempty(analysisStruct.preProc)
      analysisSubType = [analysisStruct.classifier(1:3) '_' analysisStruct.preProc{1}(1:3)];
    else
      analysisSubType = analysisStruct.classifier(1:3);
    end
    
    analysisStruct.plotOutDir = fullfile(pFolder, analysisSubType, strrep(analysisStruct.coreTitle, ' ', '_'));
    
    analysisStructAll.(analysisFileNames{ii}) = analysisStruct;
  end

end
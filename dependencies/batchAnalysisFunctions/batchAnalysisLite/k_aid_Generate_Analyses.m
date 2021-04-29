% K aid repeated
% a function which updates old analyses

rasterFile = {'D:\DataAnalysis\batchAnalysis\NeuralDecodingTB_HandDefined\NaturalSocial\rasterData\rasterData_binned_150ms_bins_50ms_sampled.mat', 'D:\DataAnalysis\batchAnalysis\NeuralDecodingTB_HandDefined\headTurnCon\rasterData\rasterData_binned_150ms_bins_50ms_sampled.mat'};
analysisTemplate = {'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses\NaturalSocial_hand\NS_Categories_AllUnits.mat', 'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses\headTurnCon_hand\HTC_Categories_AllUnits.mat'};
paradigmName = {'NS', 'HTC'};
changeOldAnalyses = false;
analysisOutDir = 'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses';

% Generate the two analyses of interest
categoriesStruct.label = 'socialCat';
categoriesStruct.label_names_to_use = {'chasing'  'fighting'  'goalDirected'  'grooming'  'idle'  'mounting'  'objects'};
socVNonSocStruct.label = 'social';
socVNonSocStruct.label_names_to_use = {'agents'  'socialInteraction'};

% Combine them for an iteritable
coreAnalysisStructs = [categoriesStruct socVNonSocStruct];

% Values to cycle through;
analyisTag = {'AllUnits', 'noHT', 'onlyHT', 'socIntSelAny', 'noSocIntSelAny', 'socIntOnset', 'socIntstimPres', 'socIntReward', 'socIntSel_noHT', 'socIntSel_onlyHT'};
analyisTagVars = {'', {{'subSel_headTurn_all', 0}}, {{'subSel_headTurn_all', 1}}, {{'socVNonSocSel_any', 1}}, {{'socVNonSocSel_any', 0}}, {{'socVNonSocSel_stimOnset', 1}}, {{'socVNonSocSel_stimPres', 1}}, {{'socVNonSocSel_reward', 1}},...
  {{'socVNonSocSel_any', 1}, {'subSel_headTurn_all', 0}}, {{'socVNonSocSel_any', 1}, {'subSel_headTurn_all', 1}}};

% Load the raster data
for par_i = 1:length(paradigmName)
  
  % Load the available raster data.
  rasterData = load(rasterFile{par_i});
  binnedLabels = rasterData.binned_labels;
  binnedSiteInfo = rasterData.binned_site_info;
  
  % Turn the binnedSiteInfo into a big logical matrix
  varFields = fields(binnedSiteInfo);
  varFields = varFields(~contains(varFields, 'binning_parameters'));
  varFieldsUniqueEntries = cell(size(varFields));
  
  for field_i = 1:length(varFields)
    % Extract data
    fieldData = binnedSiteInfo.(varFields{field_i});
    
    if size(fieldData,1) == 1
      fieldData = fieldData';
    end
    
    % Pull unique values, since each is being turned into an index
    [fieldDataUnique, ~, fieldData2LabelInd] = unique(fieldData);
    varFieldsUniqueEntries{field_i} = fieldDataUnique;
    
    % Store in the larger matrix.
    if field_i == 1
      binnedSiteInfoMatrix = fieldData2LabelInd';
    else
      binnedSiteInfoMatrix = [binnedSiteInfoMatrix; fieldData2LabelInd'];
    end
    
  end
  
  % Load the previously generated analyses to use as templates.
  tmp = load(analysisTemplate{par_i});
  
  for core_i = 1:length(coreAnalysisStructs)
    
    analysisTemplateStruct = tmp.analysisStruct;
    coreAnalysisStruct = coreAnalysisStructs(core_i);
    
    % Cycle through predetermined differences
    sites2KeepTemplate = false(size(binnedSiteInfoMatrix));
    
    % Run the initial k check
    labelPerSite = binnedLabels.(coreAnalysisStruct.label);
    label_names_to_use = coreAnalysisStruct.label_names_to_use;
    [available_sites_tmp, min_num_repeats_all_sites, ~, ~] = find_sites_with_k_label_repetitions(labelPerSite, 1, label_names_to_use);
    
    for analysis_i = 1:length(analyisTag)
      
      % Copy the Template
      analysisStruct = analysisTemplateStruct;
      siteInfoSelected = analysisStruct.site_info_selected;
      analysisStruct.sites = available_sites_tmp;
      
      % Cycle through the variables defined in the analysisStruct
      sites2Keep = sites2KeepTemplate;
      
      % Iterate and populate the sites2KeepTemplate, defining sites which meet
      % the criteria.
      for site_info_i = 1:length(siteInfoSelected)
        sites2Keep(site_info_i,:) = ismember(binnedSiteInfoMatrix(site_info_i,:), find(siteInfoSelected{site_info_i}));
      end
      
      % Check to see if anything needs to be changed.
      variables2Change = analyisTagVars{analysis_i};
      if ~isempty(variables2Change)
        
        % Modify based on variables2Change
        for var_d = 1:length(variables2Change)
          varChangeName = variables2Change{var_d}{1};
          varChangeVal = variables2Change{var_d}{2};
          
          % Find which row to change
          rowInd = strcmp(varFields, varChangeName);
          row2Change = find(rowInd);
          
          % Find the index of the values you want to keep
          val2Keep = find(varFieldsUniqueEntries{row2Change} == varChangeVal);
          
          % Go back to the binnedSiteInfoMatrix, compare against the updated
          % values, and store it.
          sites2Keep(site_info_i,:) = ismember(binnedSiteInfoMatrix(site_info_i,:), val2Keep);
        end
        
      end
      
      % Keep units which check all the boxes, store in the sites field.
      sites2Keep = sum(sites2Keep,1) == size(sites2Keep,1);
      
      % Update counts
      min_num_repeats_all_sites(~sites2Keep) = 0;
      sites_to_use = find(min_num_repeats_all_sites);
      analysisStruct.sites = sites_to_use;
      
      % Find the maximum number of k available for all the units
      new_k = min(min_num_repeats_all_sites(sites2Keep));
      [analysisStruct.num_cv_splits, analysisStruct.k_repeats_needed] = deal(new_k);
      
      % Update the Title
      originalTitle = analysisStruct.plotTitle;
      
      unitsSelected = analysisStruct.site_info_selected{strcmp(analysisStruct.site_info_items, 'UnitType')};
      if sum(unitsSelected) == 2
        unitTag = 'U&US';
      elseif sum(unitsSelected) == 1
        unitTag = 'MUA';
      end
      
      % Process the string
      removeInd = strfind(originalTitle, '(');
      titleCore = originalTitle(1:removeInd(1)-1);
      unitSubSecTag = sprintf('(%s)', analyisTag{analysis_i});
      repsTag = sprintf('(%d Reps %d %s)', new_k, length(find(sites2Keep)), unitTag);
      analysisStruct.plotTitle = [titleCore unitSubSecTag repsTag];
      
      % Create the save string
      saveFileName = sprintf('%s_%s_%s', paradigmName{par_i}, coreAnalysisStructs(core_i).label, analyisTag{analysis_i});
      
      % Save structure
      save(fullfile(analysisOutDir, saveFileName), 'analysisStruct');
      
    end
  end
  
end

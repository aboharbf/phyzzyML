function k_aid_generate_analyses()
% K aid repeated
% a function which updates old analyses

% Below are per paradigm values
paradigmName = {'NS', 'HTC'};
rasterFile = {'D:\DataAnalysis\batchAnalysis\NeuralDecodingTB_HandDefined\NaturalSocial\rasterData\rasterData_binned_150ms_bins_50ms_sampled.mat', ...
  'D:\DataAnalysis\batchAnalysis\NeuralDecodingTB_HandDefined\headTurnCon\rasterData\rasterData_binned_150ms_bins_50ms_sampled.mat'};
analysisTemplate = {'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses\NaturalSocial_hand\NS_Categories_AllUnits.mat', ...
  'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses\headTurnCon_hand\HTC_Categories_AllUnits.mat'};
analysisOutDir = {'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses\NaturalSocial',...
  'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses\headTurnCon'};

% Each paradigm iterates across the below values
categoriesStruct.label = 'socialCat';
categoriesStruct.label_names_to_use = {'chasing'  'fighting'  'goalDirected'  'grooming'  'idle'  'mounting'  'objects'};
socVNonSocStruct.label = 'social';
socVNonSocStruct.label_names_to_use = {'agents'  'socialInteraction'};
coreAnalysisStructs = [categoriesStruct socVNonSocStruct];

% analysisTags are used to name the file, and put into the title of the
% figure.
analyisTag = {'AllUnits', 'noHT', 'onlyHT', 'socIntSelAny', 'noSocIntSelAny', ...
  'socIntOnset', 'socIntstimPres', 'socIntReward', ...
  'socIntSel_noHT', 'socIntSel_onlyHT'...
  };

% Matching set of variables to change, names must match what is defined in
% the raster file.
analyisTagVars = {'', {{'subSel_headTurn_all_selInd', 0}}, {{'subSel_headTurn_all_selInd', 1}}, {{'sVns_any_selInd', 1}}, {{'sVns_any_selInd', 0}}, ...
  {{'sVns_stimOnset_selInd', 1}}, {{'sVns_stimPres_selInd', 1}}, {{'sVns_reward_selInd', 1}},...
  {{'sVns_any_selInd', 1}, {'subSel_headTurn_all_selInd', 0}}, {{'sVns_any_selInd', 1}, {'subSel_headTurn_all_selInd', 1}}...
  };

error('Need to define new analyses for bCat stuff');

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
    
    % Copy over the label, label names.
    analysisTemplateStruct.label = coreAnalysisStructs(core_i).label;
    analysisTemplateStruct.label_names_to_use = coreAnalysisStructs(core_i).label_names_to_use;
        
    % Run the initial k check
    labelPerSite = binnedLabels.(analysisTemplateStruct.label);
    [~, min_num_repeats_all_sites, ~, ~] = find_sites_with_k_label_repetitions(labelPerSite, 1, analysisTemplateStruct.label_names_to_use);
    
    % Create the sites2keep logical array for the template
    siteInfoSelected = analysisTemplateStruct.site_info_selected;
    
    % Cycle through the variables defined in the analysisStruct
    sites2KeepTemplate = false(size(binnedSiteInfoMatrix));
    
    % Iterate and populate the sites2KeepTemplate, defining sites which meet
    % the criteria.
    for site_info_i = 1:length(siteInfoSelected)
      sites2KeepTemplate(site_info_i,:) = ismember(binnedSiteInfoMatrix(site_info_i,:), find(siteInfoSelected{site_info_i}));
    end
    
    for analysis_i = 1:length(analyisTag)
      % Copy the template
      analysisStruct = analysisTemplateStruct;
      sites2KeepArray = sites2KeepTemplate;

      % Check to see if anything needs to be changed.
      variables2Change = analyisTagVars{analysis_i};
      if ~isempty(variables2Change)
        
        % Modify based on variables2Change
        for var_d = 1:length(variables2Change)
          varChangeName = variables2Change{var_d}{1};
          varChangeVal = variables2Change{var_d}{2};
          
          % Find which row to change
          row2Change = find(strcmp(varFields, varChangeName));
          
          % Find the index of the values you want to keep
          val2Keep = find(varFieldsUniqueEntries{row2Change} == varChangeVal);
          
          % Go back to the binnedSiteInfoMatrix, compare against the updated
          % values, and store it.
          sites2KeepArray(row2Change,:) = ismember(binnedSiteInfoMatrix(row2Change,:), val2Keep);
        end
        
      end
      
      % Keep units which check all the boxes, store in the sites field.
      sites2Keep = all(sites2KeepArray,1);
      
      % Update sites to be used, and the repeats used.
      analysisStruct.sites = find(sites2Keep);
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
      saveFileName = sprintf('%s_%s_%s', paradigmName{par_i}, analysisStruct.label, analyisTag{analysis_i});
      
      % Save structure
      save(fullfile(analysisOutDir, saveFileName), 'analysisStruct');
      
    end
  end
  
end
end

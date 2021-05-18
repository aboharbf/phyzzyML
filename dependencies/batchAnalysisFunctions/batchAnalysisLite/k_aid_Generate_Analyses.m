function k_aid_generate_analyses(binnedFileName, paradigmName, analysisDir)
% a function which generates analyses.
paradigmOptions = {'NaturalSocial', 'headTurnCon'};
paradigmNameTag = {'NS', 'HTC'};

par_i = strcmp(paradigmOptions, paradigmName);

if ~any(par_i)
  error('paradigm in not specified. Please check list or update k_aid_generate_analyses')
end

rasterFile = binnedFileName;
paradigmNameTag = paradigmNameTag{par_i};

if ~exist(analysisDir, 'dir')
  mkdir(analysisDir)
else
  % No need to re-do these analyses.
%   return
end

% Each paradigm iterates across the below values
categoriesStruct.label = 'socialCat';
categoriesStruct.label_names_to_use = {'chasing'  'fighting'  'goalDirected'  'grooming'  'idle'  'mounting'  'objects'};
categoriesStruct.coreTitle = 'Categories';
socVNonSocStruct.label = 'social';
socVNonSocStruct.label_names_to_use = {'agents'  'socialInteraction'};
socVNonSocStruct.coreTitle = 'Social vs non-Social Agents';
broadCatSocStruct.label = 'catBroad';
broadCatSocStruct.label_names_to_use = {'objects', 'idle', 'goalDirected', 'socialInteraction'};
broadCatSocStruct.coreTitle = 'Broad Categories';
coreAnalysisStructs = [categoriesStruct socVNonSocStruct broadCatSocStruct];

% All Generated analyses have these parameters
coreAnalysisParams.num_cv_splits = 1;
coreAnalysisParams.k_repeats_needed = 1;
coreAnalysisParams.save_extra_preprocessing_info = 0;
coreAnalysisParams.editFieldEnabled = [1 0 0];
coreAnalysisParams.num_features_to_use = 0;
coreAnalysisParams.editFieldEnabled = [1 0 0];

% All Generate analyses are generated with a combination of the features
% below.
decoders = {'max_correlation_coefficient_CL'};
decoderTags = {'MCC'};
preProcArray = {'zscore_normalize_FP'};
unitSets = {'MUA', 'U&US'};

% analysisTags are used to name the file, and put into the title of the
% figure.
analyisTag = {'AllUnits', 'noHT', 'onlyHT', 'socIntSelAny', 'noSocIntSelAny', ...
  'socIntEarly', 'socIntstimLate', 'socIntReward', ...
  'socIntSel_noHT', 'socIntSel_onlyHT'...
  'broadCatEarly', 'broadCatstimLate', 'broadCatReward', ...
  'broadCatSliding', 'broadCatAny', ...
  'broadCatAny_onlyHT', 'broadCatAny_noHT', 'nobroadCatAny_onlyHT'...
  'stimEarlySel', 'stimLateSel', 'stimRewardSel'...
  };

% Matching set of variables to change, names must match what is defined in
% the raster file.
analyisTagVars = {'', {{'subSel_headTurn_all_selInd', 0}}, {{'subSel_headTurn_all_selInd', 1}}, {{'sVns_any_selInd', 1}}, {{'sVns_any_selInd', 0}}, ...
  {{'sVns_stimEarly_selInd', 1}}, {{'sVns_stimLate_selInd', 1}}, {{'sVns_reward_selInd', 1}},...
  {{'sVns_any_selInd', 1}, {'subSel_headTurn_all_selInd', 0}}, {{'sVns_any_selInd', 1}, {'subSel_headTurn_all_selInd', 1}}...
  {{'bCat_stimEarly_selInd', 1}}, {{'bCat_stimLate_selInd', 1}}, {{'bCat_reward_selInd', 1}},...
  {{'SW_broadCatTest_selInd', 1}}, {{'bCat_any_selInd', 1}}, ...
  {{'bCat_any_selInd', 1}, {'subSel_headTurn_all_selInd', 1}}, {{'bCat_any_selInd', 1}, {'subSel_headTurn_all_selInd', 0}}, {{'bCat_any_selInd', 0}, {'subSel_headTurn_all_selInd', 1}}...
  {{'baseV_stimEarly_selInd', 1}}, {{'baseV_stimLate_selInd', 1}}, {{'baseV_reward_selInd', 1}}};

% Load the available raster data.
rasterData = load(rasterFile);
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

for unit_set_i = 1:length(unitSets)
  % Set the base state of the binnedSiteInfoMatrix
  
  % Create a template - all switches should be set to true, except for
  % UnitType - which should be either MUA or U&US
  sites2KeepTemplate = true(size(binnedSiteInfoMatrix));
  unitTypeInd = strcmp(varFields, 'UnitType');
  
  if strcmp(unitSets{unit_set_i}, 'MUA')
    sites2KeepTemplate(unitTypeInd, :) = ismember(binnedSiteInfoMatrix(unitTypeInd,:), 1);
    unitTag = 'MUA';
  else
    sites2KeepTemplate(unitTypeInd, :) = ismember(binnedSiteInfoMatrix(unitTypeInd,:), [2 3]);
    unitTag = 'UnUS';
  end
  
  for core_i = 1:length(coreAnalysisStructs)
    
    analysisTemplateStruct = struct();
    
    % Copy over the label, label names.
    analysisTemplateStruct.label = coreAnalysisStructs(core_i).label;
    analysisTemplateStruct.label_names_to_use = coreAnalysisStructs(core_i).label_names_to_use;
    analysisTemplateStruct.save_extra_preprocessing_info = 0;
    analysisTemplateStruct.num_features_to_use = 0;
    
    analysisTemplateStruct.coreTitle = sprintf('%s %s', paradigmNameTag, coreAnalysisStructs(core_i).coreTitle);
    
    % Run the initial k check
    labelPerSite = binnedLabels.(analysisTemplateStruct.label);
    [~, min_num_repeats_all_sites, ~, ~] = find_sites_with_k_label_repetitions(labelPerSite, 1, analysisTemplateStruct.label_names_to_use);
    
    for decoder_i = 1:length(decoders)
      for preProc_i = 1:length(preProcArray)
        for analysis_i = 1:length(analyisTag)
          % Copy the template
          analysisStruct = analysisTemplateStruct;
          analysisStruct.classifier = decoders{decoder_i};
          analysisStruct.preProc = preProcArray(preProc_i);
          
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
          
          % Generate the Title
          unitSubSecTag = sprintf('(%s)', analyisTag{analysis_i});
          repsTag = sprintf('(%d Reps %d %s)', new_k, length(find(sites2Keep)), unitTag);
          analysisStruct.plotTitle = strjoin({analysisStruct.coreTitle unitSubSecTag repsTag}, ' ');
          analysisStruct.plotTitle = strrep(analysisStruct.plotTitle, '_', ' ');
          
          % Create the filename.
          saveFileName = sprintf('%s_%s_%s_%s_%s', paradigmNameTag, analysisStruct.label, analyisTag{analysis_i}, unitTag, decoderTags{decoder_i});
          
          % Save structure
          save(fullfile(analysisDir, saveFileName), 'analysisStruct');
          
        end
      end
    end
  end
end

end
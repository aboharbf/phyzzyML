function k_aid_generate_analyses()
% K aid repeated
% a function which updates old analyses

% Below are per paradigm values
paradigmName = {'NS', 'HTC'};
rasterFile = {'D:\DataAnalysis\batchAnalysis\NeuralDecodingTB\NaturalSocial\rasterData\rasterData_binned_150ms_bins_50ms_sampled.mat', ...
  'D:\DataAnalysis\batchAnalysis\NeuralDecodingTB\headTurnCon\rasterData\rasterData_binned_150ms_bins_50ms_sampled.mat'};

analysisTemplate = {'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses\templateAnalysisNS.mat', ...
  'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses\templateAnalysisHTC.mat'};

analysisOutDir = {'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses\NaturalSocial',...
  'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses\headTurnCon'};

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

coreAnalysisParams.num_cv_splits = 1;
coreAnalysisParams.k_repeats_needed = 1;
coreAnalysisParams.save_extra_preprocessing_info = 0;
coreAnalysisParams.editFieldEnabled = [1 0 0];
coreAnalysisParams.num_features_to_use = 0;
coreAnalysisParams.editFieldEnabled = [1 0 0];

decoders = {'max_correlation_coefficient_CL'};
decoderTags = {'MCC'};
preProcArray = [{'zscore_normalize_FP'}; {'zscore_normalize_FP'}];
unitSets = {'MUA', 'U&US'};

% analysisTags are used to name the file, and put into the title of the
% figure.
analyisTag = {'AllUnits', 'noHT', 'onlyHT', 'socIntSelAny', 'noSocIntSelAny', ...
  'socIntOnset', 'socIntstimPres', 'socIntReward', ...
  'socIntSel_noHT', 'socIntSel_onlyHT'...
  'broadCatOnset', 'broadCatstimPres', 'broadCatReward', ...
  'broadCatSliding', ...
  };

% Matching set of variables to change, names must match what is defined in
% the raster file.
analyisTagVars = {'', {{'subSel_headTurn_all_selInd', 0}}, {{'subSel_headTurn_all_selInd', 1}}, {{'sVns_any_selInd', 1}}, {{'sVns_any_selInd', 0}}, ...
  {{'sVns_stimOnset_selInd', 1}}, {{'sVns_stimPres_selInd', 1}}, {{'sVns_reward_selInd', 1}},...
  {{'sVns_any_selInd', 1}, {'subSel_headTurn_all_selInd', 0}}, {{'sVns_any_selInd', 1}, {'subSel_headTurn_all_selInd', 1}}...
  {{'bCat_stimOnset_selInd', 1}}, {{'bCat_stimPres_selInd', 1}}, {{'bCat_reward_selInd', 1}},...
  {{'slidingWinBroadCat_selInd', 1}}};

%     {'bCat_stimOnset_selInd' }
%     {'bCat_stimPres_selInd'  }
%     {'bCat_reward_selInd'    }
%     {'slidingWinBroadCat_selInd'}];

% Load the raster data
for par_i = 1:length(paradigmName)
  
  if ~exist(analysisOutDir{par_i}, 'dir')
    mkdir(analysisOutDir{par_i})
  end
  
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
      
      analysisTemplateStruct.coreTitle = sprintf('%s %s', paradigmName{par_i}, coreAnalysisStructs(core_i).coreTitle);
      
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
            saveFileName = sprintf('%s_%s_%s_%s_%s', paradigmName{par_i}, analysisStruct.label, analyisTag{analysis_i}, unitTag, decoderTags{decoder_i});
            
            % Save structure
            save(fullfile(analysisOutDir{par_i}, saveFileName), 'analysisStruct');
            
          end
        end
      end
    end
  end
  
end
end

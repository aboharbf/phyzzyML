function analysisStructArray = generateAnalysisStructs(binnedFileName, paradigmName, analysisDir, pFolder, analysisChangeParams, params)
% a function which generates analyses.
% paradigmOptions = {'NaturalSocial', 'headTurnCon'};
% paradigmNameTag = {'NS', 'HTC'};

paradigmOptions = {'NaturalSocial', 'headTurnCon'};
paradigmNameTag = {'NS', 'HTC'};

% Function which may be of use
% ff2n(8)

par_i = strcmp(paradigmOptions, paradigmName);
subSampleRuns = false;             % A switch which enables the generation of a random subsampling of units for each analysis, in increments and numbers defined below.
unitsPerSubSample = [25:25:500 600:100:1000];
subsampleCount = 1;

if ~any(par_i)
  error('paradigm in not specified. Please check list or update k_aid_generate_analyses')
end

rasterFile = binnedFileName;
paradigmNameTag = paradigmNameTag{par_i};

if ~exist(analysisDir, 'dir')
  mkdir(analysisDir)
else
%   No need to re-do these analyses.
%   return
end

% Each paradigm iterates across the below values
categoriesStruct.labelFile = 'cat';
categoriesStruct.label = 'socialCat';
categoriesStruct.label_names_to_use = {'chasing'  'fighting'  'goalDirected'  'grooming'  'idle'  'mounting'  'objects' 'scene'};
categoriesStruct.coreTitle = 'Categories';

socCategoriesStruct.labelFile = 'socialCat';
socCategoriesStruct.label = 'socialCat';
socCategoriesStruct.label_names_to_use = {'chasing'  'fighting' 'grooming' 'mounting'};
socCategoriesStruct.coreTitle = 'Social Categories';

nonSocCategoriesStruct.labelFile = 'nonSocialCat';
nonSocCategoriesStruct.label = 'socialCat';
nonSocCategoriesStruct.label_names_to_use = {'goalDirected'  'idle'  'objects', 'scene'};
nonSocCategoriesStruct.coreTitle = 'Non-Social Categories';

socVNonSocStruct.labelFile = 'social';
socVNonSocStruct.label = 'social';
socVNonSocStruct.label_names_to_use = {'agents'  'socialInteraction'};
socVNonSocStruct.coreTitle = 'Social vs non-Social Agents';

broadCatSocStruct.labelFile = 'catBroad';
broadCatSocStruct.label = 'catBroad';
broadCatSocStruct.label_names_to_use = {'objects', 'idle', 'goalDirected', 'socialInteraction'};
broadCatSocStruct.coreTitle = 'Broad Categories';

% coreAnalysisStructs = [categoriesStruct socCategoriesStruct nonSocCategoriesStruct socVNonSocStruct broadCatSocStruct];
coreAnalysisStructs = categoriesStruct;
[coreAnalysisStructs.the_training_label_names] = deal([]);
[coreAnalysisStructs.the_test_label_names] = deal([]);
[coreAnalysisStructs.pseudoGenString] = deal('');
[coreAnalysisStructs.crossTempDecode] = deal(true);

if strcmp(paradigmName, 'NaturalSocial')
% PseudoGen related - create a version which does and does not generalize
% across instances of stimuli
coreAnalysisStructsPG = coreAnalysisStructs;
[coreAnalysisStructsPG.pseudoGenString] = deal('stimuli');
[coreAnalysisStructsPG.crossTempDecode] = deal(false);
[coreAnalysisStructs.pseudoGenString] = deal('');
[coreAnalysisStructs.crossTempDecode] = deal(true);
coreAnalysisStructs = [coreAnalysisStructs coreAnalysisStructsPG];

% New - Add for Stimuli, make sure pseudoGenString is empty
  stimuliSetAStruct.labelFile = 'stimSetA';
  stimuliSetAStruct.label = 'stimuli';
  stimuliSetAStruct.label_names_to_use = {'landscape_4003.avi', 'landscape_4004.avi', 'monkeyChasing_1111.avi', 'monkeyChasing_1112.avi', 'monkeyFighting_1121.avi', 'monkeyFighting_1123.avi', 'monkeyGoalDir_1101.avi', 'monkeyGoalDir_1102.avi', ...
    'monkeyGrooming_1141.avi', 'monkeyGrooming_1142.avi', 'monkeyIdle_1302.avi', 'monkeyIdle_1303.avi', 'monkeyMounting_1131.avi', 'monkeyMounting_1133.avi', 'objects_2101.avi', 'objects_2104.avi'};
  stimuliSetAStruct.coreTitle = 'Stim Set A';
  
  stimuliSetBStruct.labelFile = 'stimSetB';
  stimuliSetBStruct.label = 'stimuli';
  stimuliSetBStruct.label_names_to_use = {'landscape_4002.avi', 'landscape_4005.avi', 'monkeyChasing_1113.avi', 'monkeyChasing_1114.avi', 'monkeyFighting_1122.avi', 'monkeyFighting_1124.avi', 'monkeyGoalDir_1103.avi', 'monkeyGoalDir_1104.avi', ...
    'monkeyGrooming_1143.avi', 'monkeyGrooming_1144.avi', 'monkeyIdle_1301.avi', 'monkeyIdle_1305.avi', 'monkeyMounting_1132.avi', 'monkeyMounting_1134.avi', 'objects_2102.avi', 'objects_2103.avi'};
  stimuliSetBStruct.coreTitle = 'Stim Set B';
  
  stimSetStructs = [stimuliSetAStruct stimuliSetBStruct];
  [stimSetStructs.the_training_label_names] = deal([]);
  [stimSetStructs.the_test_label_names] = deal([]);
  [stimSetStructs.pseudoGenString] = deal('');
  [stimSetStructs.crossTempDecode] = deal(true);
  
  coreAnalysisStructs = [coreAnalysisStructs stimSetStructs];
end

% Add the specific binned file to use
[coreAnalysisStructs.binnedFile] = deal(binnedFileName);
[coreAnalysisStructs.paradigmName] = deal(paradigmName);

% Generate every possible cross between training and testing data, where
% for every stimulus category, you are training on stimCount - 1, then
% testing on the left out stim.
% 
% % 2 stimuli sets - find out which stimuli cluster
% rasterFileData = load(rasterFile, 'binned_labels');
% stimSet1Ind = find(cellfun(@(x) any(contains(x, 'objects_2103.avi')), rasterFileData.binned_labels.stimuli),1);
% stimSet2Ind = find(cellfun(@(x) any(contains(x, 'objects_2104.avi')), rasterFileData.binned_labels.stimuli),1);
% stimSetInd = [stimSet1Ind, stimSet2Ind];
% 
% for stimSet_i = 1:length(stimSetInd)
%   uniqueStimVec = unique([rasterFileData.binned_labels.stimuli{stimSetInd(stimSet_i)}])';
%   stimClassVec = unique(extractBefore(uniqueStimVec, '_'));
%   
%   stimByClass = cell(length(stimClassVec), 1);
%   for ii = 1:length(stimByClass)
%     stimByClass{ii} = uniqueStimVec(contains(uniqueStimVec, stimClassVec{ii}));
%   end
%   
%   % Each paradigm iterates across the below values
%   for ii = 1:2
%     % Use categoriesStruct as template, adding in a few new fields
%     stimuliStruct = categoriesStruct;
%     stimuliStruct.labelFile = sprintf('stimuliGen_stimSet%d_V%d', stimSet_i, ii);
%     stimuliStruct.pseudoGenString = 'stimuli';
%     stimuliStruct.genDS = true;
%     
%     stimuliStruct.label_names_to_use = uniqueStimVec;
% %     stimuliStruct.coreTitle = sprintf('Stimuli Gen, Stim Set %d, V%d', stimSet_i, ii);
%     logicalInd2Take = false(1,2);
%     logicalInd2Take(ii) = true;
%     
%     stimuliStruct.the_training_label_names = cellfun(@(x) x(logicalInd2Take), stimByClass, 'UniformOutput', false);
%     stimuliStruct.the_test_label_names = cellfun(@(x) x(~logicalInd2Take), stimByClass, 'UniformOutput', false);
%     
%     coreAnalysisStructs = [coreAnalysisStructs, stimuliStruct];
%     
%   end
% end

% All Generate analyses are generated with a combination of the features
% below.
decoders = {'max_correlation_coefficient_CL'};
decoderTags = {'MCC'};
preProcArray = {'zscore_normalize_FP'};
% unitSets = {'MUA', 'U&US', 'U'};
% unitSets = {'MUA', 'U'};
unitSets = {'U'};

% analysisTags are used to name the file, and put into the title of the
% figure.
% analyisTag = {'AllUnits', 'noHT', 'onlyHT', 'socIntSelAny', 'noSocIntSelAny', ...
%   'socIntEarly', 'socIntLate', 'socIntReward', ...
%   'socIntSel_noHT', 'socIntSel_onlyHT'...
%   'broadCatEarly', 'broadCatLate', 'broadCatReward', ...
%   'broadCatSliding', 'broadCatAny', ...
%   'broadCatAny_onlyHT', 'broadCatAny_noHT', 'nobroadCatAny_onlyHT'...
%   'stimEarlySel', 'stimLateSel', 'stimRewardSel'...
%   };

analyisTagDir = {'AllUnits', 'TaskMod', 'StimMod', 'StimMod_Mo', 'StimMod_Sam', 'Mo', 'Sam'};
analyisTagPlot = {'All', 'TaskMod', 'StimMod', 'StimMod M1', 'StimMod M2', 'M1', 'M2'};
  
% Matching set of variables to change, names must match what is defined in
% the raster file.
% analyisTagVars = {'', {{'subSel_headTurn_all_saccNone_selInd', 0}}, {{'subSel_headTurn_all_sacc_selInd', 0}} {{'subSel_headTurn_all_saccNone_selInd', 1}} {{'subSel_headTurn_all_sacc_selInd', 1}} {{'sVns_any_selInd', 1}}, {{'sVns_any_selInd', 0}}, ...
%   {{'sVns_stimEarly_selInd', 1}}, {{'sVns_stimLate_selInd', 1}}, {{'sVns_reward_selInd', 1}},...
%   {{'sVns_any_selInd', 1}, {'subSel_headTurn_all_selInd', 0}}, {{'sVns_any_selInd', 1}, {'subSel_headTurn_all_selInd', 1}}...
%   {{'bCat_stimEarly_selInd', 1}}, {{'bCat_stimLate_selInd', 1}}, {{'bCat_reward_selInd', 1}},...
%   {{'SW_broadCatTest_selInd', 1}}, {{'bCat_any_selInd', 1}}, ...
%   {{'bCat_any_selInd', 1}, {'subSel_headTurn_all_selInd', 1}}, {{'bCat_any_selInd', 1}, {'subSel_headTurn_all_selInd', 0}}, {{'bCat_any_selInd', 0}, {'subSel_headTurn_all_selInd', 1}}...
%   {{'baseV_stimEarly_selInd', 1}}, {{'baseV_stimLate_selInd', 1}}, {{'baseV_reward_selInd', 1}}};

analyisTagVars = {'', {{'taskModulated_selInd', 1}}, {{'taskModStim_selInd', 1}}, ...
  {{'taskModStim_selInd', 1} {'M1unit_selInd', 1}}, ...
  {{'taskModStim_selInd', 1} {'M1unit_selInd', 0}},...
  {{'M1unit_selInd', 1}}, ...
  {{'M1unit_selInd', 0}}, ...
};

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

analysisCount = 1;
for unit_set_i = 1:length(unitSets)
  % Set the base state of the binnedSiteInfoMatrix
  
  % Create a template - all switches should be set to true, except for
  % UnitType - which should be either MUA or U&US
  sites2KeepTemplate = true(size(binnedSiteInfoMatrix));
  unitTypeInd = strcmp(varFields, 'UnitType');
  
  switch unitSets{unit_set_i}
    case 'MUA'
      sites2KeepTemplate(unitTypeInd, :) = ismember(binnedSiteInfoMatrix(unitTypeInd,:), 1);
      unitTag = 'MUA';
    case 'U'
      sites2KeepTemplate(unitTypeInd, :) = ismember(binnedSiteInfoMatrix(unitTypeInd,:), 2);
      unitTag = 'U';
  end
  
  for core_i = 1:length(coreAnalysisStructs)
        
    analysisTemplateStruct = coreAnalysisStructs(core_i);
        
    analysisTemplateStruct.save_extra_preprocessing_info = 0;
    analysisTemplateStruct.num_features_to_use = 0;
    
    analysisTemplateStruct.coreTitle = sprintf('%s %s', paradigmNameTag, coreAnalysisStructs(core_i).coreTitle);
    
    % Run the initial k check
    labelPerSite = binnedLabels.(analysisTemplateStruct.label);
    [~, min_num_repeats_all_sites, ~, ~] = find_sites_with_k_label_repetitions(labelPerSite, 1, analysisTemplateStruct.label_names_to_use);
    
    for decoder_i = 1:length(decoders)
      for preProc_i = 1:length(preProcArray)
        for analysis_i = 1:length(analyisTagDir)
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
          
          % Extra step for GenDS - make sure all sites have non-0 number.
          sites2KeepArray = [sites2KeepArray; min_num_repeats_all_sites ~= 0];
          
          % Keep units which check all the boxes, store in the sites field.
          sites2Keep = all(sites2KeepArray,1);
          
          % Update sites to be used, and the repeats used.
          analysisStruct.sites = find(sites2Keep);
          new_k = min(min_num_repeats_all_sites(sites2Keep));
          [analysisStruct.num_cv_splits, analysisStruct.k_repeats_needed] = deal(new_k);
          
          % Generate the Title
          unitSubSecTag = sprintf(' (%s)', analyisTagPlot{analysis_i});
          
          % Create the filename.
          saveFileName = sprintf('%s_%s_%s_%s_%s', paradigmNameTag, analysisStruct.labelFile, analyisTagDir{analysis_i}, unitTag, decoderTags{decoder_i});
          
          % Modfiy things based on preprocessing, to create the directories
          % appropriately
          if ~isempty(analysisStruct.preProc)
            analysisSubType = [analysisStruct.classifier(1:3) '_' analysisStruct.preProc{1}(1:3)];
          else
            analysisSubType = analysisStruct.classifier(1:3);
          end
          
          if strcmp(analysisStruct.classifier, 'poisson_naive_bayes_CL')
            analysisStruct.load_data_as_spike_counts = true;
          else
            analysisStruct.load_data_as_spike_counts = false;
          end
          
          analysisStruct.plotOutDir = fullfile(pFolder, analysisSubType, strrep(analysisStruct.coreTitle, ' ', '_'));
          
          % If there is a pseudoGen label, add that to the name
          if ~isempty(analysisStruct.pseudoGenString)
            saveFileName = sprintf('%s_%s', saveFileName, analysisStruct.pseudoGenString);
            pseudoGenTag = sprintf('(pseudoGen - %s)', analysisStruct.pseudoGenString);
          else
            pseudoGenTag = '';
          end
          

          % Decide how to do cv_splits. either k-1 or 25/75
          if analysisChangeParams.expandLabelPerSplit
            analysisStruct.num_cv_splits = floor(analysisStruct.num_cv_splits/4);
            analysisStruct.num_times_to_repeat_each_label_per_cv_split = 4;
          else
            analysisStruct.num_times_to_repeat_each_label_per_cv_split = 1;
          end
          
          analysisStruct.cross_val_resample = params.cross_validator_num_resample;
          
          % Point to the correct directories
          analysisStruct.save_file_dir = fullfile(analysisStruct.plotOutDir, sprintf('results_%s', saveFileName));
          analysisStruct.shuff_file_dir = fullfile(analysisStruct.save_file_dir, 'shuff_dir');
          analysisStruct.shuffle_ds = false;
          
          % Can't be done in parafor loops
          analysisStruct.the_classifier = eval(analysisStruct.classifier);
          the_feature_preprocessors = cell(size(analysisStruct.preProc));
          for pp_i = 1:length(analysisStruct.preProc)
            the_feature_preprocessors{pp_i} = eval(analysisStruct.preProc{pp_i});
          end
          analysisStruct.the_feature_preprocessors = the_feature_preprocessors;
          analysisStruct.unitCount = length(analysisStruct.sites);

          % Create a shuffle distribution for all the comparisons
          % downstream
          shuffCount = params.null_shuffle_count;
          analysisStruct.coreAnalysis = saveFileName;           % Useful fields for later
          
          for shuff_i = 1:shuffCount
            analysisStruct.decoding_results_file = fullfile(analysisStruct.shuff_file_dir, sprintf('results_decoding_%d.mat', shuff_i));
            analysisStruct.shuffle_ds = true;

            repsTag = sprintf('(%d Reps %d %s) ', new_k, length(analysisStruct.sites), unitTag);
            analysisStruct.plotTitle = strcat(analysisStruct.coreTitle, unitSubSecTag, repsTag, pseudoGenTag);
            
            fileOutName = fullfile(analysisDir, sprintf('%s_results_shuffle_%d', saveFileName, shuff_i));
            save(fileOutName, 'analysisStruct');
            analysisStructArray{analysisCount} = fileOutName;
            analysisCount = analysisCount + 1;
          end
          analysisStruct.shuffle_ds = false;
          
          % Generate sub-samples.
          if subSampleRuns
            origOrder = analysisStruct.sites;
            for sub_i = 1:subsampleCount
              for unit_count_i = 1:length(unitsPerSubSample)
                
                if unitsPerSubSample(unit_count_i) <= length(origOrder)
                  % Find out subsample to take
                  reSampInd = sort(randperm(length(origOrder), unitsPerSubSample(unit_count_i)));
                  analysisStruct.sites = origOrder(reSampInd);
                  
                  % Set the titles
                  repsTag = sprintf('(%d Reps %d %s)', new_k, unitsPerSubSample(unit_count_i), unitTag);
                  analysisStruct.plotTitle = strcat(analysisStruct.coreTitle, unitSubSecTag, repsTag, pseudoGenTag);
                  
                  % save structure
                  saveFileNameSS = sprintf('%s_S%d_%dU', saveFileName, sub_i, unitsPerSubSample(unit_count_i));
                  analysisStruct.unitCount = unitsPerSubSample(unit_count_i);
                  analysisStruct.subInd = sub_i;
                  analysisStruct.decoding_results_file = fullfile(analysisStruct.save_file_dir, strcat('result_', saveFileNameSS, '.mat'));
                  save(fullfile(analysisDir, saveFileNameSS), 'analysisStruct');
                  analysisStructArray{analysisCount} = fullfile(analysisDir, saveFileNameSS);
                  analysisCount = analysisCount + 1;
                end
                
              end
            end
            
          else
            
            % Save structure
            analysisStruct.decoding_results_file = fullfile(analysisStruct.save_file_dir, strcat('result_', saveFileName, '.mat'));
            save(fullfile(analysisDir, saveFileName), 'analysisStruct');
            analysisStructArray{analysisCount} = fullfile(analysisDir, saveFileName);
            analysisCount = analysisCount + 1;
            
          end
          
        end
      end
    end
  end
end

end
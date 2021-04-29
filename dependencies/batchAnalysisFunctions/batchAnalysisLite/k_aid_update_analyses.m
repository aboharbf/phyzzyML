% K aid repeated
% a function which updates old analyses 

rasterFile = 'D:\DataAnalysis\batchAnalysis\NeuralDecodingTB\NaturalSocial\rasterData\rasterData_binned_150ms_bins_50ms_sampled.mat';
analysisDir = 'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML\buildAnalysisParamFileLib\NDT_analyses\NaturalSocial';
changeOldAnalyses = false;

% Load the raster data
rasterData = load(rasterFile);
binnedData = rasterData.binned_data;
binnedLabels = rasterData.binned_labels;
binnedSiteInfo = rasterData.binned_site_info;

% Load the previously generated analyses
prevAnalysisDir = dir(fullfile(analysisDir, '*.mat'));
prevAnalysisDir = fullfile({prevAnalysisDir.folder}, {prevAnalysisDir.name})';

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

% Cycle through the analysis files.
for file_i = 1:length(prevAnalysisDir)
 
  %load analysisStruct
  load(prevAnalysisDir{file_i}); 
  
  % Run the check
  k = analysisStruct.k_repeats_needed;
  label = binnedLabels.(analysisStruct.label);
  [available_sites_tmp, min_num_repeats_all_sites, num_repeats_matrix, ds.label_names_to_use] = find_sites_with_k_label_repetitions(label, 1, analysisStruct.label_names_to_use);
  
  % Site info selected
  siteInfoSelected = analysisStruct.site_info_selected;
  
  % Copy the template
  sites2Keep = sites2KeepTemplate;

  % Iterate and populate the sites2KeepTemplate, defining sites which meet
  % the criteria.
  for site_info_i = 1:length(siteInfoSelected)
    sites2Keep(site_info_i,:) = ismember(binnedSiteInfoMatrix(site_info_i,:), find(siteInfoSelected{site_info_i}));
  end
  
  % Keep units which check all the boxes
  sites2Keep = sum(sites2Keep,1) == size(sites2Keep,1);
  
  % Update the minimum number of repeats vector
  min_num_repeats_all_sites(~sites2Keep) = 0;
  
  % Get the highest number of repeats while utilizes as many units as
  % possible.
  new_k = min(min_num_repeats_all_sites(sites2Keep));
  sites_to_use = find(min_num_repeats_all_sites);
  
  % Update the analysis
  analysisStruct.sites = sites_to_use;
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
  newTitleTmp = originalTitle(1:removeInd(2)-1);
  newTitleTmp = [newTitleTmp sprintf('(%d Reps %d %s)', new_k, length(sites_to_use), unitTag)];
  analysisStruct.plotTitle = newTitleTmp;

  % save the structure back
  load(prevAnalysisDir{file_i}, 'analysisStruct');
  
end


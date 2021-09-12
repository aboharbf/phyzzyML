function report_mean_decoding_result_OneBin(analysisBatch, analysisStruct, params)

% Load the outputs from these analysisStructs
analysisOutputStructs = cellfun(@(x) load(x), analysisBatch);
analysisOutputStructs = [analysisOutputStructs.analysisStruct];
decodingResultsStructs = {analysisOutputStructs.decoding_results_file}';
decodingResultsStructs = cellfun(@(x) load(x), decodingResultsStructs);
decodingResultsStructs = [decodingResultsStructs.decoding_results];

% Extract clear variables from inputs.
if length(decodingResultsStructs) ~= 1
  decoding_results_tmp = decodingResultsStructs(1);
else
  decoding_results_tmp = decodingResultsStructs;
end

binParams = decoding_results_tmp.DS_PARAMETERS.binned_site_info.binning_parameters;
labels = decoding_results_tmp.DS_PARAMETERS.label_names_to_use; % Find Labels

if size(labels, 2) > 1
  labels = labels';
end

decoding_type = 'Zero';
switch decoding_type
  case 'Zero'
    decoding_data = [decodingResultsStructs.ZERO_ONE_LOSS_RESULTS];
  case 'Normalized'
    decoding_data = [decodingResultsStructs.NORMALIZED_RANK_RESULTS];
  case 'ROC_AUC'
    decoding_data = [decodingResultsStructs.ROC_AUC_RESULTS];
end

resultMean = mean([decoding_data.mean_decoding_results]) * 100;

% Report
fprintf('Mean Decoding for: %s = %s%% \n', analysisStruct.plotTitle, num2str(round(resultMean,2)))

end
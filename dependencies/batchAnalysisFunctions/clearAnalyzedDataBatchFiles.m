function clearAnalyzedDataBatchFiles(batchAnalysisParams)
%% Clears out all previously analzyed batch data files.

% Delete the spikePathData
delete([batchAnalysisParams.spikePathLoadParams.batchAnalysisOutput '*.mat'])

% Find and delete all the temporary files, and the preprocessed analyses.
analyzedStem = batchAnalysisParams.spikePathLoadParams.batchAnalysisOutputName;
prevFiles = dir(fullfile(batchAnalysisParams.analysisDirectory, '**', analyzedStem));
tmpFiles = dir([batchAnalysisParams.subEventPSTHParams.outputDir, filesep, '*.mat']);

prevFiles = [prevFiles; tmpFiles];
prevFiles = fullfile({prevFiles.folder}, analyzedStem)';

for ii = 1:length(prevFiles)
  delete(prevFiles{ii})
end
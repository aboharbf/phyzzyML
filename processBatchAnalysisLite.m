function spikeDataBankPath = processBatchAnalysisLite( varargin )
% Function uses batchAnalysisParamFile to coordinate the compilation of the
% desired 'spikeDataBank'. 
% Input
% - string of name of buildBatchAnalysisParamFile function
% - usePreprocessedSpikeData: a 1 or 0 indicating whether to regenerate the
% spikeDataBank (0), or use what is found (1). Default = 1.
% Output
% - paths to relevant .mat files (sliced up spikeDataBank + other
% variables). spikeDataBank files end with '_N', while the other variables
% end with 'Vars'. Vars file is 1st in cell array of file paths.

%% Parse Arguments
addpath(genpath('dependencies'));
rmpath(genpath('dependencies/mvgc_v1.0')); %note: Not sure if this is appropriate replacement for genpath_exclude. previous line caused issues in parallel runs.
addpath('buildAnalysisParamFileLib');

switch length(varargin)
  case 0
    disp('using buildBatchAnalysisLiteParamFileSocialVids');
    batchAnalysisParamFile  = 'buildBatchAnalysisLiteParamFileSocialVids';
    deleteOldAnalysis = false;
  case 1
    batchAnalysisParamFile  = varargin{1};
    deleteOldAnalysis = false;
  case 2
    batchAnalysisParamFile  = varargin{1};
    deleteOldAnalysis  = varargin{2};
  otherwise
    error('incorrect number of arguments')
end

%% Retrieve relevant parameters
addpath('buildAnalysisParamFileLib');
batchAnalysisParamFilename = feval(batchAnalysisParamFile);
batchAnalysisParams = load(batchAnalysisParamFilename);

%% If desired, clear out previously analyzed data
if deleteOldAnalysis
  clearAnalyzedDataBatchFiles(batchAnalysisParams)
end

%% Get the relevant data
spikePathFile = batchAnalysisParams.spikePathLoadParams.batchAnalysisOutput;
files = dir([spikePathFile '*.mat']);

if ~isempty(files)
  fprintf('spikePathBank found, loading... \n')
  load(spikePathFile);
else
  
  fprintf('spikePathBank Not found, generating now... \n')
  analysisDirectory = batchAnalysisParams.analysisDirectory;
  
  % Collected the two main files with data
  preprocessedList = dir([analysisDirectory filesep '**' filesep 'preprocessedData.mat']);
  analyzedList = dir([analysisDirectory filesep '**' filesep 'analyzedData.mat']);
  preprocessedPathList = fullfile({preprocessedList.folder}, {preprocessedList.name})';
  analyzedListPathList = fullfile({analyzedList.folder}, {analyzedList.name})';
  assert(length(preprocessedPathList) == length(analyzedListPathList), 'Lists arent same length, confirm every preprocessed file has an analyzed file')
  
  % Check for items which didn't make it to runAnalysis.
  analysisParamList = dir([analysisDirectory filesep '**' filesep 'AnalysisParams.mat']);
  analyzedListTmp = {analyzedList.folder}';
  paramList = {analysisParamList.folder}';
  
  % Show unanalyzed runs
  unanalyzedRuns = setdiff(paramList, analyzedListTmp);
  if ~isempty(unanalyzedRuns)
    fprintf('Found %d Unanalyzed runs.\n', length(unanalyzedRuns));
    disp(unanalyzedRuns);
%     error('All Runs must be analyzed');
  end
  fprintf('Found %d processed runs.\n', length(analyzedList));
  
  sessionList = cell(size(analyzedListPathList));
  spikePathBank = table();
  spikePathBank.analyzedDir = {preprocessedList.folder}';
  
  for ii = 1:length(analyzedListPathList)
    analyzedData = load(analyzedListPathList{ii}, 'analysisParamFilename','dateSubject', 'runNum');
    [~, sessionList{ii}] = deal(['S' analyzedData.dateSubject analyzedData.runNum]);
  end
  
  % Add session List as the row names.
  spikePathBank.Properties.RowNames = sessionList;
  
  if ~strcmp(batchAnalysisParams.monkey, 'Combo')
    monkeyInd = contains(spikePathBank.Properties.RowNames, batchAnalysisParams.monkey);
    spikePathBank = spikePathBank(monkeyInd, :);
    fprintf('spikePathBank reduces to only runs with %s \n', batchAnalysisParams.monkey);
  end
  
  % Save to folder
  save(spikePathFile, 'spikePathBank')
  
end

%Place the relevant paths into runBatchAnalysis.
runBatchAnalysisLite(spikePathBank, batchAnalysisParams)

end

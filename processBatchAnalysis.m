function spikeDataBankPath = processBatchAnalysis( varargin )
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
    disp('using buildBatchAnalysisParamFileSocialVids');
    batchAnalysisParamFile  = 'buildBatchAnalysisParamFileSocialVids';
    usePreprocessedSpikeData = 1;
  case 1
    batchAnalysisParamFile  = varargin{1};
    usePreprocessedSpikeData = 1;
  case 2
    batchAnalysisParamFile  = varargin{1};
    usePreprocessedSpikeData  = varargin{2};
  otherwise
    error('incorrect number of arguments')
end

%% Retrieve relevant parameters
addpath('buildAnalysisParamFileLib');
batchAnalysisParamFilename = feval(batchAnalysisParamFile);
load(batchAnalysisParamFilename);

%% Get the relevant data
spikeDataBaseFile = [outputDir '/' preprocessParams.spikeDataFileName];
files = dir([spikeDataBaseFile '*.mat']);

if ~isempty(files) && usePreprocessedSpikeData
  disp('spikeDataBank found, Returning paths')
  spikeDataBankPath = cell(length(files),1);
  for file_ind = 1:length(files)
    spikeDataBankPath{file_ind} = fullfile(files(file_ind).folder,files(file_ind).name);
  end
else
  fprintf('spikeDataBank Not found, generating now... \n')
  % Cycle through analyzedData.mat files, store and organize the relevant structures.
  preprocessedList = dir([analysisDirectory filesep '**' filesep 'preprocessedData.mat']);
  analyzedList = dir([analysisDirectory filesep '**' filesep 'analyzedData.mat']);
  assert(length(preprocessedList) == length(analyzedList), 'Lists arent same length, confirm every preprocessed file has an analyzed file')
  analysisParamList = dir([analysisDirectory filesep '**' filesep 'AnalysisParams.mat']);
  
  analyzedListTmp = {analyzedList.folder}';
  paramList = {analysisParamList.folder}';
  disp(setdiff(paramList, analyzedListTmp));
  
  fprintf('Found %d processed runs.\n', length(analyzedList));
  
  sessionList = cell(length(preprocessedList),1);
  spikeDataBank = struct();
  
  for ii = 1:length(preprocessedList)
    
    if mod(ii, 10) == 0
      fprintf('processing run %d... \n', ii)      
    end
    
    preprocessedData = load(fullfile(preprocessedList(ii).folder, preprocessedList(ii).name),'spikesByEvent','eventIDs','eventCategories','preAlign','postAlign');
    analyzedData = load(fullfile(analyzedList(ii).folder, analyzedList(ii).name), 'analysisParamFilename','dateSubject', 'runNum', 'groupLabelsByImage','psthByImage','psthErrByImage', 'stimStatsTable', 'subEventSigStruct', 'eyeDataStruct','spikesByEventBinned');
    analysisParams = load(fullfile(analyzedList(ii).folder, 'AnalysisParams.mat'), 'psthParams');
    
    [sessField, sessionList{ii}] = deal(['S' analyzedData.dateSubject analyzedData.runNum]);
    
    % fields from analysis params
    spikeDataBank.(sessField).start = -analysisParams.psthParams.psthPre;
    spikeDataBank.(sessField).stimDur = analysisParams.psthParams.psthImDur;
    spikeDataBank.(sessField).end = analysisParams.psthParams.psthImDur + analysisParams.psthParams.psthPost;
    
    % fields from preprocessedData
    spikeDataBank.(sessField).spikesByEvent = preprocessedData.spikesByEvent;
    spikeDataBank.(sessField).eventIDs = preprocessedData.eventIDs;
    spikeDataBank.(sessField).eventCategories = preprocessedData.eventCategories;
    spikeDataBank.(sessField).figDir = preprocessedList(ii).folder;

    % fields from analyzedData
    spikeDataBank.(sessField).dateSubject = analyzedData.dateSubject;
    spikeDataBank.(sessField).runNum = analyzedData.runNum;
    spikeDataBank.(sessField).spikesByEventBinned = analyzedData.spikesByEventBinned;
    spikeDataBank.(sessField).psthByImage = analyzedData.psthByImage;
    spikeDataBank.(sessField).psthErrByImage = analyzedData.psthErrByImage;
    spikeDataBank.(sessField).attendedObjData = analyzedData.eyeDataStruct.attendedObjData;
    spikeDataBank.(sessField).groupLabelsByImage = analyzedData.groupLabelsByImage;
    spikeDataBank.(sessField).tTest = analyzedData.stimStatsTable;
    spikeDataBank.(sessField).subEventSigStruct = analyzedData.subEventSigStruct;
    
    %spikeDataBank.(sessField).groupLabel = target; %Now in slidingANOVAparams.
  end
  spikeDataBankPath = saveSpikeDataBank(spikeDataBank, 12, 'save',outputDir);
  clearvars spikeDataBank
  nonSpikeSaveName = fullfile(outputDir, [preprocessParams.spikeDataFileName 'Vars']);
  spikeDataBankPath = [nonSpikeSaveName; spikeDataBankPath];
  save(nonSpikeSaveName)
end

%Place the relevant paths into runBatchAnalysis.
runBatchAnalysis(spikeDataBankPath)

end

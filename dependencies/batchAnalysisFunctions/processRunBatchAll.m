function processRunBatchAll()
% Function which analyzes all data from Mo and Sam.

% Things to change before next batch run:
% - Consider deleting all the figures in the folders to avoid confusion
% - create a saccade, non-saccade, and all version of events. Make sure
% eventTime indexing is robust to this change.
% - Introduce an 'analysisFinish' time into the analyzedData, so in case
% errors emerge, I can find them easier.
% - Consider stopping some testing which you know you wont use if its
% wasting time.

tic
[errorStack, errorStackMsg, files2Check] = processRunBatch('buildAnalysisParamFileSocialVids', {'naturalSocial'}, 'Mo');
MoFin = datetime;

[errorStack2, errorStackMsg2, files2Check2] = processRunBatch('buildAnalysisParamFileSocialVidsSam', {'naturalSocial'}, 'Sam');
SamFin = datetime;
toc

[errorStack3, errorStackMsg3, files2Check3] = processRunBatch('buildAnalysisParamFileSocialVids', {'headTurnCon'}, 'Mo');
MoFin2 = datetime;

[errorStack4, errorStackMsg4, files2Check4] = processRunBatch('buildAnalysisParamFileSocialVidsSam', {'headTurnCon'}, 'Sam');
SamFin2 = datetime;
toc


errorStackAll = [errorStack, errorStack2];
errorStackMsgAll = [errorStackMsg; errorStackMsg2];
files2CheckAll = [files2Check; files2Check2];

if strcmp(errorStack, 'None') && strcmp(errorStack2, 'None')
  processBatchAnalysisLite('buildBatchAnalysisLiteParamFileSocialVids', true)
else
  error('Errors occured');
end

end
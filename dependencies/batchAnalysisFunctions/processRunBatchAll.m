function processRunBatchAll()
% Function which analyzes all data from Mo and Sam.

tic
[errorStack, errorStackMsg, files2Check] = processRunBatch('buildAnalysisParamFileSocialVids', {'naturalSocial'}, 'Mo');
MoFin = datetime;

[errorStack2, errorStackMsg2, files2Check2] = processRunBatch('buildAnalysisParamFileSocialVidsSam', {'naturalSocial'}, 'Sam');
SamFin = datetime;
toc

[errorStack4, errorStackMsg4, files2Check4] = processRunBatch('buildAnalysisParamFileSocialVidsSam', {'headTurnCon'}, 'Sam');
SamFin2 = datetime;
toc

[errorStack3, errorStackMsg3, files2Check3] = processRunBatch('buildAnalysisParamFileSocialVids', {'headTurnCon'}, 'Mo');
MoFin2 = datetime;

errorStackAll = [errorStack, errorStack2];
errorStackMsgAll = [errorStackMsg; errorStackMsg2];
files2CheckAll = [files2Check; files2Check2];

if strcmp(errorStack, 'None') && strcmp(errorStack2, 'None')
  processBatchAnalysisLite('buildBatchAnalysisLiteParamFileSocialVids', true)
else
  error('Errors occured');
end

end
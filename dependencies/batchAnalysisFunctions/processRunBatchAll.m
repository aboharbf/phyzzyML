function processRunBatchAll()
% Function which analyzes all data from Mo and Sam.

tic
[errorStack, errorStackMsg, files2Check] = processRunBatch('buildAnalysisParamFileSocialVids', {'naturalSocial', 'headTurnCon'}, 'Mo');
MoFin = datetime;

[errorStack2, errorStackMsg2, files2Check2] = processRunBatch('buildAnalysisParamFileSocialVidsSam', {'naturalSocial', 'headTurnCon'}, 'Sam');
SamFin = datetime;
toc

if strcmp(errorStack, 'None') && strcmp(errorStack2, 'None')
  processBatchAnalysisLite('buildBatchAnalysisLiteParamFileSocialVids', true)
else
  error('Errors occured');
end

end
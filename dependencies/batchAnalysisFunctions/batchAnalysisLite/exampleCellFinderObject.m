function exampleCellFinderObject(selTable, params)
% A function which looks through the selTable, finds units defined as
% selective for a particular instance, and creates a list in sorted order.

tableVars = selTable.Properties.VariableNames';
pValVars = tableVars(contains(tableVars, 'epochSel') & contains(tableVars, 'pVal'));
varName = extractBetween(pValVars, 'epochSel_', '_pVal');
selIndVars = strcat(extractBefore(pValVars, 'pVal'), 'selInd');
prefStimVars = strcat(extractBefore(pValVars, 'pVal'), 'prefStim');
siteTag = strcat(string(selTable.dateSubj), string(selTable.runNum), string(selTable.channel), string(selTable.unitType));

for ii = 1:length(pValVars)
  % Identify the selective units (selInd)
  selIndVar = selTable.(selIndVars{ii});
  pValVar = selTable.(pValVars{ii});
  prefStimVar = selTable.(prefStimVars{ii});
  
  % Grab the corresponding differences (cohensD) by unit
  
  sigDiffs = pValVar(selIndVar);
  sigDiffSites = siteTag(selIndVar);
  sigDiffObj = prefStimVar(selIndVar);
  
  % Sort them based on difference size
  [sigDiffsSorted, B] = sort(sigDiffs, 'ascend');
  sigDiffSitesSorted = sigDiffSites(B);
  sigDiffObjSorted = sigDiffObj(B);
  
  % Print
  if isempty(sigDiffs)
    fprintf('==== %s Selective for %s - is empty ==== \n', params.unitTag, varName{ii})
  else
    fprintf('==== %s Selective for %s ==== \n', params.unitTag, varName{ii})
    for jj = 1:min(length(sigDiffSitesSorted), 10)
      fprintf('%s - Site %s - Pref Obj %s \n', num2str(sigDiffsSorted(jj), 3), sigDiffSitesSorted{jj}, sigDiffObj{jj})
    end
  end
  
end
end
function exampleCellFinderObject(selTable)
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
  for unitType_i = 1:2
    
    % ID unit index
    if unitType_i == 1
      unitInd = contains(selTable.unitType, 'MUA');
      tag = 'MUA';
    else
      unitInd = ~contains(selTable.unitType, 'MUA');
      tag = 'U&US';
    end
    
    sigDiffs = pValVar(selIndVar & unitInd);
    sigDiffSites = siteTag(selIndVar & unitInd);
    sigDiffObj = prefStimVar(selIndVar & unitInd);
    
    % Sort them based on difference size
    [sigDiffsSorted, B] = sort(sigDiffs, 'ascend');
    sigDiffSitesSorted = sigDiffSites(B);
    sigDiffObjSorted = sigDiffObj(B);
    
    % Print
    if isempty(sigDiffs)
      fprintf('==== %s Selective for %s - is empty ==== \n', tag, varName{ii})
    else
      fprintf('==== %s Selective for %s ==== \n', tag, varName{ii})
      for jj = 1:min(length(sigDiffSitesSorted), 10)
        fprintf('%s - Site %s - Pref Obj %s \n', num2str(sigDiffsSorted(jj), 3), sigDiffSitesSorted{jj}, sigDiffObj{jj})
      end
    end
  end
  
end
end
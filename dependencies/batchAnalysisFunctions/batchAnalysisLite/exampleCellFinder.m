function exampleCellFinder(selTable)
% A function which looks through the selTable, finds units defined as
% selective for a particular instance, and creates a list in sorted order.

tableVars = selTable.Properties.VariableNames';
cohensDVars = tableVars(contains(tableVars, 'subSel') & contains(tableVars, 'cohensD'));
varName = extractBetween(cohensDVars, 'subSel_', '_cohensD');
selIndVars = strcat(extractBefore(cohensDVars, 'cohensD'), 'selInd');
siteTag = strcat(string(selTable.dateSubj), string(selTable.runNum), string(selTable.channel), string(selTable.unitType));

for ii = 1:length(cohensDVars)
  % Identify the selective units (selInd)
  selIndVar = selTable.(selIndVars{ii});
  cohensDVar = selTable.(cohensDVars{ii});
  
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
    
    sigDiffs = cohensDVar(selIndVar & unitInd);
    sigDiffSites = siteTag(selIndVar & unitInd);
    
    % Sort them based on difference size
    [sigDiffsSorted, B] = sort(sigDiffs, 'descend');
    sigDiffSitesSorted = sigDiffSites(B);
    
    % Print
    if isempty(sigDiffs)
    fprintf('==== %s Selective for %s - is empty ==== \n', tag, varName{ii})
    else
    fprintf('==== %s Selective for %s ==== \n', tag, varName{ii})
    for jj = 1:min(length(sigDiffSitesSorted), 10)
      fprintf('%s - Site %s \n', num2str(sigDiffsSorted(jj), 3), sigDiffSitesSorted{jj})
    end
    end
  end
  
end
end
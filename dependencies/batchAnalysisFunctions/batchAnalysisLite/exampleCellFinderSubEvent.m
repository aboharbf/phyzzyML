function exampleCellFinderSubEvent(selTable, params)
% A function which looks through the selTable, finds units defined as
% selective for a particular instance, and creates a list in sorted order.

tableVars = selTable.Properties.VariableNames';
cohensDVars = tableVars(contains(tableVars, 'subSel') & contains(tableVars, 'cohensD'));
varName = extractBetween(cohensDVars, 'subSel_', '_cohensD');
selIndVars = strcat(extractBefore(cohensDVars, 'cohensD'), 'selInd');
keepIdx = contains(selIndVars, tableVars);
selIndVars = selIndVars(keepIdx);
siteTag = strcat(string(selTable.dateSubj), string(selTable.runNum), '_', string(selTable.channel), string(selTable.unitType));

for ii = 1:length(selIndVars)
  % Identify the selective units (selInd)
  selIndVar = selTable.(selIndVars{ii});
  cohensDVar = selTable.(cohensDVars{ii});
  
  % Grab the corresponding differences (cohensD)
  sigDiffs = cohensDVar(selIndVar);
  sigDiffSites = siteTag(selIndVar);
  
  % Sort them based on difference size
  [sigDiffsSorted, B] = sort(sigDiffs, 'descend');
  sigDiffSitesSorted = sigDiffSites(B);
  
  % Print
  if isempty(sigDiffs)
    fprintf('==== %s Selective for %s - is empty ==== \n', params.selParam.unitTag, varName{ii})
  else
    
    fprintf('==== %s Selective for %s ==== \n', params.selParam.unitTag, varName{ii})
    for jj = 1:min(length(sigDiffSitesSorted), 10)
      fprintf('%s - Site %s \n', num2str(sigDiffsSorted(jj), 3), sigDiffSitesSorted{jj})
    end
    
    % Depending on what you're looking to show, use a different filter
    if contains(selIndVars{ii}, 'reward')
      % open the PSTHes.
        openUnitFigs(sigDiffSitesSorted(1:2), 'saccadeRaster_*', '', params.analysisDirectory)
      % openUnitPSTHes(sortedUnitNames(end:-1:end-units2Report+1), batchAnalysisParams.analysisDirectory)
    end
    
  end

end
end
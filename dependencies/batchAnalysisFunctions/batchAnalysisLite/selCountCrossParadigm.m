function selCountCrossParadigm(spikePathBank, selTablePerRun, batchAnalysisParams);
% function which determines which sites/kinds of selectivity exists across
% paradigms. 

paradigmList = unique(spikePathBank.paradigmName);
runNames = extractAfter(spikePathBank.Properties.RowNames, 'S');

% Site identification - for the sake of saying which MUA are the same site
% across paradigms
for run_i = 1:length(selTablePerRun)
  if ~isempty(selTablePerRun{run_i})
    % This siteInfoVec will vary depending on question. to find actual runs
    % at the same site, leave off channel. To do more granular looks at
    % units which show some selectivity across paradigms, include it.
    siteInfoChVec = [selTablePerRun{run_i}.dateSubj selTablePerRun{run_i}.gridHole selTablePerRun{run_i}.channel selTablePerRun{run_i}.unitType selTablePerRun{run_i}.recDepth];
    siteInfoVec = [selTablePerRun{run_i}.dateSubj selTablePerRun{run_i}.gridHole selTablePerRun{run_i}.recDepth];
    
    allSitesCh = join(siteInfoChVec, '_');
    allSites = join(siteInfoVec, '_');
    
    selTablePerRun{run_i}.siteChID = allSitesCh;
    selTablePerRun{run_i}.siteID = allSites;
    
  end
end

% for each paradigm, make a list of unique sites
siteListPerParadigm = cell(size(paradigmList));
for par_i = 2:length(paradigmList)
  pInd = strcmp(spikePathBank.paradigmName, paradigmList{par_i});
  selTableParadigmPerRun = selTablePerRun(pInd);
  
  % Combine across tables
  selTableParadigm = vertcat(selTableParadigmPerRun{:});
  siteListPerParadigm{par_i} = unique(selTableParadigm.siteID);
  
end

% What NS sites have HTC, HTC, or both?
NS_HTC = intersect(siteListPerParadigm{2}, siteListPerParadigm{3});
NS_HTI = intersect(siteListPerParadigm{2}, siteListPerParadigm{4});
NS_HTC_HTI = intersect(NS_HTC, NS_HTI);


% Find the Natural Social sites with selectivity of interest
pInd = strcmp(spikePathBank.paradigmName, 'NaturalSocial');
selTableNS = selTablePerRun(pInd);
% Combine across tables
selTableNS = vertcat(selTableNS{:});
selTableNS = expandSelTableComboEvents(selTableNS, batchAnalysisParams.selParam);

% Find the Natural Social sites with selectivity of interest
pInd = strcmp(spikePathBank.paradigmName, 'headTurnCon');
selTableHTC = selTablePerRun(pInd);
% Combine across tables
selTableHTC = vertcat(selTableHTC{:});
selTableHTC = expandSelTableComboEvents(selTableHTC, batchAnalysisParams.selParam);

NS_Sites = selTableNS.siteChID;
HTC_Sites = selTableHTC.siteChID;
[sitesInBoth, NS_ind, HTC_ind] = intersect(NS_Sites, HTC_Sites);

crossParCheck = {'subSel_headTurn_all', 'socIntSel_any', 'subSel_reward', 'subSel_rewardAbsent', 'saccSel', 'fixationSel', 'socIntSel_reward'};

for cross_i = 1:length(crossParCheck)
  
  NS_sel = (selTableNS.(crossParCheck{cross_i}) ~= 0);
  HTC_sel = (selTableHTC.(crossParCheck{cross_i}) ~= 0);
  
  % Which sites exist in both
  NS_sel_both = NS_sel(NS_ind);
  HTC_sel_both = HTC_sel(HTC_ind);
  
  selBothInd = NS_sel_both & HTC_sel_both;
  
  % Which sites are selective in each
  NS_sel_count = sum(NS_sel);
  HTC_sel_count =  sum(HTC_sel);
  
  NS_sel_rate = sum(NS_sel)/length(NS_sel);
  HTC_sel_rate =  sum(HTC_sel)/length(HTC_sel);
  chanceRate = NS_sel_rate * HTC_sel_rate;
   
  % jeez
  sitesSelBoth = sitesInBoth(selBothInd);
  chanceCount = length(sitesInBoth) * chanceRate;
  
  fprintf('Comparing the field %s, %d sites in both, chance count %d \n', crossParCheck{cross_i}, length(sitesSelBoth), round(chanceCount))
  
end


end


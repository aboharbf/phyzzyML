function stimulusGrid(spikePathBank, params)

tmp = load(params.eventDataPath);
eventListPlot = strrep(tmp.eventData.Properties.VariableNames, '_', ' ');

allStimuliVec = unique(vertcat(spikePathBank.stimuli{strcmp(spikePathBank.paradigmName, 'NaturalSocial')}));
stimPerRun = spikePathBank.stimuli(strcmp(spikePathBank.paradigmName, 'NaturalSocial'));
stimSets = {stimPerRun{2}; stimPerRun{3}};
dataMat = zeros(2,2);

for set_i = 1:length(stimSets)
  eventData = tmp.eventData(stimSets{set_i},:);
  % eventData = tmp.eventData;
  allStimuliVec = eventData.Properties.RowNames;
  tokens2Find = {'.avi', '_\d{3}', 'monkey', 'human'};
  tokens2Rep = {'', '', 'm', 'h'};
  allStimuliVecNames = regexprep(allStimuliVec, tokens2Find, tokens2Rep);
  allStimuliVecNames = strrep(allStimuliVecNames, '_', ' ');  % For the Animation tags.
  
  % Counts of events per stimulus, specifically headTurning
  monkeyId = contains(allStimuliVec, 'monkey');
  allStimuliVec = allStimuliVec(monkeyId);
  goalInd = find(contains(allStimuliVec, 'Goal'));
  idleInd = find(contains(allStimuliVec, 'Idle'));
  socInd = find(~contains(allStimuliVec, 'Goal') & ~contains(allStimuliVec, 'Idle'));
  
  allStimuliVec = allStimuliVec([socInd; idleInd; goalInd]);
  
  eventData = eventData(allStimuliVec, :);
  eventList = eventData.Properties.VariableNames;
  headTurnInd = contains(eventList, 'headTurn');
  eventList = eventList(headTurnInd);
  
  eventCountMat = cellfun(@(x) size(x, 1), table2cell(eventData(:,headTurnInd)));
  eventPerStim = sum(eventCountMat, 2);
  eventsPerSocial = eventPerStim(1:length(socInd));
  eventsPerNonSocial = eventPerStim(length(socInd)+1:end);
  
  fprintf('Events in Social %d, Events in Non-Social %d \n', sum(eventsPerSocial), sum(eventsPerNonSocial));
  
  % Plot
  figTitle = 'Events in Stimuli';
  figH = figure('units', 'normalized', 'position', [0.3542    0.1398    0.4573    0.7657]);
  imagesc(eventCountMat)
  figH.Children(1).XTick = 1:size(eventCountMat, 2);
  figH.Children(1).XTickLabel = strrep(eventList, '_', ' ');
  figH.Children(1).YTick = 1:size(eventCountMat, 1);
  figH.Children(1).YTickLabel = strrep(allStimuliVec, '_', ' ');
  title(figTitle)
  colorbar()
  
  dataMat(:,set_i) = [sum(eventsPerSocial); sum(eventsPerNonSocial)];
end

figTitle = 'Events in Stimuli';
figH = figure('Name', figTitle, 'NumberTitle','off','units','normalized', 'outerposition', [.3 .3 .4 .6]);
colNamesPlot = {'Stimulus Set 1', 'Stimulus Set 2'};
X = categorical(colNamesPlot);
X = reordercats(X,colNamesPlot);
bar(X, dataMat');
figH.Children(1).FontSize = 16;
ylabel('Head Turning Count');
legend('Head Turning Events in Social Videos', 'Head Turning Events in Agent Controls');

end
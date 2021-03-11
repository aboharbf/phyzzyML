function spikePathBank = stimulusStatistics(spikePathBank, params)
% Stimuli Presentation count and 'Novelty' Related Information.
%Code below creates a single large vector of stimuli used, and uses this to
%create individual vectors containing which viewing of the stimulus this
%represent (i.e. 'this run represents the 10th viewing of X.avi'). It also
%appends a dateTime vector to each structure related to how long since the
%last recording day.

runList = spikePathBank.Properties.RowNames;

tok2Find = {'S20', 'Mo00'};
tok2Rep = {'', 'R'};
runListLabels = regexprep(runList, tok2Find, tok2Rep);

% Extract the eventIDs field, generate a cell array of unique stimuli
[allStimCatVecPerRun, allCategoryVecPerRun, dateSubjectVec] = spikePathLoad(spikePathBank, {'eventIDs', 'categoryList', 'dateSubject'}, params.spikePathLoadParams);
allStimuliVec = unique(vertcat(allStimCatVecPerRun{:}));
allCategoryVec = unique(horzcat(allCategoryVecPerRun{:}))';

% Attach the paradigm to the table for later reference
paradigmIndex = nan(size(allStimCatVecPerRun));
paradigmName = cell(size(allStimCatVecPerRun));

for run_i = 1:length(allStimCatVecPerRun)
  
  runStimList = allStimCatVecPerRun{run_i};
  if any(contains(runStimList, 'headTurnCon'))
    paradigmIndex(run_i) = 2;
    paradigmName{run_i} = 'headTurnCon';
  elseif any(contains(runStimList, 'headTurnIso'))
    paradigmIndex(run_i) = 3;
    paradigmName{run_i} = 'headTurnIso';
  elseif any(contains(runStimList, 'monkey'))
    paradigmIndex(run_i) = 1;
    paradigmName{run_i} = 'NaturalSocial';
  elseif any(contains(runStimList, 'Otis'))
    paradigmIndex(run_i) = 4;
    paradigmName{run_i} = 'FamiliarFace';
  end

end

tokens2Find = {'.avi', '_\d{3}', 'monkey', 'human'};
tokens2Rep = {'', '', 'm', 'h'};
allStimuliVecNames = regexprep(allStimuliVec, tokens2Find, tokens2Rep);
allStimuliVecNames = strrep(allStimuliVecNames, '_', ' ');  % For the Animation tags.

% Produce matrix (N stim * M runs) which gives 0 for non present stim, count of stim presentation otherwise.
stimLogicalArray = zeros(length(allStimuliVec),length(runList));
for run_ind = 1:length(runList)
  stimLogicalArray(:,run_ind) = ismember(allStimuliVec, allStimCatVecPerRun{run_ind});
end

% Turn that matrix into a matrix of nth presentation.
csStimLogicalArray = cumsum(stimLogicalArray,2);
csStimLogicalArray(~stimLogicalArray) = 0;

% When was a stimulus first seen? Index of runList where first presentation took place.
[firstStimPresInd, ~] = find(csStimLogicalArray' == 1);

% Sort the previous count array by the first presentations
[~, newInd] = sort(firstStimPresInd);
csStimLogicalArraySorted = csStimLogicalArray(newInd, :);
allStimuliVecNames = allStimuliVecNames(newInd, :);

% If this is run, it should be saved
params.figStruct.saveFig = 1;
params.figStruct.closeFig = 0;
sliceStart = [1 40 80];
sliceEnd = [39 79 107];

% Generate image figure for this using XYGrid
for ii = 1:length(sliceStart)
  figTitle = sprintf('StimRunGrid_%d', ii);
  h = figure('NumberTitle', 'off', 'Name', figTitle, 'units', 'normalized', 'outerposition', [0 0 1 1]);
  XYGrid(h, csStimLogicalArraySorted(sliceStart(ii):sliceEnd(ii), :), allStimuliVecNames(sliceStart(ii):sliceEnd(ii)), runListLabels, params.xyStimParams);
  saveFigure(params.outDir, figTitle, [], params.figStruct, [])
end

% Generate image figure for eventData using XYGrid
tmp = load(params.eventDataPath);
eventListPlot = strrep(tmp.eventData.Properties.VariableNames, '_', ' ');
try
  eventData = tmp.eventData(allStimuliVec,:);
  
  for ii = 1:length(sliceStart)
    figTitle = sprintf('StimEventGrid_%d', ii);
    h = figure('NumberTitle', 'off', 'Name', figTitle, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    XYGrid(h, eventData(sliceStart(ii):sliceEnd(ii), :), allStimuliVecNames(sliceStart(ii):sliceEnd(ii)), eventListPlot, params.xyEventparams)
    saveFigure(params.outDir, figTitle, [], params.figStruct, [])
  end
  
  % Simple subEvent visualization, which happen when.
  eventMat = generateEventImage(tmp.eventData, 2800); % Hard coded stim length.
  eventMat = eventMat(:, 1:2800, :);  % Remove excess due to events persisting longer than stim was shown.
  for event_i = 1:length(eventListPlot)
    for ii = 1:length(sliceStart)
      figTitle = sprintf('Time of %s event occurances - %d', eventListPlot{event_i}, ii);
      h = figure('NumberTitle', 'off', 'Name', figTitle, 'units', 'normalized', 'outerposition', [0 0 1 1]);
      params.xyEventparams.plotTitle = figTitle;
      XYGrid(h, eventMat(sliceStart(ii):sliceEnd(ii), :, event_i), allStimuliVecNames(sliceStart(ii):sliceEnd(ii)), [], params.xyEventparams)
      saveFigure(params.outDir, figTitle, [], params.figStruct, [])
    end
  end
catch
  disp('Not all events represented in eventData table');
end

% Append a dateTime to each field with the time in days since the last
% recording. Add the relevant slice of the larger csStimLogicalArray.
allDateTimeVec = NaT(size(runList));
for run_ind = 1:length(runList)
  allDateTimeVec(run_ind) = datetime(extractBetween(dateSubjectVec{run_ind},1,8),'InputFormat','yyyyMMdd');     %Generate 'daysSinceLastRec' for each field.
end

% find unique recording dates, and the distance between them. Add these
% to the spikeDataBank later.
uniqueDateTimeVec = unique(allDateTimeVec);
daysSinceLastRec = [1000; days(diff(uniqueDateTimeVec))];

% use the dateTime and stimulus presentation matrix to find out how long
% in days takes place before a particular showing of a stimulus.
daysSinceLastPres = zeros(size(stimLogicalArray));
for stim_ind = 1:size(stimLogicalArray,1)
  presentationInd = logical(stimLogicalArray(stim_ind,:)); %When was the stim shown
  daysSinceLastPres(stim_ind,presentationInd) = [1000; days(diff(allDateTimeVec(presentationInd)))]; %Duration between those dates in days
end


[daysSinceLastPresVec, stimPresCountVec, stimPresArrayVec] = deal(cell(size(runList)));
daysSinceLastRecVec = nan(size(runList));
for run_ind = 1:size(stimLogicalArray,2)
  [~, big2SmallInd] = ismember(allStimCatVecPerRun{run_ind}, allStimuliVec);
  daysSinceLastPresVec{run_ind} = daysSinceLastPres(big2SmallInd,run_ind);
  stimPresCountVec{run_ind} = csStimLogicalArray(big2SmallInd,run_ind);
  stimPresArrayVec{run_ind} = csStimLogicalArray(:,run_ind);
  daysSinceLastRecVec(run_ind) = daysSinceLastRec(allDateTimeVec(run_ind) == uniqueDateTimeVec);
end

% Make Catagories column vectors
allCategoryVecPerRun = cellfun(@(x) x', allCategoryVecPerRun, 'UniformOutput', false);

% Fill out spikePathBank with generated values.
spikePathBank.stimuli = allStimCatVecPerRun;
spikePathBank.categories = allCategoryVecPerRun;
spikePathBank.paradigmInd = paradigmIndex;
spikePathBank.paradigmName = paradigmName;
spikePathBank.dateTime = allDateTimeVec;
spikePathBank.daysSinceLastPres = daysSinceLastPresVec;
spikePathBank.stimPresCount = stimPresCountVec;
spikePathBank.stimPresArray = stimPresArrayVec;
spikePathBank.daysSinceLastRec = daysSinceLastRecVec;

end

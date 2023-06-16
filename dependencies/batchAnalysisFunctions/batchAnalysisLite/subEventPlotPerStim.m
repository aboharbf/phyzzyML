function subEventPlotPerStim()
% Function which turns the eventData table into an image, where each row is
% a stimulus, each column is a time bin, and each row contains 4 different
% lines demarking the 4 different events in the eventData.

eventDataFile = 'C:\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories\eventData.mat';
outputDir = 'C:\OneDrive\Lab\Documents\Figures April 2020';

figStruct = struct();
figStruct.saveFig = 1;                
figStruct.closeFig = 0;               
figStruct.exportFig = 1;             
figStruct.saveFigData = 0;            

load(eventDataFile); % Puts eventData into the workspace.

% Stim in table - Make the table just stimuli in the 'NaturalSocial'
% paradigm.
stimInTable = eventData.Properties.RowNames;
stimInTable = extractBefore(stimInTable, '.avi');
stimInTable = stimInTable(1:68);
stimInTable = stimInTable(~contains(stimInTable, 'Following') & ~contains(stimInTable, 'Foraging') & ~contains(stimInTable, 'Holding') & ~contains(stimInTable, 'human') & ~contains(stimInTable, 'landscape') & ~contains(stimInTable, '5'));
stimInTable = stimInTable(~contains(stimInTable, 'human') & ~contains(stimInTable, 'landscape') & ~contains(stimInTable, '5'));
stimInTable = stimInTable(~contains(stimInTable, 'monkeyIdle_1304'));
stimInTable = [stimInTable; 'monkeyIdle_1305'];

% Extract the data from the table, use generateEventImage to make into
% matrix.
eventDataStim = eventData(strcat(stimInTable,'.avi'), :);
eventList = eventDataStim.Properties.VariableNames(1:4);
eventDataStim = eventDataStim(:, eventList);
eventList = strrep(eventList, '_', ' ');
eventMat = generateEventImage(eventDataStim, 2800);

% Make plot Labels.
stimInTablePlot = strcat(extractBetween(stimInTable, 'monkey', '_'), cellfun(@(x) x(end), stimInTable));
stimInTable = strrep(stimInTable, '_', '');

% Resort data and labels
goalDirInd = contains(stimInTablePlot, 'GoalDir');
idleInd = contains(stimInTablePlot, 'Idle');
socInd = ~goalDirInd & ~idleInd;

resortInd = [find(socInd); find(goalDirInd); find(idleInd)];
stimInTablePlot = stimInTablePlot(resortInd);
eventMat = eventMat(resortInd, :, :);

figTitle = sprintf('Sub event occurances in Stimuli');
eventMat2Figure(eventMat, stimInTablePlot, eventList, figTitle)

saveFigure(outputDir, figTitle, [], figStruct, []);

end
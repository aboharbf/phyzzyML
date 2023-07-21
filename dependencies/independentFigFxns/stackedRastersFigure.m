function stackedRastersFigure()
% a function which stacks rasters, and create a line PSTH at the top. for
% every raster submitted, this raster adds a raster to the bottom of the
% stack, and another line to the line plot at the top.

% Targets to load
epochNames = ["Fixation", "Stimulus Early", "Stimulus Late", "Reward"];
epochLabelXVals = [-850, 1, 1000, 2850];
targetDirectories = ['D:\DataAnalysis\20201201Mo\Basic\004', "D:\DataAnalysis\20201206Mo\Basic\001"...
  "D:\DataAnalysis\20201204Mo\Basic\003", "D:\DataAnalysis\20210630Sam\Basic\003"];
targetUnitLabels = [23 2; 63 2; 14 2; 13 2];

figOutDir = "D:\batchAnalysis\selCount";

% Generate Unit Names
targetDirParts = cellfun(@(x) strsplit(x, filesep), targetDirectories, 'UniformOutput', false);
subjDate = cellfun(@(x) x(end-2), targetDirParts);
runNum = cellfun(@(x) x(end), targetDirParts);
chanString = strcat('Ch', string(targetUnitLabels(:,1)));
unitString = strcat('U', string(targetUnitLabels(:,2)-1));
unitName = strcat(subjDate', runNum', chanString, unitString);

% For every entry above...
[spikeData, psthParamsPerRaster] = deal(cell(length(targetDirectories),1));

for raster_i = 1:length(targetDirectories)
  
  % Grab the file with the spikeData.
  if raster_i == 1
    tmp = load(fullfile(targetDirectories{raster_i}, 'preprocessedData.mat'), 'spikesByEventFixAlign'); %loads 'spikesByEvent'
    spikesByEvent = tmp.spikesByEventFixAlign;
  else
    load(fullfile(targetDirectories{raster_i}, 'preprocessedData.mat'), 'spikesByEvent'); %loads 'spikesByEvent'
  end
  load(fullfile(targetDirectories{raster_i}, 'AnalysisParams.mat'), 'psthParams'); %loads 'spikesByEvent'
  load(fullfile(targetDirectories{raster_i}, 'AnalysisParams.mat'), 'spikeAlignParams'); % loads 'lfpAlignParams
%   load(fullfile(targetDirectories{raster_i}, 'AnalysisParams.mat'), 'lfpAlignParams'); % loads 'lfpAlignParams

  psthParamsPerRaster{raster_i} = psthParams;

  % Stack all the spikes for the selected target unit together
  chan_i = targetUnitLabels(raster_i, 1);
  unit_i = targetUnitLabels(raster_i, 2);
  spikeDataTmp = [];
  for event_i = 1:length(spikesByEvent)
    spikeDataTmp =[spikeDataTmp; spikesByEvent{event_i}{chan_i}{unit_i}];
  end
  
  % Add to the stack - pretending each raster is actually activity for a
  % different event, on chan 1 unit 1 - for the sake of later functions.
  spikeData{raster_i}{1}{1} = spikeDataTmp;
  
end

% Shifts the spikes in the fix aligned first example cell back.
fixData = spikeData{1}{1}{1};
for trial_i = 1:length(fixData)
  fixData(trial_i).times = fixData(trial_i).times - psthParams.psthPre;
end
spikeData{1}{1}{1} = fixData;

% Once all the spikeData is collected, generate a PSTH from it.
psthEmptyByStim = initNestedCellArray(spikeData, 'zeros', [1, 1]);
[psthByStim, psthErrByStim] = calcStimPSTH(spikeData, psthEmptyByStim, 1, psthParams, spikeAlignParams);

% Generate Period shading
ITI = 0;%psthParams.ITI - 100;
pre = psthParams.psthPre - ITI;
ImDur = psthParams.psthImDur;
post = psthParams.psthPost;

totalTime =  ITI + pre + ImDur + post;
% peakValues = round(max(psthByStim{1}{1}(:)));
peakValues = 20;

periodImage = zeros(peakValues, totalTime, 3);
periodWindows = [1 ITI; ITI+1 ITI+pre; ITI+pre+1 ITI+pre+500; ITI+pre+500 + 1 ITI+pre+ImDur; ITI+pre+ImDur + 1 ITI+pre+ImDur+post];
periodColor = [0.7 0.7 0.7;...
                0.8 0.6 0.6;...
                0.4 0.4 0.8;...
                0.6 0.8 0.6;...
                0.4 0.4 0.4;];
unitColor = [0.8 0 0; ...
              0 0 0.8;...
              0 0.6 0; ...
              0 0 0;];

for epoch_i = 1:length(periodWindows)
  for layer_i = 1:3
    periodImage(:, periodWindows(epoch_i,1):periodWindows(epoch_i,2), layer_i) = deal(periodColor(epoch_i,layer_i));
  end
end

% Plotting
figTitle = sprintf('Units with Epoch selective activity');
figH = figure('Name', figTitle, 'NumberTitle','off', 'Units', 'normalized', 'Position', [0.3443 0.1176 0.4937 0.7185]);

psthParams.lineProps.col = mat2cell(unitColor, ones(1, size(unitColor, 1)), 3);
psthAxes = subplot(3,1,1);
pImageHand = image(periodImage);
pImageHand.XData = [-psthParams.psthPre psthParams.psthImDur + psthParams.psthPost];
set(psthAxes,'Ydir','normal')

% plotPSTH(psthByStim{1}{1}, psthErrByStim{1}{1}, psthAxes, psthParams, 'line', figTitle, unitName);
plotPSTH(psthByStim{1}{1}, [], psthAxes, psthParams, 'line', figTitle, unitName);

% Add text labels
peakVals = max(psthByStim{1}{1}, [], 2);
for label_i = 1:length(peakVals)
  xVal = epochLabelXVals(label_i);
  yVal = peakValues - 1;
  epochLabel = strsplit(epochNames{label_i}, ' ');
  textCount = 0;
  for ii = 1:length(epochLabel)
    text(xVal, yVal + textCount, epochLabel{ii}, 'FontWeight', 'bold');
    textCount = textCount - 1.5; 
  end
end

rasterAxes = subplot(3,1,2:3);
raster(spikeData, unitName, psthParams, 0, 1, 1, unitColor);

% Clean up Legends
psthAxes.FontSize = 10;
rasterAxes.FontSize = 15;
legendH = findobj(rasterAxes.Parent.Children, 'Type', 'Legend');
delete(legendH(1));

% Remove some labels, align things well.
psthAxes.XLabel.String = '';
psthAxes.XTickLabel = '';
psthAxes.Units = 'normalized';
rasterAxes.Units = 'normalized';
rasterAxes.Position(1) = psthAxes.Position(1);
rasterAxes.Position(3) = psthAxes.Position(3);
rasterAxes.Position(2) = psthAxes.Position(2) - rasterAxes.Position(4);

% Export the figure
figStruct.saveFig = 1;      % save the figure in its output directory.           
figStruct.closeFig = 0;     % close the figure once it is saved
figStruct.exportFig = 0;    % export figure using export_fig.
figStruct.exportFigSvg = 1;    % export figure using export_fig.
figStruct.saveFigData = 0;  % save data with the figure.
figStruct.noOverWrite = 1;  % If a figure is already there, don't make it again.

saveFigure(figOutDir, 'epochUnits_PSTHandRaster', [], figStruct, [])




function figh = plot_per_label_accuracy(analysisBatch, analysisStruct, params)
% a function which generates a plot containing:
% - single trace per label used in the decoder
% - a mean trace of the overall accuracy
% - vertical lines at important events (denoted below)
% - Relabels X-axis (with numbers below).
% Inputs are
% - decoding_results, ds generated by neural decoding toolbox
% - analysisStruct, generated by k_aid scripts.

% Load the outputs from these analysisStructs
analysisOutputStructs = cellfun(@(x) load(x), analysisBatch);
analysisOutputStructs = [analysisOutputStructs.analysisStruct];
decodingResultsStructs = {analysisOutputStructs.decoding_results_file}';
decodingResultsStructs = cellfun(@(x) load(x), decodingResultsStructs);
decodingResultsStructs = [decodingResultsStructs.decoding_results];

% Extract clear variables from inputs.
if length(decodingResultsStructs) ~= 1
  decoding_results_tmp = decodingResultsStructs(~[analysisOutputStructs.shuffle_ds]);
else
  decoding_results_tmp = decodingResultsStructs;
end

binParams = decoding_results_tmp.DS_PARAMETERS.binned_site_info.binning_parameters;
labels = decoding_results_tmp.DS_PARAMETERS.label_names_to_use; % Find Labels

if size(labels, 2) > 1
  labels = labels';
end

decoding_type = 'Zero';
switch decoding_type
  case 'Zero'
    decoding_data = [decodingResultsStructs.ZERO_ONE_LOSS_RESULTS];
  case 'Normalized'
    decoding_data = [decodingResultsStructs.NORMALIZED_RANK_RESULTS];
  case 'ROC_AUC'
    decoding_data = [decodingResultsStructs.ROC_AUC_RESULTS];
end

% Set up the other variables
binArray = 1:length(binParams.the_bin_start_times);

points_to_label = params.plotParams.points_to_label;
points_for_lines = params.plotParams.points_for_lines;

p_val_threshold = params.p_val_threshold;
sig_bins = params.sig_bins;
shift = binParams.alignment_event_time;                % The time before the stimulus the data extends to.
bin_start_times = binParams.the_bin_start_times - (shift + 1);   % 

plotP = params.plot_per_label_acc;
justMean = plotP.justMean;
plotMean = plotP.plotMean;
chanceAtBottom = plotP.chanceAtBottom;
plotError = plotP.plotError;
plotEachLabel = plotP.plotEachLabel;
groupNames = plotP.groupNames;
groups = plotP.groups;
sig_bar_pos = plotP.sig_bar_pos;
sig_color = plotP.sig_color;

% Identify what the x axis will represent.
binLabel = plotP.binLabel;
switch binLabel
  case 'End'
    offSet = params.binWidth;
  case 'Middle'
    offSet = params.binWidth/2;
end
the_bin_labels = bin_start_times + offSet;

bins_to_label = interp1(the_bin_labels, binArray, points_to_label);
x_for_lines = interp1(the_bin_labels, binArray, points_for_lines);

points_to_label = points_to_label(~isnan(bins_to_label));
bins_to_label = bins_to_label(~isnan(bins_to_label));

xMin = round(interp1(the_bin_labels, binArray, points_to_label(1) - 50));


if isnan(x_for_lines(1))
  x_for_lines(1) = binArray(1);
end

if isnan(x_for_lines(2))
  x_for_lines(2) = binArray(end);
end

if isnan(xMin)
  xMin = binArray(1);
end

assert(~any(isnan([bins_to_label x_for_lines xMin])), 'NaNs on plotting variables');

% Extract results
% if length(decoding_data) ~= 1
  % Initialize arrays
  gridSize = size(decoding_data(1).mean_decoding_results);
  if analysisOutputStructs(1).crossTempDecode
    diagInd = logical(eye(gridSize(1)));
  else
    diagInd = true(gridSize);
  end
  correctLineMeanStack = nan(gridSize(1), length(decoding_data));
  correctLineSTDStack = nan(gridSize(1), length(decoding_data));
  
  % Cycle through and collect needed data
  for ii = 1:length(decoding_data)
    correctLineMeanStack(:, ii) = decoding_data(ii).mean_decoding_results(diagInd);
    correctLineSTDStack(:, ii) = decoding_data(ii).stdev.over_CVs_combined_over_resamples(diagInd);
  end
  
  correctLineMean = mean(correctLineMeanStack, 2);
  correctLineSTD = mean(correctLineSTDStack, 2);
  
% end

if size(correctLineMean,1) == size(correctLineMean, 2)
  % This means Cross temporal decoding is on.
  correctLineMean = correctLineMean(logical(eye(size(correctLineMean, 1))));
end

% Extract a per label correct percentage for each label.
tmp = [decoding_data.confusion_matrix_results];
if length(tmp) ~= 1
  confMat = cat(4, tmp.confusion_matrix);
else
  confMat = tmp.confusion_matrix;
end

confMat = confMat./sum(confMat,1); % Normalize
confMat = mean(confMat, 4);
correctLineStack = nan(size(confMat,1), size(confMat,3));

for bin_i = 1:size(correctLineStack, 2)
  infoTmp = confMat(:, :, bin_i);
  correctLineStack(:, bin_i) = infoTmp(logical(eye(length(infoTmp))));
end

% Check if plotting error was requested, despite not having many lines.
if size(correctLineStack, 1) == 1 && plotError
  warning('There is only 1 decoding result, no error line to plot')
  plotError = 0;
end

if length(labels) == 2
  warning('Binary Decoder forced to plotting of just Mean')
  justMean = 1;
end

figTitle = sprintf('Decoding for %s', analysisStruct.plotTitle);
figh = figure('Name', figTitle, 'units', 'normalized', 'outerposition',[.1 .1 0.8 0.8]);
CoreLabel = strsplit(analysisStruct.plotTitle, ' ');
plotTitle = sprintf('Decoding for %s', CoreLabel{2});

axesh = axes(figh);
axesh.FontSize = 20;
title(plotTitle)
hold on
linePlotHandles = gobjects(length(labels) + 2, 1)';

% Matrix for line types.
switch analysisStruct.labelFile
  case 'cat'
    colorNum = [1 0 0; 0.8 0 0; 0.6 0 0; 0.4 0 0; 0 1 0; 0 0.8 0; 0 0 1; 0 0 0.6;]; % Segment colors;
    colorLabel = {'chasing', 'fighting', 'mounting', 'grooming', 'goalDirected', 'idle', 'objects', 'scene'}; % Segment colors;
    
    % Sort them
    [~, sortInd] = ismember(colorLabel, labels);
    colors = colorNum;
    correctLineStack = correctLineStack(sortInd ,:);
    labels = labels(sortInd);
    lineStyle = repmat({'-'}, [8,1]);
    lineWidth = 3;

  case {'stimSetA', 'stimSetB'}
    
    colorNum = [1 0 0; 0.8 0 0; 0.6 0 0; 0.4 0 0; 0 1 0; 0 0.8 0; 0 0 1; 0 0 0.6;]; % Segment colors;
    colorLabel = {'chasing', 'fighting', 'mounting', 'grooming', 'goalDirected', 'idle', 'objects', 'scene'}; % Segment colors;
    
    % Resort labels, data
    [labels, resortInd] = resortAndTrimLabels(labels); % Function resorts into order seen in colorLabel
    correctLineStack = correctLineStack(resortInd ,:);
        
    % Color and line style index
    colorIndex = reshape([1:8; 1:8], [16, 1]);
    lineStyle = repmat({'-'; ':'}, [8,1]);
    colors = colorNum(colorIndex,:);
    lineWidth = 2;

end

if justMean
  % Initialize this if only plotting.
  plotLabels = {};
else
  % Plot each curve individually.
  for ii = 1:size(correctLineStack, 1)
    linePlotHandles(ii) = plot(1:length(correctLineStack(ii,:)), correctLineStack(ii,:), 'linewidth', 3, 'color', colors(ii,:), 'lineStyle', lineStyle{ii});
  end
  plotLabels = labels;
end

% Plot the mean and chance
sigBarLabel = sprintf('Significant regions (>%s%%)', num2str((1 - p_val_threshold) * 100));
if plotMean || justMean
  linePlotHandles(length(plotLabels) + 1) = plot(correctLineMean, 'linewidth', 6, 'color', 'k', 'lineStyle', '-');
  linePlotHandles(length(plotLabels) + 1).Tag = 'All Label Mean';
  linePlotHandles(length(plotLabels) + 2) = plot(xlim(), [1/length(labels) 1/length(labels)], 'linewidth', 3, 'color', 'k', 'lineStyle', '--');
  allLabels = [plotLabels; 'All Label Mean'; 'Theoretical chance'; sigBarLabel];
else
  linePlotHandles(length(plotLabels) + 2) = plot(xlim(), [1/length(labels) 1/length(labels)], 'linewidth', 3, 'color', 'k');
  allLabels = [plotLabels; 'Theoretical chance'; sigBarLabel];
end

if strcmp(sig_bar_pos, 'top')
   % Top 10%
   yLims = [axesh.Position(2)+axesh.Position(4) - axesh.Position(4)/10, axesh.Position(4)/10];
else
   yLims = [axesh.Position(2), axesh.Position(4)/10];
end

decodingAx = gca;

% Add the legend
% legend(linePlotHandles, allLabels, 'AutoUpdate', 'off', 'location', 'northeastoutside')
legH = legend(allLabels, 'AutoUpdate', 'off', 'location', 'northeastoutside');
legH.FontSize = 14;

ylabel('Decoding Accuracy')
xlabel(sprintf('Bin %s time', binLabel))
xticks(bins_to_label);
xticklabels(points_to_label);

% Plot Significant p values
[sigBarImgAxs, sigBarHands] = add_bars_to_plots([], [], {sig_bins(xMin:end)}, sig_color, {sig_bins(xMin:end)}, yLims);
linePlotHandles = [linePlotHandles, sigBarHands];

if chanceAtBottom && min(linePlotHandles(1).YData) > 0.2
  ylim([floor(min(correctLineStack(:))*10)/10, 1]);
else
  ylim([0, 1]);
end
xlim([xMin, length(the_bin_labels)]);
for ii = 1:length(x_for_lines)
  plot([x_for_lines(ii), x_for_lines(ii)], ylim(), 'linewidth', 4, 'color', [0.2, 0.2, 0.2])
end

% Re-adjust significance bar to align with slightly changed size of plot.
for ii = 1:length(sigBarImgAxs)
  sigBarImgAxs(ii).Position(1:3) = decodingAx.Position(1:3);
  sigBarImgAxs(ii).Position(2) = sigBarImgAxs(ii).Position(2)*1.1;
end

% Save the figure
figData.lineTraces = correctLineMean;
figData.lineSTD = correctLineSTD;
figData.meanTrace = mean(correctLineMean);
figData.sig_bins = sig_bins;
figData.analysisStruct = analysisStruct; % Save the analysisStruct with the figure
figData.the_bin_start_times = bin_start_times;
figData.the_bin_labels = the_bin_labels;

if justMean
  switchString = 'justMean';
elseif plotMean == 0
  switchString = 'AllTrace';
else
  switchString = 'AllTrace + Mean';
end

saveFigure(analysisStruct.plotOutDir, sprintf('1. %s - %s', figTitle, switchString), figData, params.figStruct, [])
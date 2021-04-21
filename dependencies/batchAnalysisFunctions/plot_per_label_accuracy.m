function figh = plot_per_label_accuracy(decoding_results, ds, analysisStruct, params)
% a function which generates a plot containing:
% - single trace per label used in the decoder
% - a mean trace of the overall accuracy
% - vertical lines at important events (denoted below)
% - Relabels X-axis (with numbers below).
% Inputs are
% - decoding_results, ds generated by neural decoding toolbox
% - analysisStruct, generated by k_aid app.

points_to_label = params.plotParams.points_to_label;
points_for_lines = params.plotParams.points_for_lines;
the_bin_start_times = 1:params.stepSize:(params.end_time - params.binWidth  + 1);
p_val_threshold = params.p_val_threshold;
sig_bins = params.sig_bins;
sig_color = params.plot_per_label_acc.sig_color;
shift = params.plotParams.shift;

justMean = params.plot_per_label_acc.justMean;
plotMean = params.plot_per_label_acc.plotMean;
chanceAtBottom = params.plot_per_label_acc.chanceAtBottom;
plotError = params.plot_per_label_acc.plotError;
plotEachLabel = params.plot_per_label_acc.plotEachLabel;
groupNames = params.plot_per_label_acc.groupNames;
groups = params.plot_per_label_acc.groups;
sig_bar_pos = params.plot_per_label_acc.sig_bar_pos;

figStruct = params.figStruct;

figTitle = sprintf('Per Label accuracy trace for %s', analysisStruct.plotTitle);
figh = figure('Name', figTitle, 'units', 'normalized', 'outerposition',[.1 .1 0.8 0.8]);

axesh = axes(figh);
axesh.FontSize = 24;
title(figTitle)

hold on

% Find Unique Labels, extract them from confusionMatrix, plot them
if isempty(ds.label_names_to_label_numbers_mapping)
  labels = ds.label_names_to_use;
else
  labels = ds.label_names_to_label_numbers_mapping(ds.label_names_to_use);
end
if size(labels, 2) > 1
  labels = labels';
end

% Cycle through decoding results
confusionMatSize = size(decoding_results{1}.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix);
confusionMat = nan([confusionMatSize, length(decoding_results)]);

% Create trace w/ mean and std - * Need to explore this structure further
% to understand error bars.
% cvCount = size(decoding_results{1}.ZERO_ONE_LOSS_RESULTS.decoding_results, 1);
% meanTracePerCat = squeeze(mean(decoding_results{1}.ZERO_ONE_LOSS_RESULTS.decoding_results, 1));
% stdTracePerCat = decoding_results{1}.ZERO_ONE_LOSS_RESULTS.stdev.over_CVs_combined_over_resamples;

for dec_i = 1:length(decoding_results)
  confusionMatTmp = decoding_results{dec_i}.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix;
  confusionMat(:,:,:, dec_i) = confusionMatTmp ./ sum(confusionMatTmp, 1);
end
correctLineStack = zeros([length(labels), confusionMatSize(3), length(decoding_results)]);

% Collect each curve
for dec_i = 1:length(decoding_results)
  for jj = 1:length(labels)
    correctTrace = squeeze(confusionMat(jj,jj,:, dec_i))';
    correctLineStack(jj,:, dec_i) = correctTrace;
  end
end

% Create the plotting results
correctLineMean = mean(correctLineStack, 3);
correctLineSTD = std(correctLineStack, 0, 3);

% Check if plotting error was requested, despite not having many lines.
if size(correctLineStack, 3) == 1 && plotError
  warning('There is only 1 decoding result, no error line to plot')
  plotError = 0;
end

if length(labels) == 2
  warning('Binary Decoder forced to plotting of just Mean')
  justMean = 1;
end

% Create the traces and labels to plot
% - Either just the mean
% - the mean across individual preset groups
% - Individual traces 

if justMean
  % Initialize this if only plotting.
  plotLabels = {};
elseif ~plotEachLabel
  
  % Plot means of preset groups.
  groupsCount = length(groupNames);
  linePlotHandles = gobjects(groupsCount + 2, 1)';
  meanTraces2Plot = zeros(groupsCount, size(correctLineMean,2));
  meanErr2Plot = zeros(groupsCount, size(correctLineSTD,2));
  meanTraceLabels = cell(1, groupsCount);
  
  for group_i = 1:length(groups2Plot);
    % Find the traces to add
    [groupPresent, groupPresentInd] = intersect(labels, groups{group_i});
    % Combine them
    meanTraces2Plot(group_i,:) = mean(correctLineMean(groupPresentInd, :, :));
    meanErr2Plot(group_i,:) = mean(correctLineSTD(groupPresentInd, :, :));
    meanTraceLabels{group_i} = strjoin(groupPresent,' ');
  end
  
  %Plot them
  if plotError
    mseb(1:size(meanTraces2Plot,2), meanTraces2Plot, meanErr2Plot);
  else
    linePlotHandles(1:length(groups2Plot)) = plot(meanTraces2Plot', 'linewidth', 3);
  end
  plotLabels = meanTraceLabels';
  
else
  % Plot each curve individually.
  linePlotHandles = gobjects(length(labels) + 2, 1)';
  
  if plotError
    mseb(1:size(correctLineMean, 2), correctLineMean, correctLineSTD);
  else
    linePlotHandles(1:end-2) = plot(correctLineMean', 'linewidth', 3);
  end
  
  plotLabels = labels;
  
end

% Plot the mean and chance
sigBarLabel = sprintf('Significant regions (>%s%%)', num2str((1 - p_val_threshold) * 100));
if plotMean || justMean
  linePlotHandles(length(plotLabels) + 1) = plot(mean(correctLineMean), 'linewidth', 5, 'color', 'k');
  linePlotHandles(length(plotLabels) + 1).Tag = 'All Label Mean';
  linePlotHandles(length(plotLabels) + 2) = plot(xlim(), [1/length(labels) 1/length(labels)], 'linewidth', 3, 'color', 'b');
  allLabels = [plotLabels; 'All Label Mean'; 'Theoretical chance'; sigBarLabel];
else
  linePlotHandles(length(plotLabels) + 2) = plot(xlim(), [1/length(labels) 1/length(labels)], 'linewidth', 3, 'color', 'b');
  allLabels = [plotLabels; 'Theoretical chance'; sigBarLabel];
end

if strcmp(sig_bar_pos, 'top')
   % Top 10%
   yLims = [axesh.Position(2)+axesh.Position(4) - axesh.Position(4)/10, axesh.Position(4)/10];
else
   yLims = [axesh.Position(2), axesh.Position(4)/10];
end

decodingAx = gca;

% Plot Significant p values
[sigBarImgAxs, sigBarHands] = add_bars_to_plots([], [], {sig_bins}, sig_color, {sig_bins}, yLims);

linePlotHandles = [linePlotHandles, sigBarHands];

% Add the legend
% legend(linePlotHandles, allLabels, 'AutoUpdate', 'off', 'location', 'northeastoutside')
legend(allLabels, 'AutoUpdate', 'off', 'location', 'northeastoutside')

% Label axes correctly
the_bin_start_times_shift = the_bin_start_times - shift;
bins_to_label = interp1(the_bin_start_times_shift, 1:length(the_bin_start_times_shift), points_to_label);
x_for_lines = interp1(the_bin_start_times_shift, 1:length(the_bin_start_times_shift), points_for_lines);

ylabel('Decoding Accuracy')
xlabel('Bin start time')
xticks(bins_to_label);
xticklabels(points_to_label);

if chanceAtBottom && min(linePlotHandles(1).YData) > 0.2
  ylim([floor(min(correctLineStack(:))*10)/10, 1]);
else
  ylim([0, 1]);
end
xlim([1, length(the_bin_start_times)]);
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
figData.analysisStruct = analysisStruct; % Save the analysisStruct with the figure

saveFigure(analysisStruct.plotOutDir, ['1. ' figTitle], figData, figStruct, [])
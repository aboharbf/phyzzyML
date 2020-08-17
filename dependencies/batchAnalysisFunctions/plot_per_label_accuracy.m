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

figTitle = sprintf('Per Label accuracy trace for %s', analysisStruct.plotTitle);
figh = figure('Name', figTitle, 'units', 'normalized', 'outerposition',[.1 .1 0.8 0.8]);
axesh = axes(figh);
axesh.FontSize = 16;

title(figTitle)
hold on

% Find Unique Labels, extract them from confusionMatrix, plot them
labels = ds.label_names_to_label_numbers_mapping(ds.label_names_to_use);
if size(labels, 2) > 1
  labels = labels';
end
confusionMat = decoding_results.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix;
confusionMat = confusionMat ./ sum(confusionMat, 1);
correctLineStack = zeros(length(labels), size(confusionMat,3));
for jj = 1:length(labels)
  correctTrace = squeeze(confusionMat(jj,jj,:))';
  correctLineStack(jj,:) = correctTrace;
  plot(correctTrace, 'linewidth', 3)
end

% Plot the mean and chance, Add the legend
plot(mean(correctLineStack), 'linewidth', 5, 'color', 'r')
plot(xlim(), [1/length(labels) 1/length(labels)], 'linewidth', 3, 'color', 'b')
legend([labels; {'All Label Mean'}; {'Theoretical chance'}], 'AutoUpdate', 'off', 'location', 'northeastoutside')

% Label the axes, calculate actual bin starts
ylabel('Decoding Accuracy')
xlabel('Bin start time')

% Label x axes
the_bin_start_times_shift = the_bin_start_times - 800;
bins_to_label = interp1(the_bin_start_times_shift, 1:length(the_bin_start_times_shift), points_to_label);
x_for_lines = interp1(the_bin_start_times_shift, 1:length(the_bin_start_times_shift), points_for_lines);
xticks(bins_to_label);
xticklabels(points_to_label);
ylim([0, 1]);
for ii = 1:length(x_for_lines)
  plot([x_for_lines(ii), x_for_lines(ii)], ylim(), 'linewidth', 4, 'color', 'k')
end


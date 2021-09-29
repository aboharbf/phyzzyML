function generate_TCT_plot(analysisStruct, saved_results_name, params)
% Code which generates TCT matrix using neural decoding toolbox, and
% modifies axes to better delinates important times in trial.

% use core code to generate TCT matrix
plot_obj = plot_standard_results_TCT_object(analysisStruct.decoding_results_file);
plot_obj.saved_results_structure_name = saved_results_name;
plot_obj.display_TCT_movie = 0;
plot_obj.TCT_figure_number = length(findobj('type','figure')) + 1; % Open a new figure.
plot_obj.plot_results;  % plot the TCT matrix and a movie showing if information is coded by a dynamic population code

% Figure titles
params.figTitle = sprintf('Cross Temporal Decoding of %s', analysisStruct.plotTitle);
if params.addTCTSigShading
  sigStr = sprintf(', %s%% threshold', num2str(100 - (params.p_val_threshold * 100)));
  params.figTitle = horzcat([params.figTitle, sigStr]);
end

% Generate a title
title(params.figTitle);

% Generate vectors for new X and Y axes
tmp = load(analysisStruct.decoding_results_file);
binningParams = tmp.decoding_results.DS_PARAMETERS.binned_site_info.binning_parameters;

points_to_label = params.plotParams.points_to_label;
points_for_lines = params.plotParams.points_for_lines;
the_bin_start_times = binningParams.the_bin_start_times - 1;
the_bin_start_times_shift = the_bin_start_times - binningParams.alignment_event_time;

% Identify what the x axis will represent.
binLabel = params.plot_per_label_acc.binLabel;
switch binLabel
  case 'End'
    offSet = params.binWidth;
  case 'Middle'
    offSet = params.binWidth/2;
end
the_bin_start_times_shift = the_bin_start_times_shift + offSet;

bins_to_label = interp1(the_bin_start_times_shift, 1:length(the_bin_start_times_shift), points_to_label);
lines_to_add = interp1(the_bin_start_times_shift, 1:length(the_bin_start_times_shift), points_for_lines);

% Change x axis, set to manual
xticks(bins_to_label);
xticklabels(points_to_label);
xlim(xlim())

% Change y axis, set to manual
yticks(bins_to_label);
yticklabels(points_to_label);
ylim(ylim())
hold on

% Add vertical and horizontal line
for line_i = 1:length(lines_to_add)
  lineCoord = [lines_to_add(line_i), lines_to_add(line_i)];
  plot(lineCoord, ylim(), 'color', 'k')
  plot(xlim(), lineCoord, 'color', 'k')
end

% Make the correct size, increase font.
j = gcf;
j.Units = 'Normalized';
j.Position = [0 0.05 0.55 0.7];
ax = gca;
ax.FontSize = 14;
xlabel('Test time (ms)', 'FontSize', 20)
ylabel('Train time (ms)', 'FontSize', 20)

% Significance testing
if params.addTCTSigShading
  sigImg = (ax.Children(end).CData) > (params.decoding_threshold * 100);
  sigImg = double(sigImg);
  
  if params.removeSmallPatches
    sigImgPatches = bwconncomp(sigImg);
    patchSizes = cellfun('length', sigImgPatches.PixelIdxList)';
    patchRemoveInd = patchSizes < params.removeSmallPatch_cutOff;
    patchesRemove = sigImgPatches.PixelIdxList(patchRemoveInd);
    patchRemoveInd = vertcat(patchesRemove{:});
    sigImg(patchRemoveInd) = 0;
  end
  
  sigImg(sigImg == 0) = 0.6;
  ax.Children(end).AlphaData = sigImg;
end

% Save Figure
saveFigure(analysisStruct.plotOutDir, ['2. ' params.figTitle], analysisStruct, params.figStruct, [])

end
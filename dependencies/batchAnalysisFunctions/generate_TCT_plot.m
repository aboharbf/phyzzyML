function generate_TCT_plot(analysisStruct, save_file_name, saved_results_name, params)
% Code which generates TCT matrix using neural decoding toolbox, and
% modifies axes to better delinates important times in trial.

params.figTitle = sprintf('Cross Temporal Decoding of %s', analysisStruct.plotTitle);
if params.addTCTSigShading
  sigStr = sprintf(', %s%% threshold', num2str(100 - (params.p_val_threshold * 100)));
  params.figTitle = horzcat([params.figTitle, sigStr]);
end

% use core code to generate TCT matrix
plot_obj = plot_standard_results_TCT_object(save_file_name);
plot_obj.saved_results_structure_name = saved_results_name;
plot_obj.display_TCT_movie = 0;
plot_obj.TCT_figure_number = length(findobj('type','figure')) + 1; % Open a new figure.
plot_obj.plot_results;  % plot the TCT matrix and a movie showing if information is coded by a dynamic population code

% Generate a title
title(params.figTitle);

% Generate vectors for new X and Y axes
points_to_label = params.plotParams.points_to_label;
points_for_lines = params.plotParams.points_for_lines;
the_bin_start_times = 1:params.stepSize:(params.end_time - params.binWidth  + 1);
the_bin_start_times_shift = the_bin_start_times - params.plotParams.shift;

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
saveFigure(pFolder, ['2. ' params.figTitle], analysisStruct, params.figStruct, [])

end
function figh = plot_per_label_accuracy_combined(decoding_dir_list, params)
% a function which generates a plot containing:
% - single trace per label used in the decoder
% - a mean trace of the overall accuracy
% - vertical lines at important events (denoted below)
% - Relabels X-axis (with numbers below).
% Inputs:
% - decoding_dir_list, a set of directories containing the decoding
% results to be combined.
% - analysisStruct, generated by k_aid_generate_analyses.
% - params, details below

% Figure settings

figStruct.saveFig = 1;      % save the figure in its output directory.
figStruct.closeFig = 0;     % close the figure once it is saved
figStruct.exportFig = 1;    % export figure using export_fig.
figStruct.saveFigData = 1;  % save data with the figure.
figStruct.noOverWrite = 0;  % If a figure is already there, don't make it again.

fontPlotSize = 21;

% Consistent variables across analyses
endTimeFile = load('D:\DataAnalysis\batchAnalysis\NeuralDecodingTB\NaturalSocial\max_zsc\NS_Broad_Categories\results_NS_catBroad_AllUnits_MUA_MCC\decoding_results_run1.mat');
end_time = endTimeFile.decoding_results.DS_PARAMETERS.binned_site_info.binning_parameters.end_time;
shift = endTimeFile.decoding_results.DS_PARAMETERS.binned_site_info.binning_parameters.alignment_event_time;
binWidth = 150;
stepSize = 50;
points_to_label = [-300, 0, 500, 1000, 1500, 2000, 2500, 3000];
points_for_lines = [0, 2800];

% Info which aids in plotting - Legend is created as text hovering over
% traces - The info below is the coordinate, the alignment, and color.
positionInfo = {{50, 0.9, 'left', [1 0 0]};
  {1600, 0.8, 'center', [0 0 1]};
  {2750, 0.9, 'right', [0 0 0]};
  {1600, 0.9, 'center', [0.5 0.5 0.5]};
  {50, 0.8, 'left', [0 0 1]};
  {2750, 0.8, 'right', [0.5 0.5 0.5]}};

points_for_text_x = unique(cellfun(@(x) x{1}, positionInfo)); % Unique entries in the list above, for the sake of converting into the proper numbers for plotting.

% Label axes correctly - Prepare data for it
offSet = binWidth;
the_bin_start_times = (1:stepSize:(end_time - binWidth  + 1)) + offSet;
binCount = length(the_bin_start_times);
the_bin_start_times_shift = the_bin_start_times - shift;
bins_to_label = interp1(the_bin_start_times_shift, 1:binCount, points_to_label);
bins_for_text = interp1(the_bin_start_times_shift, 1:binCount, points_for_text_x);
x_for_lines = interp1(the_bin_start_times_shift, 1:binCount, points_for_lines);
xMin = round(interp1(the_bin_start_times_shift, 1:binCount, points_to_label(1) - 100));


% Combined Plot Types - Below are the types of plots to make for each
% group. There are 4 of each.

positionPerLabel = {[1 2]}; % Position info for each Label
subSlicePlotTitle = {'Combined', 'Combined'}; % Labels that show up on plot
plotSubTypes2Combine = {'socIntSelAny','socIntSelAny'}; % The Identifier tag to be used to grab the correct file.

% Per Plot -
coreDirectory = {'D:\DataAnalysis\batchAnalysis\NeuralDecodingTB\NaturalSocial\max_zsc', 'D:\DataAnalysis\batchAnalysis\NeuralDecodingTB\headTurnCon\max_zsc'};
coreDirectoryParadigm = {'NS', 'HTC'};

% Iterate across
plotCoreNameTag = {'socInt'};
fileCoreName = {'broadCat', 'socInt', 'broadCat', 'stim'};

templateDecodingString = 'D:\DataAnalysis\batchAnalysis\NeuralDecodingTB\NaturalSocial\max_zsc\%s\results_NS_%s_AllUnits_MUA_MCC\decoding_results_run1.mat';


for par_i = 1%:length(coreDirectory)
  
  plotOutDir = fullfile(coreDirectory{par_i}, 'combo_plot_directory_HardCode');
  unitType = {'UnUS'};
  
  plotComboTitle = {'Social vs Non-Social Category classification'};
  plotCoreNames = {{'Non-Social Categories', 'Social Categories'}};
  plotCoreFileTag = {{'nonSocialCat', 'socialCat'}};
  subSlicePossibleTitles = {{'Non-Social Categories', 'Social Categories'}};

% Identify all data files for this paradigm.
  plotDataFiles = dir(fullfile(coreDirectory{par_i}, '**', '*_data.mat'));
  plotDataFiles = fullfile({plotDataFiles.folder}, {plotDataFiles.name})';
  
  for unit_i = 1:length(unitType)
    
    unitIndex = contains(plotDataFiles, unitType{unit_i});
    [unitCountPerTrace, meanTraceStack, sigTraceStack] = deal(cell(size(plotCoreNames)));
    figDataStructs = [];
    
    for group_i = 1:length(plotCoreNames)
      
      for plot_i = 1:length(plotCoreNames{group_i})
        
        % Identify the combined paradigm plot name, and identify the correct
        % files.
        parPlotName = strjoin([coreDirectoryParadigm(par_i), plotCoreNames{group_i}(plot_i)], ' ');
        sub_plot_ind = contains(plotDataFiles, parPlotName) & unitIndex & contains(plotDataFiles, plotSubTypes2Combine{plot_i});
        analysisDirName = strrep(parPlotName, ' ', '_');  % Used later for template - may lead to problems if combined types dont have same category count (this case they do).
        
        figDataStructs = [figDataStructs, load(plotDataFiles{find(sub_plot_ind,1)})];
        
      end
      
      % Collect Result from files
      figData = [figDataStructs.figData];
      analysisStructs = [figData.analysisStruct];
      
      % Collect variables from the data files
      unitCountPerTrace = cellfun('length', {analysisStructs.sites});
      meanTraceStack = vertcat(figData.meanTrace);
      sigTraceStack = vertcat(figData.sig_bins);
      
      % Initiate the plot
      figTitle = sprintf('Combined Accuracy traces for %s - %s %s', parPlotName, plotCoreNameTag{1}, unitType{unit_i});
      figh = figure('Name', figTitle, 'units', 'normalized', 'outerposition', [0.1042 0.1370 0.6130 0.6843]);
      axesh = axes(figh);
      
      % Plotting
      % Plot the significant slices
      sigData = deal(meanTraceStack);
      sigData(~sigTraceStack) = deal(nan);
      
      % Processing the trace file tags into labels
      hold on
      legendLabels = plotCoreNames{group_i};
      traceInfo = positionInfo(positionPerLabel{group_i});
      
      % Plot significant traces wider then non-significant
      for trace_i = 1:size(sigData,1)
        plot(meanTraceStack(trace_i, :), 'linewidth', 2, 'color', traceInfo{trace_i}{end});
        plot(sigData(trace_i,:), 'linewidth', 5, 'color', traceInfo{trace_i}{end});
      end
      
      % Check if the labels need to be shifted up
      plot_text_y = cellfun(@(x) x{2}, traceInfo);
      bins_for_text_plot = cellfun(@(x) bins_for_text(x{1} == points_for_text_x), traceInfo);
      if max(sigData(:)) > 0.8
        plot_text_y = plot_text_y + 0.05;
      end
      
      % Add legend
      for leg_i = 1:length(legendLabels)
        legendText = [legendLabels{leg_i} ' (' num2str(unitCountPerTrace(leg_i)) ')'];
        textHandle = text(bins_for_text_plot(leg_i), plot_text_y(leg_i), legendText, 'FontSize', fontPlotSize, 'color', traceInfo{leg_i}{end} , 'HorizontalAlignment', traceInfo{leg_i}{3}, 'FontWeight', 'bold');
      end
      
      ylabel('Decoding Accuracy')
      xlabel('Bin End time')
      xticks(bins_to_label);
      xticklabels(points_to_label);
      title(figTitle)
      
      ylimSize = ylim();
      ylim([ylimSize(1), 1]);
      
      % Calibrate the x axis
      xlim([xMin, length(the_bin_start_times)]);
      for ii = 1:length(x_for_lines)
        plot([x_for_lines(ii), x_for_lines(ii)], ylim(), 'linewidth', 4, 'color', [0.2, 0.2, 0.2])
      end
      axesh.FontSize = fontPlotSize;
      
      % Create the needed theoretical chance line
      templateStruct = load(sprintf(strrep(templateDecodingString, '\', '\\'), analysisDirName, plotCoreFileTag{group_i}{end}));
      objectCount = templateStruct.decoding_results.CV_PARAMETERS.num_unique_labels;
      plot(xlim(), [1/objectCount 1/objectCount], 'linewidth', 3, 'color', 'k');
      
      saveFigure(plotOutDir, figTitle, [], figStruct, [])
    end
  end
end

end
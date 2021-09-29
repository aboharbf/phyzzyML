function saveFigure(outDir, filename, figData, figStruct, footer)
% function which saves present figure object, along with its data,
% according to parameters defined in the figStruct. Also adds text to the
% lower left of a figure.
% Input Arguments are:
%   outDir - save path for the figure
%   filename - name of saved figure
%   figData - data to be saved w/ the figure.
%   figStruct - structure containing switches defined below for saving,
%   exporting, and data saving w/ figure.
%   footer - text to add to the bottom figure.

% Make the outDir if it doesn't exist
if ~exist(outDir,'dir')
  mkdir(outDir)
end

if ~isempty(footer)
  ax = axes('Position',[0 0 1 .05], 'Visible','off');
  text(ax, .025, .5, footer, 'fontsize', 8);
end
if figStruct.saveFig     
  savefig(fullfile(outDir, strcat(filename, '.fig')));
end
if figStruct.exportFig
  export_fig(fullfile(outDir, strcat(filename, '.png')),'-m1.2','-transparent','-opengl');
end
if figStruct.saveFigData && ~isempty(figData)
  save(fullfile(outDir,[filename,'_data.mat']),'figData');
end
if figStruct.closeFig
  close()
end

end


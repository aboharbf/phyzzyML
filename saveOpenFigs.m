function saveOpenFigs()

h =  findobj('type','figure');
saveDir = 'C:\OneDrive\Lab\Documents\DisserationFigs\saccadePostEx';
figCoreName = 'uniqueCellEx2';

figStruct.saveFig = 1;      % save the figure in its output directory.           
figStruct.closeFig = 0;     % close the figure once it is saved
figStruct.exportFig = 1;    % export figure using export_fig.
figStruct.saveFigData = 0;  % save data with the figure.
figStruct.noOverWrite = 1;  % If a figure is already there, don't make it again.

for fig_i = 1:length(h)
  figure(h(fig_i));
  if ~isempty(h(fig_i).Name)
    saveFigure(saveDir, sprintf('%s - %s - %d', figCoreName, strrep(h(fig_i).Name, '/', ' of '), fig_i), [], figStruct, [])
  else
    saveFigure(saveDir, sprintf('%s - %d', figCoreName, fig_i), [], figStruct, [])
  end
end

winopen(saveDir);

end

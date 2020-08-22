function exportOpenFigures()
% A function which checks all the figures which are currently openned and
% applies the 'exportFig' function to them.

tmp = pwd;
addpath(genpath(fileparts(tmp)));
path2save = 'H:\Analyzed\batchAnalysis\NeuralDecodingTB';

figHandles = findobj('Type', 'figure');

for ii = 1:length(figHandles)
  set(0, 'currentfigure', figHandles(ii))
  figName = figHandles(ii).Name;
  if isempty(figName)
    figName = figHandles(ii).Children(2).Title.String;
  end
  export_fig(fullfile(path2save, strcat(figName, '.png')),'-m1.2','-transparent','-opengl');
  
end
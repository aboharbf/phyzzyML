function [ ] = saveFigure(outDir, filename, figData, figStruct, footer)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(footer)
  ax = axes('Position',[0 0 1 .05], 'Visible','off');
  text(ax, .025,.5,footer, 'fontsize',12);
end
if figStruct.saveFig     
  savefig(fullfile(outDir, strcat(filename, '.fig')));
end
if figStruct.exportFig
  export_fig(fullfile(outDir, strcat(filename, '.png')),'-m1.2','-transparent','-opengl');
end
if figStruct.saveFigData
  save(fullfile(outDir,[filename,'_data.mat']),'figData');
end
if figStruct.closeFig
  close()
end

end


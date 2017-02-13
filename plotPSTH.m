function [ psthAxes ] = plotPSTH( psthArray, psthAxes, psthPre, psthPost, psthImDur, plotType, psthTitle, ylabels, psthColormap )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

xrange= [-psthPre psthImDur+psthPost]; 
nrows = size(psthArray,1);
yrange= [1 nrows];
title(psthTitle); 
if strcmp(plotType,'color')
  caxis([min(min(psthArray)),max(max(psthArray))]);
  if isempty(psthAxes)
    psthAxes = imagesc(xrange, yrange, psthArray);
  else
    imagesc(psthAxes, xrange, yrange, psthArray);
  end
  h = colorbar; 
  ylabel(h,'Hz','FontSize',14);
  colormap(psthColormap);
  ylimits= ylim();
  yRange = ylimits(2) - ylimits(1);  %todo: remove egregious use of this variable twice, but with different caps
  hold on
  draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',4);
  draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',4);
  set(gca,'YTick',linspace(ylimits(1)+yRange/(2*nrows),ylimits(2)-yRange/(2*nrows),nrows),'YTicklabel',ylabels,...
    'box','off','TickDir','out','FontSize',14,'TickLength',[.012 .012],'LineWidth',1.3);
else
  error('psth line plot not yet implemented');
end
xlabel('Time from stimulus onset (ms)', 'FontSize',14);
hold off

Output.DEBUG('min of psth %s: %d\n',psthTitle,(min(min(psthArray))));
Output.DEBUG('max of psth %s: %d\n',psthTitle,(max(max(psthArray))));
end


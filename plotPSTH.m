function [ psthAxes,  colorBarh, vertLineHands, legendH] = plotPSTH(psthArray, psthError, psthAxes, psthParams, plotType, psthTitle, ylabels)
% Function which accepts data and generate a PSTH. Accepts the following
% inputs.
% psthArray: stimuli * bin in size.
% psthError: used when plotType is 'line' to created shaded region of
% error.
% psthAxes: a handle for axes. if not defined, gca are used.
% psthParams: contains the timing variables below, as well as lineprops for
% line plot PSTH and the colormap for the standard one.
% plotType: 'color' or 'line'. line uses mseb().
% psthTitle: the title above the axes.
% ylabels: a 'stimuli * 1' cell array, with the name for each stimulus.
% used for the legend of the line plot, and the y axis labels for the color
% PSTH.
% Outputs:
% handles for the psthAxes, colorBar, or vertical lines.

psthPre = psthParams.psthPre;
psthImDur = psthParams.psthImDur;
psthPost = psthParams.psthPost;

xrange= [-psthPre psthImDur+psthPost]; 
nrows = size(psthArray,1);
yaxis = [1 nrows];
switch plotType
  case 'color'
    try
      caxis([min(min(psthArray)),max(max(psthArray))]);
    catch
      disp('start error message');
      disp(psthArray);
      disp([min(min(psthArray)),max(max(psthArray))]);
      assert(0,'failed in plotPSTH');
      return
    end
    
    if isempty(psthAxes)
      psthAxes = imagesc(xrange, yaxis, psthArray);
    else
      imagesc(psthAxes, xrange, yaxis, psthArray);
    end
    colorBarh = colorbar;
    ylabel(colorBarh,'Firing Rate [Hz]','FontSize',14);
    if isfield(psthParams, 'colormap')
      colormap(psthParams.colormap);
    end
    ylimits= ylim();
    yRange = ylimits(2) - ylimits(1);
    set(gca,'YTick',linspace(ylimits(1)+yRange/(2*nrows),ylimits(2)-yRange/(2*nrows),nrows),'YTicklabel',ylabels,...
      'box','off','TickDir','out','FontSize',14,'TickLength',[.012 .012],'LineWidth',1.3);
  case 'line'
    colorBarh = [];
    if isempty(psthError)
      psthError = zeros(size(psthArray));
    end
    
    if isfield(psthParams, 'lineProps')
      lineProps = psthParams.lineProps;
    else
      lineProps = [];
    end
    
    if ~isempty(psthAxes)
      axes(psthAxes);
    end
    
    psthAxes = mseb(xrange(1):xrange(2), psthArray, psthError, lineProps);
    
    xlim(xrange);
    ylim(ylim()); % Shifts auto ylim to manual, preventing lines below from expanding the window.
    legendH = legend(ylabels, 'location', 'northeastoutside', 'AutoUpdate', 'off');
end

% Deliniate stimulus presentation period
hold on
if (psthPre+psthPost)/psthImDur > 20
  stimDurLineWidth = 0.1;
else
  stimDurLineWidth = 4;
end

vertLineColor = [0.5, 0.5, 0.5];
vertLineHands = gobjects(2,1);
vertLineHands(1) = line([0 0], ylim(),'Color',vertLineColor,'LineWidth',stimDurLineWidth);
vertLineHands(2) = line([psthImDur psthImDur], ylim(),'Color',vertLineColor,'LineWidth',stimDurLineWidth);

xlabel('Time from stimulus onset (ms)', 'FontSize',14);
title(psthTitle);
hold off

Output.DEBUG('min of psth %s: %d\n',psthTitle,(min(min(psthArray))));
Output.DEBUG('max of psth %s: %d\n',psthTitle,(max(max(psthArray))));

end


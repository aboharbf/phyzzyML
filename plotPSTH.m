function [psthAxes, colorBarh, vertLineHands, legendH] = plotPSTH(psthArray, psthError, psthAxes, psthParams, plotType, psthTitle, ylabels)
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

if isfield(psthParams, 'addLegend')
  addLegend = psthParams.addLegend;
else
  addLegend = true;
end

xrange= [-psthPre psthImDur+psthPost]; 
nrows = size(psthArray,1);
yaxis = [1 nrows];

switch plotType
  case 'color'
%     minVal = min(psthArray, [], 'all');
%     maxVal = max(psthArray, [], 'all');
%     try
%       caxis([minVal, maxVal]);
%     catch
%       error('Error using caxis, min = %d, max = %d', minVal, maxVal);
%     end
    
    if isempty(psthAxes)
      psthAxes = gca;
    end
    
    imagesc(psthAxes, xrange, yaxis, psthArray);
          
    psthAxes.YLim = [yaxis(1) - 0.5, yaxis(2) + 0.5];
    psthAxes.XLim = [ -psthPre, psthImDur + psthPost];
    
    colorBarh = colorbar(psthAxes);
    ylabel(colorBarh,'Firing Rate [Hz]','FontSize',14);
    if isfield(psthParams, 'colormap')
      colormap(psthParams.colormap);
    end
    ylimits = ylim(psthAxes);
    yRange  = ylimits(2) - ylimits(1);
    set(psthAxes,'YTick',linspace(ylimits(1)+yRange/(2*nrows),ylimits(2)-yRange/(2*nrows),nrows),'YTicklabel',ylabels,...
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
    else
      psthAxes = gca;
    end
    
    patchHandles = mseb(xrange(1):xrange(2), psthArray, psthError, lineProps);
    
    xlim(xrange);
    ylim(ylim()); % Shifts auto ylim to manual, preventing lines below from expanding the window.
    ylabel('Firing Rate (Hz)')

    % Assumption - if I've added more labels than there are lines, I want
    % the legend to include these labels for the sake of objects to be
    % generated later.
    if size(psthArray,1) < length(ylabels)
      hold on
      dummyLines = length(ylabels) - size(psthArray,1);
      dummyLinHands = gobjects(dummyLines,1);
      for ii = 1:dummyLines     
        dummyLinHands(ii) = plot(0, 0, 'color', 'k', 'linewidth', 0.5);
      end
    end
    if addLegend
      legendH = legend(ylabels, 'location', 'northeastoutside', 'AutoUpdate', 'off');
    end
    
end

% Deliniate stimulus presentation period
hold(psthAxes, 'on')
if (psthPre+psthPost)/psthImDur > 20
  stimDurLineWidth = 0.1;
else
  stimDurLineWidth = 4;
end

vertLineColor = [0.5, 0.5, 0.5];
vertLineHands = gobjects(2,1);
vertLineHands(1) = line(psthAxes, [0 0], ylim(psthAxes), 'Color', vertLineColor, 'LineWidth', stimDurLineWidth);
vertLineHands(2) = line(psthAxes, [psthImDur psthImDur], ylim(psthAxes),'Color',vertLineColor,'LineWidth',stimDurLineWidth);

xlabel(psthAxes, 'Time from stimulus onset (ms)', 'FontSize',14);
title(psthAxes, psthTitle);
hold(psthAxes, 'off')

Output.DEBUG('min of psth %s: %d\n',psthTitle,(min(min(psthArray))));
Output.DEBUG('max of psth %s: %d\n',psthTitle,(max(max(psthArray))));

end


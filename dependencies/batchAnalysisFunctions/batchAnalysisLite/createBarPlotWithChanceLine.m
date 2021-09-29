function figH = createBarPlotWithChanceLine(varargin)

if nargin == 6
  [barLabels, dataMat, alpha, unitCount, figTitle, legendEntries] = deal(varargin{:});
elseif nargin == 7
  [barLabels, dataMat, alpha, unitCount, figTitle, legendEntries, params] = deal(varargin{:});
end

% params contain colors - use the labels to select colors
if exist('params', 'var')
  [test, colorInd] = ismember(params.labels, legendEntries);
  colorMat = params.colors(test, :);
  colorInd = colorInd(~colorInd == 0);
  colorMat = colorMat(colorInd,:);
else
  colorMat = [];
end

figH = figure('Name', figTitle, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [.3 .3 .65 .77]);
X = categorical(barLabels);
X = reordercats(X,barLabels);
barh = bar(X, dataMat, 1);
ylimSize = ylim();

if length(barh) == 1
  barh.BarWidth = 0.5;
  legLoc = 'northeastoutside';
elseif length(legendEntries) > 4
  legLoc = 'northeastoutside';
else
  legLoc = 'northeast';
end

for bar_i = 1:length(barh)
  if ~isempty(colorMat)
    barh(bar_i).FaceColor = colorMat(bar_i,:);
  elseif all(barh(bar_i).FaceColor == [0 0 1])
    barh(bar_i).FaceColor = [0 0.25 1];
  end
  
  xtips1 = barh(bar_i).XEndPoints;
  ytips1 = barh(bar_i).YEndPoints;
  
  if bar_i ~= 1
%     diffY = abs(ytipsLast - ytips1);
%     if any(diffY < 1) && size(dataMat,1) >= 3
%       changeInd = diffY < 1;
%       ytips1(changeInd) = ytips1(changeInd) + ((ylimSize(2) * 0.05) );
%     end
  end
   
  unitTag = string(barh(bar_i).YData);
  labels2 = string(round(barh(bar_i).YData/unitCount,3) * 100);
  percentTag = strcat('(', labels2, ')');
  
  if ~isempty(legendEntries)
    barh(bar_i).DisplayName = legendEntries{bar_i};
  else
    barh(bar_i).DisplayName = 'Count';
  end
  
  textH = text(xtips1, ytips1, unitTag, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);
  textH = text(xtips1, ytips1, percentTag, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 14);
  ytipsLast = ytips1;
end

% Resize and shape the figure
YlimN = figH.Children.YLim;
figH.Children.YLim(2) = YlimN(2) * 1.2;

xlim(xlim());
hold on
if alpha ~= 0
  chanceUnitCount = round(sum(unitCount)*(alpha/2));
  chanceLine = plot([X(1) X(end)], [chanceUnitCount chanceUnitCount], 'Color', 'red', 'LineStyle', '--', 'linewidth', 3);
  chanceLine.DisplayName = sprintf('Chance (%s)', num2str(alpha/2));
end

figH.Children.FontSize = 16;

title(figTitle);
ylabel('Unit Count')
legend('Location', legLoc)

% Make sure the rightmost bar is not obscured with the legend. Prevent this
% by reserving the top 33% for the legend
barHeights = [barh.YData];
yLimits = ylim();
if any(barHeights >= yLimits(2) * 0.75)
  ylim([0 yLimits(2)*1.1])
end

end
function figH = createBarPlotWithChanceLine(barLabels, dataMat, alpha, unitCount, figTitle, legendEntries)

figH = figure('Name', figTitle, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [.3 .3 .5 .6]);

X = categorical(barLabels);
X = reordercats(X,barLabels);
barh = bar(X, dataMat, 1);
ylimSize = ylim();

if size(dataMat,2) == 1
  barh.BarWidth = 0.5;
end

for bar_i = 1:length(barh)
  xtips1 = barh(bar_i).XEndPoints;
  ytips1 = barh(bar_i).YEndPoints;
  
  if bar_i ~= 1
    diffY = abs(ytipsLast - ytips1);
    if any(diffY < 1) && size(dataMat,1) >= 3
      changeInd = diffY < 1;
      ytips1(changeInd) = ytips1(changeInd) + ((ylimSize(2) * 0.05) );
    end
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
  textH = text(xtips1, ytips1, percentTag, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10);
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
xlabel('Preferred Epoch')
ylabel('Unit Count')
legend()

end
function createBarPlotWithChanceLine(barLabels, dataMat, alpha, unitCount, figTitle, legendEntries)

figH = figure('Name', figTitle, 'NumberTitle','off','units','normalized', 'outerposition', [.3 .3 .5 .6]);

X = categorical(barLabels);
X = reordercats(X,barLabels);
barh = bar(X, dataMat);
for bar_i = 1:length(barh)
  xtips1 = barh(bar_i).XEndPoints;
  ytips1 = barh(bar_i).YEndPoints;
  labels1 = string(barh(bar_i).YData);
  labels2 = string(round(barh(bar_i).YData/unitCount,3));
  labels3 = strcat(labels1, '(', labels2, ')');
  
  if ~isempty(legendEntries)
    barh(bar_i).DisplayName = legendEntries{bar_i};
  else
    barh(bar_i).DisplayName = 'Count';
  end
  
  text(xtips1, ytips1, labels3, 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom')
end

% Resize and shape the figure
YlimN = figH.Children.YLim;
figH.Children.YLim(2) = YlimN(2) * 1.2;

xlim(xlim());
hold on
chanceUnitCount = round(sum(unitCount)*(alpha/2));

chanceLine = plot([X(1) X(end)], [chanceUnitCount chanceUnitCount], 'Color', 'red', 'LineStyle', '--', 'linewidth', 3);
chanceLine.DisplayName = sprintf('Chance (%s)', num2str(alpha/2));
figH.Children.FontSize = 16;

title(figTitle);
xlabel('Preferred Epoch')
ylabel('Unit Count')
legend()

end
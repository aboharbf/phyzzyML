function selectivityPerEventBarGraphs(selTableParadigm, paradigm, barPlotParams)

plotEye = 1;
selCheck = barPlotParams.selCheck;
UnitTypes = barPlotParams.UnitTypes;
UnitTypePlot = barPlotParams.UnitTypePlot;
colNamePoss = barPlotParams.colNamePoss;
colNamePlotAll = barPlotParams.colNamePlotAll;
alpha = 0.05;                                         % Alpha that the runs were done at.

% Generate bar plots showing total counts in each category
for unitType_i = 1:length(UnitTypes)
  % Filter the table per unit
  unitInd = contains(selTableParadigm.unitType, UnitTypes{unitType_i});
  selTableParadigmUnit = selTableParadigm(unitInd, :);
  
  if ~isempty(selTableParadigmUnit)
    % Plots 1 - fixation/saccade
    unitCount = sum(unitInd);
    chanceUnitCount = round(sum(unitCount) * alpha);
    saccSelCount = sum(selTableParadigmUnit.saccSel ~= 0);
    fixSelCount = sum(selTableParadigmUnit.fixationSel ~= 0);
    
    % Reward processing - take the rewardAbsent vec, and use the other reward
    % paradigm to fill in parts where that paradigm wasn't used.
    rewardVec = selTableParadigmUnit.subSel_rewardAbsent;
    rewardAbsMissing = isnan(rewardVec);
    rewardVec(rewardAbsMissing) = selTableParadigmUnit.subSel_reward(rewardAbsMissing);
    rewardVec = sum(rewardVec ~= 0);
    
    % If this row is present, and isn't all NaNs, use it.
    
    %   soccSacCount = sum(selTableParadigmUnit.saccSel_socInt ~= 0);
    %   events = {'Directional Saccade', 'Fixation', 'Reward', 'Saccades During SocInt'};
    %   dataMat = [saccSelCount; fixSelCount; rewardVec; soccSacCount];
    
    %   events = {'Directional Saccade', 'Fixation', 'Reward'};
    %   dataMat = [saccSelCount; fixSelCount; rewardVec];
    
%     events = {'Directional Saccade', 'Fixation', 'Reward'};
%     dataMat = [saccSelCount; fixSelCount; rewardVec];
    
    events = {'Fixation', 'Reward'};
    dataMat = [fixSelCount; rewardVec];
    
    % Plot
    figTitle = sprintf('%s activity selective for non-Stimulus Events during %s', UnitTypePlot{unitType_i}, paradigm);
    figH = figure('Name', figTitle, 'NumberTitle','off','units','normalized', 'outerposition', [.3 .3 .4 .6]);
    
    % Plot
    X = categorical(events);
    X = reordercats(X,events);
    barh = bar(X, dataMat);
    
    % Detail the Plot
    legend('Unit Count', 'Location','northeastoutside');
    for bar_i = 1:length(barh)
      xtips1 = barh(bar_i).XEndPoints;
      ytips1 = barh(bar_i).YEndPoints;
      
      labels1 = string(barh(bar_i).YData);
      labels2 = string(round(barh(bar_i).YData/unitCount,3));
      labels3 = strcat(labels1, '(', labels2, ')');
      
      text(xtips1, ytips1, labels3, 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom')
    end
    
    axH = findobj(figH.Children, 'Type', 'Axes');
    YlimN = axH.YLim;
    axH.YLim(2) = YlimN(2) * 1.2;
    xlim(xlim());
    
    figH.Children(2).FontSize = 16;
    xlabel('Event')
    ylabel('Unit Count')
    
    % Add line for chance number of units.
    hold on
    chanceLine = plot([X(1) X(end)], [chanceUnitCount chanceUnitCount], 'Color', 'red', 'LineStyle', '--', 'linewidth', 3);
    chanceLine.DisplayName = sprintf('Chance (%s)', num2str(alpha));
    
    %   title(figTitle);
    
    saveFigure(barPlotParams.outputDir, figTitle, [], barPlotParams.figStruct, [])
  end
end

end
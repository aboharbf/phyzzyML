function corrPlotSimple(selTable)
% Take a selTable, extract certain columns, and see if they are correlated.

% Reward processing - take the rewardAbsent vec, and use the other reward
% paradigm to fill in parts where that paradigm wasn't used.
rewardVec = selTable.subSel_rewardAbsent;
rewardAbsMissing = isnan(rewardVec);
rewardVec(rewardAbsMissing) = selTable.subSel_reward(rewardAbsMissing);

selTable.reward = rewardVec;

% variableNames = selTable.Properties.VariableNames;
variableNames = {'socIntSel_stimOnset', 'socIntSel_stimPres', 'reward', 'saccSel', 'fixationSel'};

correlationMat = selTable{:, variableNames};
correlationInfo = corr(correlationMat, 'Rows', 'pairwise');

% Create an array that create
indexForPlot = meshgrid(1:6, 1:6);
points2Correlate = cell(size(indexForPlot));

for row_i = 1:length(variableNames)
  for col_i = 1:length(variableNames)
    if row_i == col_i
      data2Store = selTable.(variableNames{row_i});
      
    else
      data2Store = [selTable.(variableNames{row_i}), selTable.(variableNames{col_i})];
      
    end    
    
    % keep
    keepInd = logical(sum(data2Store ~= 0, 2));
    
    % Store
    points2Correlate{row_i,col_i} = data2Store(keepInd,:);
    
    
  end
end

% Plot
figure();

for row_i = 1:length(variableNames)
  for col_i = 1:length(variableNames)
    % Plot index
    subplot(length(variableNames), length(variableNames), length(variableNames) * (row_i - 1) + col_i)
    
    if row_i == col_i
      histogram(points2Correlate{row_i, col_i}, 50);
    else
      scatH = scatter(points2Correlate{row_i, col_i}(:,1), points2Correlate{row_i, col_i}(:,2));
      title(sprintf('R^2 = %s', num2str(correlationInfo(row_i, col_i), 2)));
    end
    
    if row_i == length(variableNames)
      xlabel(variableNames{col_i})
    end
    
    if col_i == length(variableNames)
      ylabel(variableNames{col_i})
    end
    
  end
end

end
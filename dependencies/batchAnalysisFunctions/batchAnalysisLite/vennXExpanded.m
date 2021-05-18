function vennXExpanded(sigInd, figTitle, colNames)
% A wrapper for vennX, which makes it more straight forward to turn an
% array of logical values into the appropriate cross sections. 

[unitCount, colInCount] = size(sigInd);

if colInCount == 2
  sigIndTmp = nan(unitCount, 3);
  
  sigIndTmp(:,1) = sigInd(:,1) & ~sigInd(:,2); % A Count
  sigIndTmp(:,2) = sigInd(:,1) & sigInd(:,2);  % AB Count
  sigIndTmp(:,3) = ~sigInd(:,1) & sigInd(:,2); % B Count
  
  % Chance rates.
  chanceAB = (sum(sigInd(:,1))/unitCount) * (sum(sigInd(:,2))/unitCount) * unitCount;

  chanceCounts = round(chanceAB);
  chanceTitleInds = [1 2];

elseif colInCount == 3
  sigIndTmp = nan(unitCount, 7);
  
  sigIndTmp(:,1) = sigInd(:,1) & ~sigInd(:,2) & ~sigInd(:,3); % A Count
  sigIndTmp(:,2) = sigInd(:,1) & sigInd(:,2) & ~sigInd(:,3);  % AB Count
  sigIndTmp(:,3) = ~sigInd(:,1) & sigInd(:,2) & ~sigInd(:,3); % B Count
  sigIndTmp(:,4) = ~sigInd(:,1) & sigInd(:,2) & sigInd(:,3);  % BC Count
  sigIndTmp(:,5) = ~sigInd(:,1) & ~sigInd(:,2) & sigInd(:,3); % C Count
  sigIndTmp(:,6) = sigInd(:,1) & ~sigInd(:,2) & sigInd(:,3);  % AC Count
  sigIndTmp(:,7) = sigInd(:,1) & sigInd(:,2) & sigInd(:,3);   % ABC Count
  
  % Chance rates.
  chanceAB = (sum(sigInd(:,1))/unitCount) * (sum(sigInd(:,2))/unitCount) * unitCount;
  chanceBC = (sum(sigInd(:,2))/unitCount) * (sum(sigInd(:,3))/unitCount) * unitCount;
  chanceAC = (sum(sigInd(:,1))/unitCount) * (sum(sigInd(:,3))/unitCount) * unitCount;
  
  chanceCounts = round([chanceAB, chanceBC, chanceAC]);
  chanceTitleInds = [1 2; 2 3; 1 3];
else
  error('Input must have 2 or 3 columns')
end

sigCountVenn = sum(sigIndTmp);

% Plot
% Venn Diagram Code
vennX(sigCountVenn, 0.01);
figH = gcf;
figH.Children.FontSize = 13;
figH.Units = 'normalized';
if colInCount == 3
  figH.Position = [0.4 0.2 .46, .62];
elseif colInCount == 2
  figH.Position = [0.4    0.2    0.46    0.4691];
end
figH.NumberTitle = 'off';

figH.Name = figTitle;
title(figTitle)

tmpAx = gca;
numberHandles = tmpAx.Children(1:end-1);
for num_i = 1:length(numberHandles)
  if strcmp(numberHandles(num_i).String, '0')
    delete(numberHandles(num_i));
  else
    numberHandles(num_i).FontSize = 15;
  end
end

% Add Labels
text(0.06, 0.95, colNames{1}, 'units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');
text(0.856, 0.95, colNames{2}, 'units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');
if length(sigCountVenn) == 7
  text(0.63, 0.045, colNames{3}, 'units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');
end

labelCount = 0;
for ii = 1:size(chanceTitleInds,1)
  text2Insert = sprintf('Chance, %s + %s, %d', colNames{chanceTitleInds(ii, 1)}, colNames{chanceTitleInds(ii, 2)}, chanceCounts(ii));
  tmpX = text(0.016, 0.12 - (0.04 * labelCount), text2Insert, 'units', 'normalized', 'FontSize', 10);
  labelCount = labelCount + 1;
end
function vennXExpanded(sigInd, figTitle, colNames)
% A wrapper for vennX, which makes it more straight forward to turn an
% array of logical values into the appropriate cross sections. 

scrambleText = false;
[unitCount, colInCount] = size(sigInd);
numberSize = 24;
labelSize = 20;
chanceLabelSize = 12;

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
  chanceA = round((sum(sigInd(:,1))/unitCount) * 100);
  chanceB = round((sum(sigInd(:,2))/unitCount) * 100);
  chanceC = round((sum(sigInd(:,3))/unitCount) * 100);
  
  chanceAB = (sum(sigInd(:,1))/unitCount) * (sum(sigInd(:,2))/unitCount) * unitCount;
  chanceBC = (sum(sigInd(:,2))/unitCount) * (sum(sigInd(:,3))/unitCount) * unitCount;
  chanceAC = (sum(sigInd(:,1))/unitCount) * (sum(sigInd(:,3))/unitCount) * unitCount;
  
  chanceCounts = round([chanceAB, chanceBC, chanceAC]);
  chanceTitleInds = [1 2; 2 3; 1 3];
else
  error('Input must have 2 or 3 columns')
end

% Calculate error bars
if scrambleText
  countPerEvent = sum(sigInd);
  eventOverlapMeanSD = nan(size(chanceTitleInds));
  for chance_i = 1:size(chanceTitleInds, 1)
    scrambSet = nan(1000,1);
    for scramb_i = 1:1000
      % Randomly grab X of the total population
      scrambMat = zeros(unitCount,2);
      for uni_p = 1:2
        scrambMat(randperm(unitCount, countPerEvent(chanceTitleInds(chance_i, uni_p))), uni_p) = deal(1);
      end    

      scrambSet(scramb_i) = sum(sum(scrambMat,2) == 2);

    end
    eventOverlapMeanSD(chance_i, 1) = mean(scrambSet);
    eventOverlapMeanSD(chance_i, 2) = std(scrambSet);
  end

  % Once more for the 3 way
  if size(sigInd,2) == 3
    scrambSet = nan(1000,1);
    for scramb_i = 1:1000
      % Randomly grab X of the total population
      scrambMat = zeros(unitCount,3);
      for uni_p = 1:3
        scrambMat(randperm(unitCount, countPerEvent(uni_p)), uni_p) = deal(1);
      end

      scrambSet(scramb_i) = sum(sum(scrambMat,2) == 3);

    end
    eventOverlapMeanSD(4,1) = mean(scrambSet);
    eventOverlapMeanSD(4,2) = std(scrambSet);
  end

end

sigCountVenn = sum(sigIndTmp);

% Plot
% Venn Diagram Code
vennX(sigCountVenn, 0.01);
figH = gcf;
figH.Children.FontSize = 13;
figH.Units = 'normalized';
if colInCount == 3
  figH.Position = [0.4 0.2 .53, .72];
elseif colInCount == 2
  figH.Position = [0.4 0.2 .5  .5];
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
    numberHandles(num_i).FontSize = numberSize;
  end
end

if scrambleText
  plotTags = {'A', 'B', 'C'};
  colNames = strcat(colNames, '^', plotTags(1:length(colNames)));
end

% Add Labels
text(0.068, 0.95, colNames{1}, 'units', 'normalized', 'FontSize', labelSize, 'FontWeight', 'bold');
text(0.65, 0.95, colNames{2}, 'units', 'normalized', 'FontSize', labelSize, 'FontWeight', 'bold');
if length(sigCountVenn) == 7
  text(0.53, 0.045, colNames{3}, 'units', 'normalized', 'FontSize', labelSize, 'FontWeight', 'bold');
end

if scrambleText
  labelCount = 0;
  for ii = 1:size(chanceTitleInds,1)
    text2Insert = sprintf('%s + %s, %s +/- %s', plotTags{chanceTitleInds(ii, 1)}, plotTags{chanceTitleInds(ii, 2)}, ...
      num2str(eventOverlapMeanSD(ii,1),4), num2str(eventOverlapMeanSD(ii,2),2));
    tmpX = text(0.016, 0.12 - (0.04 * labelCount), text2Insert, 'units', 'normalized', 'FontSize', chanceLabelSize);
    labelCount = labelCount + 1;
  end

  if size(eventOverlapMeanSD,1) == 4
    text2Insert = sprintf('%s + %s + %s, %s +/- %s', plotTags{1}, plotTags{2}, plotTags{3}, ...
      num2str(eventOverlapMeanSD(4,1),4), num2str(eventOverlapMeanSD(4,2),2));
    tmpX = text(0.016, 0.12 - (0.04 * labelCount), text2Insert, 'units', 'normalized', 'FontSize', chanceLabelSize);
  end
end

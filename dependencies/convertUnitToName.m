function stringName = convertUnitToName(unitNumber, totalUnitCount, outputType)
% Assumes phyzzy convention, returns string with appropriate Label
% outputType - 1 = Named labels, which can be used in figures, 2 = US, U or M.

stringName = cell(length(unitNumber), 1);

for ii = 1:length(unitNumber)
  
  switch outputType
    case 1
      if unitNumber(ii) == totalUnitCount
        stringName{ii} = 'MUA';
      elseif unitNumber(ii) == 1
        stringName{ii} = 'US';
      else
        stringName{ii} = sprintf('U%d', unitNumber(ii)'-1);
      end
    case 2
      if unitNumber(ii) == totalUnitCount
        stringName{ii} = 'M';
      elseif unitNumber(ii) == 1
        stringName{ii} = 'US';
      else
        stringName{ii} = 'U';
      end
  end
    
end

% for a bit of backwards compatibility/ease
if length(stringName) == 1
  stringName = stringName{1};
end

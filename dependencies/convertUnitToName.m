function stringName = convertUnitToName(unitNumber, totalUnitCount, outputType)
% Assumes phyzzy convention, returns string with appropriate Label
% outputType - 1 = Named labels, which can be used in figures, 2 = US, U or M.

switch outputType
  case 1
    if unitNumber == totalUnitCount
      stringName = 'MUA';
    elseif unitNumber == 1
      stringName = 'Unsorted';
    else
      stringName = sprintf('U%d', unitNumber'-1);
    end
  case 2
    if unitNumber == totalUnitCount
      stringName = 'M';
    elseif unitNumber == 1
      stringName = 'US';
    else
      stringName = 'U';
    end
end
    
end
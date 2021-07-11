function selTable = initializeSelTable(figStruct, dateSubject, runNum, gridHole, recDepth)

channelUnitNames = figStruct.channelUnitNames;
channelNames = figStruct.channelNames;

chanCount = length(channelUnitNames);
gridHoleVec = cell(chanCount,1);
gridHoles = strsplit(gridHole, ';');
gridHoles = strtrim(gridHoles);

if contains(dateSubject, 'Mo')
  % I used 32 channel probes in Mo.
  switch chanCount
    case 32
      [gridHoleVec(:)] = deal(gridHoles(1));
    case 33
      [gridHoleVec{1:32}] = deal(gridHoles{1});
      [gridHoleVec{33}] = deal(gridHoles{2});
    case 64
      [gridHoleVec{1:32}] = deal(gridHoles{1});
      [gridHoleVec{33:64}] = deal(gridHoles{2});
    otherwise
      for ii = 1:length(gridHoles)
        gridHoleVec(ii) = gridHoles(ii);
      end
  end
  
else
  % I use 16 channel probes in Sam.
  switch chanCount
    case 16
      [gridHoleVec(:)] = deal(gridHoles(1));
    case 32
      [gridHoleVec(1:16)] = deal(gridHoles(1));
      [gridHoleVec(17:32)] = deal(gridHoles(2));
    otherwise
      for ii = 1:length(gridHoles)
        gridHoleVec(ii) = gridHoles(ii);
      end
  end
end
gridHoleVecTmp = string(gridHoleVec);

% Identify Channels, initialize cell table.
unitPerChan = cellfun('length', channelUnitNames);

channelVec = arrayfun(@(chanInd) repmat(channelNames(chanInd), [unitPerChan(chanInd), 1]), 1:length(channelNames), 'UniformOutput', false);
channelVec = string(vertcat(channelVec{:}));

gridHoleVec = arrayfun(@(chanInd) repmat(gridHoleVecTmp(chanInd), [unitPerChan(chanInd), 1]), 1:length(gridHoleVecTmp), 'UniformOutput', false);
gridHoleVec = string(vertcat(gridHoleVec{:}));

unitTypeVec = string([channelUnitNames{:}]');
dateSubjVec = repmat(string(dateSubject), [length(channelVec), 1]);
runNumVec = repmat(string(runNum), [length(channelVec), 1]);

recDepthVec = repmat(string(recDepth), [length(channelVec), 1]);

selTable = table(dateSubjVec, runNumVec, gridHoleVec, channelVec, unitTypeVec, recDepthVec);
selTable.Properties.VariableNames = extractBefore(selTable.Properties.VariableNames, 'Vec');

end
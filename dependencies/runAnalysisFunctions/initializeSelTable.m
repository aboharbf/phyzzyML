function selTable = initializeSelTable(figStruct, dateSubject, runNum, gridHole, recDepth, spikesByChannel)

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

meanWaveFormVec = cell(size(recDepthVec));

figStruct.closeFig = 1;
figStruct.exportFig = 1;
figStruct.saveFig = 0;

for chan_i = 1:length(spikesByChannel)
  % Identify unique units
  unitsOnChan = unique(spikesByChannel(chan_i).units)';
  unitsOnChan = unitsOnChan(unitsOnChan > 0);
  
  % Find their mean wave forms
  for unit_i = unitsOnChan
    % Find the mean waveform
    unitIndex = spikesByChannel(chan_i).units == unit_i;
    waveFormUnit = spikesByChannel(chan_i).waveforms(unitIndex,:);
    randSamp = randperm(size(waveFormUnit,1), min(size(waveFormUnit, 1), 1000));
    meanWaveFormTrace = mean(waveFormUnit);
    
    % Plotting - dont make it if its already there
    waveFormTitle = sprintf('%s U%d waveForms', figStruct.channelNames{chan_i}, unit_i);
    waveFormFigPath = fullfile(figStruct.figDir, strrep(waveFormTitle, ' ', '_'));
    
    if ~exist(strcat(waveFormFigPath, '.png'), 'file')
      figure('Name', waveFormTitle, 'NumberTitle', 'off');
      plot(waveFormUnit(randSamp,:)'); hold on; plot(meanWaveFormTrace, 'linewidth', 3, 'color', 'k');
      title(waveFormTitle);
      xlabel('Sample');
      ylabel('Amplitude (mV)');
      xlim([1 length(meanWaveFormTrace)])
      
      % Save
      saveFigure(figStruct.figDir, strrep(waveFormTitle, ' ', '_'), [], figStruct, figStruct.figTag);
      
    end
    
    % Slot in
    unitTableIndex = strcmp(unitTypeVec, sprintf('U%d', unit_i)) & strcmp(channelVec, figStruct.channelNames{chan_i});
    meanWaveFormVec{unitTableIndex} = meanWaveFormTrace;

  end
  
end

selTable = table(dateSubjVec, runNumVec, gridHoleVec, channelVec, unitTypeVec, recDepthVec, meanWaveFormVec);
selTable.Properties.VariableNames = extractBefore(selTable.Properties.VariableNames, 'Vec');

end
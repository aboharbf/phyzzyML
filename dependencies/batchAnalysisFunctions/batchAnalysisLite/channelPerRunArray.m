function spikePathBank = channelPerRunArray(spikePathBank);

% Extract the run list
runList = extractAfter(spikePathBank.Row, 'S');
chanPerRun = spikePathBank.chanCount;
selTablePerRun = spikePathBank.selTable;

channelsInRunArray = cell(size(selTablePerRun));
% For each of the entries
for ii = 1:length(runList)
  % Open the selTable and identify unique channel count
  channelsInRun = unique(double(extractAfter(selTablePerRun{ii}.channel, 'Ch')));
  
  % Make sure it matches the spikePathBank count before storing
  if length(channelsInRun) == chanPerRun(ii)
    channelsInRunArray{ii} = channelsInRun;
  else
    error('mismatch in selTable w/ spikePathBank');
  end
end

spikePathBank.channelNumbers = channelsInRunArray;


end
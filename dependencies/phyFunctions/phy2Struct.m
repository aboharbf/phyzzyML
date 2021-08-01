function [tmpSpikes, unitLabels]  = phy2Struct(dateSubjNum, params)
% a function which converts the typical set of .tsv files produced by Phy
% GUI into a .mat file.

addpath(genpath(params.phyParams.phyPath));
addpath(genpath(params.phyParams.spikesDir));
plotUnits = false;

% Find the directory
binFile = dir(fullfile(params.phyParams.ephysBinVolume, '**', [dateSubjNum, '.bin']));
binFile = fullfile(binFile.folder, binFile.name);
myKsDir = fileparts(binFile);
[~, B, ~] = fileparts(myKsDir);
sp = loadKSdir(myKsDir);

% Identify present .tsv files
clusInfoFiles = dir(fullfile(myKsDir, 'cluster_info.tsv'));
clusInfoFiles = string(fullfile(clusInfoFiles.folder, clusInfoFiles.name));
assert(~isempty(clusInfoFiles), 'cluster_info file missing in %s, make sure sorted w/ Phy', myKsDir)
clusInfoStruct= tdfread(clusInfoFiles);
wfs = params.phyParams.waveFormSize;

% Parameters for gwf.
gwfparams.dataDir = myKsDir;              % KiloSort/Phy output folder
gwfparams.fileName = strcat(B, '.bin');   % .dat file containing the raw
gwfparams.dataType = 'int16';             % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = sp.n_channels_dat;        % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = -(wfs/2)+1:(wfs/2);     % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 1e3;                     % Number of waveforms per unit to pull out
gwfparams.inputType = 'dat';              % string, 'dat' (default), 'ns5'

gwfparams.spikeTimes =    sp.st * sp.sample_rate;  % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = sp.clu; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

% Retrieve all the waveforms, do this by setting the retrieve count to be
% the maximum number of spikes in any cluster
[unitIDs, ~, ic] = unique(sp.clu);
a_counts = accumarray(ic,1);
numUnits = size(unitIDs,1);
gwfparams.nWf = max(a_counts);                     % Number of waveforms per unit to pull out

wf = getWaveForms(gwfparams);

% Process waveForm means a bit - find out the largest amplitude, and
waveFormMeans = wf.waveFormsMean;
mainChanPerUnit = zeros(length(wf.unitIDs), 1);
for unit_i = 1:size(waveFormMeans,1)
  
  % Identify channel the unit is 'actually' on.
  unitData = squeeze(waveFormMeans(unit_i,:,:));
  maxPeak = max(abs(unitData), [], 2);
  [~, mainChanPerUnit(unit_i)] = max(maxPeak);
  
  if plotUnits
    % Identify spikes belonging to this unit
    unitSpikeInd = sp.clu == wf.unitIDs(unit_i);
    figure()
    hold on
    plot(squeeze(wf.waveForms(unit_i, :, mainChanPerUnit(unit_i), :))');   % Plot all the traces for the unit
    plot(squeeze(waveFormMeans(unit_i, mainChanPerUnit(unit_i), :)), 'linewidth', 3, 'color', 'k')   % Plot the mean
    title(sprintf('Unit %d, Ch %d', unit_i, mainChanPerUnit(unit_i)));
  end
  
end

% Generate output vectors in rough state
[tmpSpikes.TimeStamp, tmpSpikes.Electrode, tmpSpikes.Unit, tmpSpikes.Waveform] = deal([]);
tmpSpikes.Unit = sp.clu;
tmpSpikes.TimeStamp = sp.st * sp.sample_rate;
for unit_i = 1:length(wf.unitIDs)
  % Identify units
  unitInd = sp.clu == wf.unitIDs(unit_i);
  % Pull waveforms, and stack
  waveFormsUnit = squeeze(wf.waveForms(unit_i, 1:sum(unitInd), mainChanPerUnit(unit_i), :));
  tmpSpikes.Waveform = [tmpSpikes.Waveform; waveFormsUnit];
  tmpSpikes.Electrode = [tmpSpikes.Electrode; repmat(mainChanPerUnit(unit_i), [sum(unitInd),1])];
end

% Fix 1 - Go through channel by channel and rename units according to their
% channel
for ii = 1:length()
  
end

% Outputs of interest
unitLabels = strtrim([string(clusInfoStruct.group), clusInfoStruct.KSLabel]);
channelPerUnit = mainChanPerUnit;

end
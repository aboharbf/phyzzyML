function  [spikeDataBank, stimPSTH] = noveltyAnalysis(spikeDataBank, stimPSTH, meanPSTHStruct, frameFiringStruct, params, figStruct)
% Function seeks to analyze whether most active PSTHes are enriched from
% novel runs.
% Inputs:
% - spikeDataBank
% - stimPSTH - the output of meanPSTH, which has indexing information in
% meanPSTHStruct.IndStructs.
disp('Starting Novelty Analysis...');

if ~exist(params.outputDir, 'dir')
  mkdir(params.outputDir)
end

% Plot 1 - RasterStack for stimuli 


% Figure 1 - iterate through stimuli, grabbing mean and max values across
% presentations in order, plot them on a single line. Mark dates where long
% recording dates were taken with vertical lines.


end

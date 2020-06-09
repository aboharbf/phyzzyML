%% Neural Decoding Toolbox Testing

% Add Phyzzy
addpath(genpath('C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML'));

paramFile = 'buildBatchAnalysisParamFileSocialVids';
paramPath = eval(paramFile);
tmp = load(paramPath);
addpath(genpath('C:\OneDrive\Lab\ESIN_Ephys_Files\Analysis\phyzzyML'))

NeuralDecodingTB([], tmp.NDTParams)
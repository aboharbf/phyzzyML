function varargout = spikePathLoad(spikePathBank, fields2Retrieve, params)
% A function which loads variables across different .mat files, and returns
% the output
% Inputs:
% - spikePathBank: has column 'analyzedDir', pointing to the analysis
% directory for each run
% - fields2Retrieve: a cell array with the name of the field(s) to be
% retrieved from each run.
% - params: a param file containing.

analyzedDirVec = 1:length(spikePathBank.analyzedDir);

varargout = cell(size(fields2Retrieve));
[varargout{:}] = deal(cell(length(analyzedDirVec),1));

for ii = 1:length(fields2Retrieve)
  variableFile = params.files(params.fileInds(strcmp(params.fileVars, fields2Retrieve{ii})));
  assert(~isempty(variableFile), "file associated w/ '%s' not found, check spelling or update ParamFile to specify variable's file", fields2Retrieve{ii});
  fileList = fullfile(spikePathBank.analyzedDir, variableFile);
  
  % Load those variables from each file
  for jj = analyzedDirVec
      tmp = load(fileList{jj}, fields2Retrieve{ii});
      assert(~isempty(fields(tmp)), 'file %s is missing variable %s', fileList{jj}, fields2Retrieve{ii});
      varargout{ii}{jj} = tmp.(fields2Retrieve{ii});
  end
  
end
  
end
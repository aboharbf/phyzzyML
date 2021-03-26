function varargout = spikePathLoad(spikePathBank, fields2Retrieve, params)
% A function which loads variables across different .mat files, and returns
% the output
% Inputs:
% - spikePathBank: has column 'analyzedDir', pointing to the analysis
% directory for each run
% - fields2Retrieve: a cell array with the name of the field(s) to be
% retrieved from each run.
% - params: a param file containing the files, variables, and indicies.

% Check that all the variables are in the list
varsNotRep = setdiff(fields2Retrieve, params.fileVars);
assert(isempty(varsNotRep), 'Variables set to retrieve are not all represented in the variable back, variables are spelled corrected in inputs and paramFile.')

analyzedDirVec = 1:length(spikePathBank.analyzedDir);

varPresent = false(length(analyzedDirVec), length(fields2Retrieve));

% Compile variables which each file needs to have extracted
filePerField2Retrieve = cellfun(@(x) params.fileInds(strcmp(params.fileVars, x)), fields2Retrieve);
files2Load = unique(filePerField2Retrieve);
fileNames2Load = params.files(files2Load)';
[varLoadStrings, variables2Load] = deal(cell(length(files2Load), 1));
for string_i = 1:length(varLoadStrings)
  tmp = fields2Retrieve(filePerField2Retrieve == files2Load(string_i));
  variables2Load{string_i} = string(strcat(tmp));
  varLoadStrings{string_i} = strjoin(["'" strjoin(variables2Load{string_i}, "', '") "'"], '');
end

% Initialize Outputs.
varargout = cell(size(fields2Retrieve));
[varargout{:}] = deal(cell(length(analyzedDirVec),1));

% Create string padding
strPads = repmat("'", [length(spikePathBank.analyzedDir), 1]);

for file_i = 1:length(fileNames2Load)
  fileList = fullfile(spikePathBank.analyzedDir, fileNames2Load{file_i});
  fileList = join([strPads fileList strPads], '');
  
  % Load those variables from each file
  for jj = analyzedDirVec
    eval(sprintf('tmp = load(%s, %s);', fileList{jj}, varLoadStrings{file_i}));
    if ~isempty(fields(tmp))
      % find what was loaded to avoid string processing further.
      loadedVars = fields(tmp);
      
      % sort them into their variable slots.
      [vars2Sort, varInd] = intersect(fields2Retrieve, loadedVars);
      
      % Quick Check for errors
      missingVars = setdiff(variables2Load{file_i}, vars2Sort);
      if ~isempty(missingVars)
        for miss_i = 1:length(missingVars)
          missInd = strcmp(fields2Retrieve, missingVars(miss_i));
          varPresent(jj, missInd) = true;
        end
      end
      
      % cycle through variables, storing them correctly
      for var_i = 1:length(vars2Sort)
        varargout{varInd(var_i)}{jj} = tmp.(vars2Sort{var_i});
      end
    else
      
      % Add these to missing var mat.
      missingVars = variables2Load{file_i};
      for miss_i = 1:length(missingVars)
        missInd = strcmp(fields2Retrieve, missingVars(miss_i));
        varPresent(jj, missInd) = true;
      end

    end
    
  end
end

if any(any(varPresent))
  disp(spikePathBank.analyzedDir(any(varPresent,2)));
  error('Not all variables in all runs, check files');
end
  
end
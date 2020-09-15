function output = saveSpikeDataBank(structData, parts, SaveOrLoad, outDir)
% Save or loads large structs into slices. used by preprocessBatchAnalysis
% and batchAnalysis for intermittent saving.
% - structData: The large struct
% - parts: the number of parts to divide it into.
% - SaveOrLoad: whether to 'save' or 'load' the struct.
switch SaveOrLoad
  case 'save'
    if ~exist(outDir, 'dir')
      mkdir(outDir);
    end
    
    if isempty(parts)
      previousSlices = dir(fullfile(outDir, 'spikeDataBank_*.mat'));
      parts = length(previousSlices);
    end
    
    %saves variable in N parts.
    fieldsToStore = fields(structData);
    fieldCount = length(fieldsToStore);
    fieldsPerFile = ceil(fieldCount/parts);
    
    starts = round(linspace(1, fieldCount, parts+1));
    ends = [starts(2:end-1)-1 starts(end)];
    starts = starts(1:end-1);
    
    output = cell(ceil(fieldCount/fieldsPerFile),1);
    
    for part_i = 1:parts
      % Remove the data not to be included in this slice
      fieldInd = zeros(fieldCount, 1);
      fieldInd(starts(part_i):ends(part_i)) = 1;
      spikeDataBankSlice = rmfield(structData, fieldsToStore(~fieldInd));
    
      % Generate the Path
      output{part_i} = fullfile(outDir, sprintf('spikeDataBank_%d.mat', part_i));

      % Save the file
      assert(isvarname('spikeDataBankSlice'), "spikeDataBankSlice doesn't exist")
      fprintf('Saving %s now... \n', output{part_i});
      save(output{part_i},'spikeDataBankSlice')
      clear spikeDataBankSlice
    end
    
  case 'load'
    fprintf('loading spikeDataBank...\n')
    spikeDataBankTmp = struct();
    filesToCombine = dir(fullfile(outDir, ['spikeDataBank_*.mat']));
    
    % Arrange correctly. In cases over 10, sorting is messed up.
    sliceNum = extractBetween([{filesToCombine.name}'], '_', '.');
    sliceNum = cellfun(@(x) str2num(x), sliceNum, 'UniformOutput', 0);
    sliceNum = [sliceNum{:}];
    [~, newInd] = sort(sliceNum);
    filesToCombine = filesToCombine(newInd);
    
    % Loop through files and concatonate them into output struct.
    for ii = 1:length(filesToCombine)
      tmp = load(fullfile(filesToCombine(ii).folder, filesToCombine(ii).name));
      fieldsToCombine = fields(tmp.spikeDataBankSlice);
      for field_ind = 1:length(fieldsToCombine)
        spikeDataBankTmp.(fieldsToCombine{field_ind}) = tmp.spikeDataBankSlice.(fieldsToCombine{field_ind});
      end
    end
    output = spikeDataBankTmp;
end
end
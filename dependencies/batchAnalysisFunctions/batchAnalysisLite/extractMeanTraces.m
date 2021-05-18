function [perNameTraces, perNameErrTraces, counts, nameVector] = extractMeanTraces(perRunVector, referenceVector, nameVector, allMeanSwitch, removeEmpty)
% Iterates through perRunVector, assuming a 
% structure, and looks for traces associated w/ each entry in the
% nameVector. Follows phzzy 'unit' convention (Unsorted, Units, MUA).
% Input:
% - perRunVector: contents {N}{chan}{unit}(stim * ms)
% - nameVector: a cell array of names represented in members of
% referenceVector, M in length.
% - referenceVector: an N*1 cell array, containing elements of nameVector.
% - allMeanSwitch: True/false to concatonate all mean at the bottom.

% Output:
% perNameTraces: mean traces collected for every member of nameVector (M*3
% or M+1 * 3 depending on allMeanSwitch.
% counts: the number of traces used to generate the means in perNameTraces.

% Extract unique entries
if isempty(nameVector)
  nameVector = unique(vertcat(referenceVector{:}));
end

if size(nameVector,1) == 1
  nameVector = nameVector';
end

maxTraceCount = 2000;
% Initialize an array as being far larger than it needs to be.
perNameTraces = cell(length(nameVector), 3);
[perNameTraces{:}] = deal(nan(maxTraceCount, size(perRunVector{1}{1}{1}, 2)));

for stim_i = 1:length(nameVector) 
  for run_i = 1:length(referenceVector)
    % Find out if the run has the stimuli
    stimRunInd = strcmp(nameVector{stim_i}, referenceVector{run_i});
    
    % If it does, cycle through the activity and collect the right vector.
    if any(stimRunInd)
      for chan_i = 1:length(perRunVector{run_i})
        unitCount = length(perRunVector{run_i}{chan_i});
        
        % 0 spike channels have only 1 unit, skip it
        if unitCount == 1
          continue
        end
        
        for unit_i = 1:unitCount
          
          % If there is activity to sort
          if ~isempty(perRunVector{run_i}{chan_i}{unit_i})
          
            % Classify unit as Unsorted, MUA, or Unit.
            if unit_i == 1
              group_i = 1;
            elseif unit_i == unitCount
              group_i = 3;
            else
              group_i = 2;
            end
          
          % Put the trace in the right cell
          
          % Find the next empty space in the cell, and add them there
          traces2Add = perRunVector{run_i}{chan_i}{unit_i}(stimRunInd, :);
          indices2Write = find(isnan(perNameTraces{stim_i, group_i}(:,1)), size(traces2Add, 1));
          perNameTraces{stim_i, group_i}(indices2Write,:) = traces2Add;
          
          end
        end
      end      
    end
  end
  
  % Get rid of excess NaN Padding here
  for group_i = 1:3
    perNameTraces{stim_i, group_i} = perNameTraces{stim_i, group_i}(~isnan(perNameTraces{stim_i, group_i}(:,1)),:);
  end
  
end

% Get rid of excess cells
if allMeanSwitch
  grandMean = cell(1,3);
  for group_i = 1:3
    grandMean{1, group_i} = vertcat(perNameTraces{:, group_i});
  end
  
  perNameTraces = [perNameTraces; grandMean];
  if any(contains(nameVector, '.avi'))
    nameVector = extractBefore(nameVector, '.avi'); 
  end
  nameVector = [nameVector; 'Grand Mean'];
end

% Process the collected data into means and counts
counts = cellfun(@(traces) size(traces, 1), perNameTraces);
perNameErrTraces = cellfun(@(traces) nanstd(traces)/sqrt(size(traces,1)), perNameTraces, 'UniformOutput', false);
perNameTraces = cellfun(@(traces) nanmean(traces, 1), perNameTraces, 'UniformOutput', false);

% check for empty cells, as this messes w/ later plotting
if any(cellfun('isempty', perNameTraces)) & removeEmpty
  keepInd = ~cellfun('isempty', perNameTraces(:,1));
  perNameTraces = perNameTraces(keepInd, :);
  perNameErrTraces = perNameErrTraces(keepInd, :);
  counts = counts(keepInd, :);
  nameVector = nameVector(keepInd, :);
end

% Remove underscores from titles
nameVector = strrep(nameVector, '_', '');

end
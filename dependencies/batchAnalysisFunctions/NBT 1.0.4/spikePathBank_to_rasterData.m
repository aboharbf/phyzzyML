function  spikePathBank_to_rasterData(spikePathBank, rasterDataPath, params)
% Turns the spikeDataBank struct of processBatchAnalysis into the format
% needed to process by Neural Decoding Toolbox.

% Creates file with the following 3 fields
% raster_data - 1 * unit cell array, each cell containing trial * bins.
% raster_labels - struct with each field contain 1*trial labels. Labels
% like 'stimulus' or 'position' go here. These are currently generated
% based on the stimParamFile used by phyzzy combined with plotIndex. The
% labels and their fieldNames are defined in params.plotIndParams.
% raster_site_info - file*1 labels/numbers - session_ID, channel, unit,
% site.

%runData.spikesByEventBinned{stim}{chan}{unit} --> needs to be turned to
%rasterDataUnit

if ~exist(rasterDataPath, 'dir')
  mkdir(rasterDataPath);
end

% Load the eventData for attaching significance activity to units. This is
% generated during processRunBatch.
if exist(params.subEventBatchStructPath, 'file')
  tmp = load(params.subEventBatchStructPath);
  sigUnitGrid = tmp.subEventBatchStruct.sigUnitGrid;
  eventList = tmp.subEventBatchStruct.eventList;
end

% Load the spike data to be turned into raster data.
% [spikesByEventBinned, psthParams, selTable] = spikePathLoad(spikePathBank, {'spikesByEventBinned', 'psthParams', 'selTable'}, params.spikePathLoadParams);
% selTable = vertcat(selTable{:});
% runList = spikePathBank.Properties.RowNames;

[spikesByEventBinned, psthParams] = spikePathLoad(spikePathBank, {'spikesByEventBinned', 'psthParams'}, params.spikePathLoadParams);
selTable = vertcat(spikePathBank.selTable{:});
runList = spikePathBank.Properties.RowNames;

% epochSelInd = contains(selTable.Properties.VariableNames', 'BaseV') & contains(selTable.Properties.VariableNames', 'PVal');
% eventSelInd = contains(selTable.Properties.VariableNames', 'Sel_');
% removeInd = strcmp(selTable.Properties.VariableNames', 'saccSel_socVNonSoc_diff');
% var2AddInd = epochSelInd | eventSelInd & ~removeInd;
% 
% vars2Add = selTable.Properties.VariableNames(var2AddInd);
% vars2Add = vars2Add(~contains(vars2Add, 'broadCategoriesSel'));
% vars2Add = vars2Add(~contains(vars2Add, 'saccSel_broadCategories_diff'));

% Take sub selection here. give them shorter names for defining.
vars2AddOrig2 = selTable.Properties.VariableNames(contains(selTable.Properties.VariableNames', 'selInd'));
vars2AddOrig = [
  {'baseV_Fix_selInd'                          }
  {'baseV_stimEarly_selInd'                    }
  {'baseV_stimLate_selInd'                     }
  {'baseV_reward_selInd'                       }
  {'subSel_headTurn_all_selInd'                }
  {'subSel_rewardCombo_selInd'                 }
  {'epochSel_socVNonSoc_stimEarly_selInd'      }
  {'epochSel_socVNonSoc_stimLate_selInd'       }
  {'epochSel_socVNonSoc_reward_selInd'         }
  {'epochSel_socVNonSoc_any_selInd'            }
  {'epochSel_broadCategories_stimEarly_selInd' }
  {'epochSel_broadCategories_stimLate_selInd'  }
  {'epochSel_broadCategories_reward_selInd'    }
  {'epochSel_broadCategories_any_selInd'       }
  {'slidingWin_broadCatTest_selInd'            }
  {'slidingWin_categoriesTest_selInd'          }];

vars2Add = strrep(vars2AddOrig, 'epochSel_socVNonSoc', 'sVns');
vars2Add = strrep(vars2Add, 'epochSel_broadCategories', 'bCat');
vars2Add = strrep(vars2Add, 'slidingWin', 'SW');

% vars2Added = [    {'baseV_Fix_selInd'          }
%     {'baseV_stimEarly_selInd'    }
%     {'baseV_stimLate_selInd'     }
%     {'baseV_reward_selInd'       }
%     {'subSel_headTurn_all_selInd'}
%     {'subSel_rewardCombo_selInd' }
%     {'sVns_stimEarly_selInd'     }
%     {'sVns_stimLate_selInd'      }
%     {'sVns_reward_selInd'        }
%     {'sVns_any_selInd'           }
%     {'bCat_stimEarly_selInd'     }
%     {'bCat_stimLate_selInd'      }
%     {'bCat_reward_selInd'        }
%     {'bCat_any_selInd'           }
%     {'SW_broadCatTest_selInd'    }
%     {'SW_categoriesTest_selInd'  }];
  
countPerVar = zeros(length(vars2Add),1);

for run_i = 1:length(runList)
  
  % print a message the the data is being binned (and add a dot for each file that has been binned
  curr_bin_string = ['Processing Runs: ' num2str(run_i) ' of ' num2str(length(runList))];
  if run_i == 1
    disp(curr_bin_string);
  else
    fprintf([repmat(8,1,bin_str_len) curr_bin_string]);
  end
  bin_str_len = length(curr_bin_string);
  
  % Collect data for this particular run
  runData = struct();
%   runData.start = -psthParams{run_i}.psthPre + psthParams{run_i}.ITI + params.fixShorten;
  runData.start = -psthParams{run_i}.psthPre - params.fixShorten;
  runData.end = psthParams{run_i}.psthImDur + psthParams{run_i}.psthPost;
  padSize = psthParams{run_i}.movingWin(1)/2;
  
  % Full Length = psthPre + psthImDur + psthPost + psthParams{run_i}.movingWin(1)
  
  runData.spikesByEventBinned = spikesByEventBinned{run_i};
  runData.eventIDs = spikePathBank.stimuli{run_i};
  
  % Indicies for activity, accounting for padding, ITI
  vecLength = size(runData.spikesByEventBinned{1}{1}{1},2);
  theoreticalLength = psthParams{run_i}.psthPre + psthParams{run_i}.psthImDur + psthParams{run_i}.psthPost + psthParams{run_i}.movingWin(1) + 1;
  assert(vecLength == theoreticalLength, 'Something is off w/ vector lengths');

  % Determine the slice of activity to be used.
%   startTime = psthParams{run_i}.psthPre - params.fixShorten + padSize;
  startTime = padSize;
  endTime = vecLength - padSize;
  dataLength = endTime - startTime + 1;
  
  % find out trials per stimuli to generate useful indicies
  trialsPerStim = cellfun(@(x) size(x{1}{1},1), runData.spikesByEventBinned);
  trialsPerStimStart = cumsum([1; trialsPerStim(1:end-1)]);
  trialsPerStimEnd = cumsum(trialsPerStim);
  trialsPerUnit = sum(trialsPerStim);
  stimCount = 1:length(trialsPerStim);
  
  % Generate the raster file fields which are the same for every unit.
  % raster_labels - labels that exist per trial. generate them per stimuli,
  % then repmat them.
  
  % Generate a vector for stimuli for each trial.
  stimuliVecTmp = cellfun(@(x) convertCharsToStrings(x), runData.eventIDs);
  stimuliVecTmp = arrayfun(@(x) repmat(stimuliVecTmp(x), [trialsPerStim(x), 1]), stimCount, 'UniformOutput', 0)';
  stimuliVec = vertcat(stimuliVecTmp{:});
  
  % Generate a vector for raster labels, with 0/1 for binary choices, or an
  % index for a larger set.
  params.plotIndParams.plotLabels = params.rasterParams.plotIndParams.plotLabels;
  rasterLabelTmp = plotIndex(runData.eventIDs, params.plotIndParams);
  rasterLabelTmp = arrayfun(@(x) repmat(rasterLabelTmp(x, :), [trialsPerStim(x), 1]), stimCount, 'UniformOutput', 0)';
  rasterLabel = vertcat(rasterLabelTmp{:});
    
  % Assign to raster labels.
  raster_labels = struct();
  raster_labels.stimuli = stimuliVec';
  %   raster_labels.stimPresCount = stimPresCountPerTrial';
  
  % iterate through rasterLabels, creating a string array which grants the
  % label of 'target' or 'non-target' to each trial.
  for label_i = 1:length(params.rasterParams.rasterLabels)
    % Convert the indices into strings
    label = params.rasterParams.plotIndParams.plotLabels{label_i};
    entry = cell(size(rasterLabel, 1), 1);
    [entry{:}] = deal('None');
    
    if iscell(label)
      % Use entries as indices to the provided labels
      entryInd = rasterLabel(:,label_i);
      for lab_i = 1:length(label)
        entry(entryInd == lab_i) = label(lab_i);
      end
    else
      % use single label, add 'Non' to make 2nd.
      entryInd = logical(rasterLabel(:,label_i));
      entry(entryInd) = {label};
      entry(~entryInd) = {['Non-' label]};
    end
    
    % Store them
    raster_labels.(params.rasterParams.rasterLabels{label_i}) = entry;
  end
  
  % raster_site_info
  raster_site_info.session_ID = runList{run_i};
  raster_site_info.alignment_event_time = abs(runData.start);
  
  % Construct raster_data and save file per unit.
  selTableRows = selTable(strcmp(selTable.dateSubj, runList{run_i}(2:end-3)) & strcmp(selTable.runNum, runList{run_i}(end-2:end)), :);
  channels = strcat('Ch', string(sort(double(extractAfter(unique(selTableRows.channel), 'Ch')))));
  
  for chan_i = 1:length(runData.spikesByEventBinned{1})
    
    % Units were not labeled as expected. Not every channel has unsorted
    % activity.
    chLabel = channels(chan_i);
    selTableChRows = selTableRows(strcmp(selTableRows.channel, chLabel), :);
    ULabelList = selTableChRows.unitType;
    
    for unit_i = 1:length(runData.spikesByEventBinned{1}{chan_i})
      
      % Actual ULabel
      ULabel = ULabelList{unit_i};
      
      % Cycle through the stim, store into raster_data.
      raster_data = zeros(trialsPerUnit, dataLength);
      for stim_i = 1:length(runData.spikesByEventBinned)
        raster_data(trialsPerStimStart(stim_i):trialsPerStimEnd(stim_i),:) = runData.spikesByEventBinned{stim_i}{chan_i}{unit_i}(:, startTime:endTime);
      end
      
      % Generate the raster_site_info
      raster_site_info.UnitType = convertUnitToName(unit_i, length(runData.spikesByEventBinned{1}{chan_i}), 2);
      
      rasterFileName = strjoin([runList(run_i), chLabel, ULabel], '_');
      
      % use selTable to add selectivities for individual events.
      gridRow = selTableChRows(unit_i, :);
      for event_i = 1:length(vars2Add)
        raster_site_info.(vars2Add{event_i}) = gridRow.(vars2AddOrig{event_i});
        countPerVar(event_i) = countPerVar(event_i) + sum(raster_site_info.(vars2Add{event_i}));
      end
      
      % If desired, remove data and accompanying labels if the labels are
      % empty (Specifically, in my FamilirFace paradigm, Empty cage
      % stimuli should have no labels).
      if unit_i == 1 && chan_i == 1
        removeInd = false(size(raster_data(:, 1)));
        for label_i = 1:length(params.rasterParams.rasterLabels)
          emptyEntries = cellfun('isempty', raster_labels.(params.rasterParams.rasterLabels{label_i}));
          removeInd = emptyEntries | removeInd;
        end
        
        % Remove elements in each and cut the data
        if any(removeInd)
          raster_labels.stimuli = raster_labels.stimuli(~removeInd);
          for label_i = 1:length(params.rasterParams.rasterLabels)
            raster_labels.(params.rasterParams.rasterLabels{label_i}) = raster_labels.(params.rasterParams.rasterLabels{label_i})(~removeInd);
          end
        end
        
      end
      raster_data = raster_data(~removeInd,:);

      % Save file in output folder
      save(fullfile(rasterDataPath, rasterFileName), 'raster_data', 'raster_labels', 'raster_site_info')
      
    end
  end
  
end

end

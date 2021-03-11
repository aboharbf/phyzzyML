function [spikePathBank, params] = pullPSTH(spikePathBank, params)
% Function which combines stimulus presentations across all runs in the spikeDataBank.
% Inputs include spikeDataBank and list of parameters.
disp('Generating Normalized/Thresholded version of PSTHes, if desired...');

psthParams = params.meanPSTHParams;
processedVarFileInd = length(params.spikePathLoadParams.files) + 1;

if psthParams.normalize || psthParams.rateThreshold
  
  % extract the variables
  [psthByImage, psthErrByImage, psthByCategory, psthErrByCategory, psthParamsPerRun] = spikePathLoad(spikePathBank, {'psthByImage', 'psthErrByImage', 'psthByCategory', 'psthErrByCategory', 'psthParams'}, psthParams.spikePathLoadParams);
  runList = spikePathBank.Properties.RowNames;
  procFileName = params.spikePathLoadParams.batchAnalysisOutputName;
  procFileList = fullfile(spikePathBank.analyzedDir, procFileName);
  
  % For every variable above, normalize/threshold as needed, and save it to a
  % file in the local analysis folder.
  tic
  for run_i = 1:length(runList)
    
    % Initialize structures to save, processed versions.
    [psthByImageProc, psthErrByImageProc] = deal(initNestedCellArray(psthByImage{run_i}, 'NaN'));
    [psthByCategoryProc, psthErrByCategoryProc] = deal(initNestedCellArray(psthByCategory{run_i}, 'NaN'));
    
    for chan_i = 1:length(psthByImage{run_i})
      for unit_i = 1:length(psthByImage{run_i}{chan_i})
        % Retrieve correct PSTH from run
        unitActivity = psthByImage{run_i}{chan_i}{unit_i};
        unitErr = psthErrByImage{run_i}{chan_i}{unit_i};
        unitCategoryAct = psthByCategory{run_i}{chan_i}{unit_i};
        unitCategoryErr = psthErrByCategory{run_i}{chan_i}{unit_i};
        
        % Do not include activity which does not meet 1 Hz Threshold.
        if psthParams.rateThreshold ~= 0
          trialTime = size(unitActivity,2)/1000;
          averageTrialSpikes = sum(mean(unitActivity));
          minHz = psthParams.rateThreshold * trialTime;
          if averageTrialSpikes < minHz
            [unitActivity, unitErr] = deal(nan(size(unitActivity)));
            [unitCategoryAct, unitCategoryErr] = deal(nan(size(unitCategoryAct)));
          end
        end
        
        if psthParams.normalize == 1
          % If desired, Z score PSTHs here based on fixation period activity.
          runPsthParams = psthParamsPerRun{run_i};
          tmp = unitActivity(:,1:abs(runPsthParams.ITI));
          tmp = reshape(tmp, [size(tmp,1) * size(tmp,2), 1]);
          fixMean = mean(tmp); %Find activity during fixation across all stim.
          fixSD = std(tmp);
          if  fixMean ~= 0 && fixSD ~= 0
            unitActivity = ((unitActivity - fixMean)/fixSD);
            unitErr = ((unitErr - fixMean)/fixSD);
          end
          
          % Not sure if it is possible for this to yield different number,
          % but will do seperately just in case.
          tmp = unitCategoryAct(:,1:abs(runPsthParams.psthPre));
          tmp = reshape(tmp, [size(tmp,1) * size(tmp,2), 1]);
          fixMean = mean(tmp); %Find activity during fixation across all stim.
          fixSD = std(tmp);
          if  fixMean ~= 0 && fixSD ~= 0
            unitCategoryAct = ((unitCategoryAct - fixMean)/fixSD);
            unitCategoryErr = ((unitCategoryErr - fixMean)/fixSD);
          end
        end
        
        psthByImageProc{chan_i}{unit_i} = unitActivity;
        psthErrByImageProc{chan_i}{unit_i} = unitErr;
        psthByCategoryProc{chan_i}{unit_i} = unitCategoryAct;
        psthErrByCategoryProc{chan_i}{unit_i} = unitCategoryErr;
        
      end
    end
    
    save(procFileList{run_i}, 'psthByImageProc', 'psthErrByImageProc', 'psthByCategoryProc', 'psthErrByCategoryProc');
    
  end
  toc
  
  % Update the parameters passed in, so downstream functions will know how
  % to find these variables.
  procVarList = {'psthByImageProc', 'psthErrByImageProc', 'psthByCategoryProc', 'psthErrByCategoryProc'};
  params.spikePathLoadParams.files = [params.spikePathLoadParams.files, procFileName];
  params.spikePathLoadParams.fileVars = [params.spikePathLoadParams.fileVars, procVarList];
  params.spikePathLoadParams.fileInds = [params.spikePathLoadParams.fileInds, ones(1,length(procVarList))*processedVarFileInd];
  
end

if params.meanPSTHParams.removeRewardEpoch
  
  if params.meanPSTHParams.normalize
    varCheck = 'psthByImageProcNoReward';
  else
    varCheck = 'psthByImageNoReward';
  end
  
  % check to see if the reward-less PSTH exists
  if ~any(strcmp(params.spikePathLoadParams.fileVars, varCheck))
    % Find out which variables to process
    variables2Extract = {'psthByImage', 'psthErrByImage', 'psthByCategory', 'psthErrByCategory', 'psthParams'};
    if params.meanPSTHParams.normalize
      variables2Extract(1:4) = strcat(variables2Extract(1:4), 'Proc');
    else
      
    end
    variables2Save(1:4) = strcat(variables2Extract(1:4), 'NoReward');
    
    % Collect the vectors of activity from each analysis directory.
    [psthByImage, psthErrByImage, psthByCategory, psthErrByCategory, psthParamsPerRun] = spikePathLoad(spikePathBank, variables2Extract, params.spikePathLoadParams);
    analyzedDataBatchPaths = fullfile(spikePathBank.analyzedDir, params.spikePathLoadParams.batchAnalysisOutputName);
    
    for run_i = 1:length(psthByImage)
      % Identify the period of the activity to keep.
      nonRewardPeriod = 1:1+(psthParamsPerRun{run_i}.psthPre+psthParamsPerRun{run_i}.psthImDur);
      
      % Initialize new vectors for shorter vectors
      [A, B] = deal(initNestedCellArray(psthByImage{run_i}, 'NaN'));
      [C, D] = deal(initNestedCellArray(psthByCategory{run_i}, 'NaN'));
      
      % Cycle through and shorten the vectors.
      for chan_i = 1:length(psthByImage{run_i})
        for unit_i = 1:length(psthByImage{run_i}{chan_i})
          A{chan_i}{unit_i} = psthByImage{run_i}{chan_i}{unit_i}(:, nonRewardPeriod);
          B{chan_i}{unit_i} = psthErrByImage{run_i}{chan_i}{unit_i}(:, nonRewardPeriod);
          C{chan_i}{unit_i} = psthByCategory{run_i}{chan_i}{unit_i}(:, nonRewardPeriod);
          D{chan_i}{unit_i} = psthErrByCategory{run_i}{chan_i}{unit_i}(:, nonRewardPeriod);
        end
      end
      
      tmpVar = {'A', 'B', 'C', 'D'};
      
      % Save the new variables to the output file in each directory
      for var_i = 1:length(variables2Save)
        eval(sprintf('%s = %s;', variables2Save{var_i}, tmpVar{var_i}))
        if ~exist(analyzedDataBatchPaths{run_i})
          save(analyzedDataBatchPaths{run_i}, variables2Save{var_i})
        else
          save(analyzedDataBatchPaths{run_i}, variables2Save{var_i}, '-append')
        end
      end
    end
    
    % Update the param files, save to the original.
    if ~any(strcmp(params.spikePathLoadParams.fileVars, varCheck))
      params.spikePathLoadParams.files = [params.spikePathLoadParams.files, params.spikePathLoadParams.batchAnalysisOutputName];
    end
    params.spikePathLoadParams.fileVars = [params.spikePathLoadParams.fileVars, variables2Save];
    params.spikePathLoadParams.fileInds = [params.spikePathLoadParams.fileInds, ones(1,length(variables2Save))*processedVarFileInd];
    batchAnalysisParams = params;
    spikePathFile = batchAnalysisParams.spikePathLoadParams.batchAnalysisOutput;
    save(spikePathFile, 'batchAnalysisParams', '-append')
    
  end
end

end

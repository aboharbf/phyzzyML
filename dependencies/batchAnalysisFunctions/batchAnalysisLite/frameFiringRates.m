function [spikeDataBank, frameFiringStruct]  = frameFiringRates(spikeDataBank, params, figStruct)
% Generates means and distributions based on what is being looking at.
disp('Starting Frame Firing Rate Analysis...');

% Generate variables needed
runList = fields(spikeDataBank);

if ~exist(params.outputDir, 'dir')
  mkdir(params.outputDir)
end

% Step 1 - generate bin times and spike rates.
if ~isfield(spikeDataBank.(runList{end}), 'frameFiringRates')
  for run_ind = 1:length(runList)
    runStruct = spikeDataBank.(runList{run_ind});
    if length(unique(cellfun('length',runStruct.attendedObjData.frameStartInd))) == 1 % If all the videos are the same number of frames
      starts = runStruct.attendedObjData.frameStartInd{1}'+params.delay;
      ends = runStruct.attendedObjData.frameEndInd{1}'+params.delay;
      spikeDataBank.(runList{run_ind}).frameTimes = [starts,ends];
      spikeDataBank.(runList{run_ind}).frameFiringRates = cell(length(starts),1);
      frameFiringRatesTmp = cell(length(starts),1);
      for frame_ind = 1:length(starts)
        % spikeDataBank.run.framingFiringRates{epoch/bin}{channel}{unit}{stim}
        [frameFiringRatesTmp{frame_ind}, ~, ~] = spikeCounter(spikeDataBank.(runList{run_ind}).spikesByEvent, starts(frame_ind), ends(frame_ind));
      end
    else
      error('Different stimuli contain different number of frames. Implement different spikeCounting method.')
    end
    % Step 1.5 - Rearrange frameFiringRatesTmp to match attendedObjVect data {channel}{unit}{stimuli}(trials*frames)
    frameCount = length(frameFiringRatesTmp);
    trialCount = length(frameFiringRatesTmp{1}{1}{1}{1}.counts);
    frameFiringRates = initNestedCellArray(frameFiringRatesTmp{1},'zeros',[trialCount, frameCount]);
    for chan_ind = 1:length(frameFiringRatesTmp{1})
      for unit_ind = 1:length(frameFiringRatesTmp{1}{chan_ind})
        frameFiringRatesUnitTmp = initNestedCellArray(runStruct.attendedObjData.attendedObjVect, 'zeros', [size(runStruct.attendedObjData.attendedObjVect{1})],1);
        for frame_ind = 1:length(frameFiringRatesTmp)
          for stim_ind = 1:length(frameFiringRatesTmp{frame_ind}{chan_ind}{unit_ind})
            if frame_ind == 1
              frameFiringRatesUnitTmp{stim_ind} = zeros(size(runStruct.attendedObjData.attendedObjVect{stim_ind}));
            end
            if params.useRates
            frameFiringRatesUnitTmp{stim_ind}(:,frame_ind) = frameFiringRatesTmp{frame_ind}{chan_ind}{unit_ind}{stim_ind}.rates;
            else
              frameFiringRatesUnitTmp{stim_ind}(:,frame_ind) = frameFiringRatesTmp{frame_ind}{chan_ind}{unit_ind}{stim_ind}.counts;
            end
          end
        end
        frameFiringRates{chan_ind}{unit_ind} = frameFiringRatesUnitTmp;
      end
    end
    spikeDataBank.(runList{run_ind}).frameFiringRates = frameFiringRates;
    disp(run_ind)
  end
  saveSpikeDataBank(spikeDataBank, [], 'save',fileparts(params.outputDir));
else
  fprintf('Frame firing rates already calculated., continuing... \n');
end

% Step 2 - Go through each run, sum frame counts/rates associated with the
% specific objects being attended.
objList = spikeDataBank.(runList{1}).attendedObjData.objList;
objListBroad = {'Face','Body','Hand','GazeFollow','Object','Bkg'};

if params.broadLabels
  objListPlot = objListBroad;
else
  objListPlot = objList;
end

% Step 3 - Generate vectors with data from all runs, organized as follows
% from ObjFrameFiringRates{channel}{unit}{stim}{obj}{dataType} to
% objFrameFiringRatesTotal{stim}{obj}{groupingType}{dataType}
groupingType = {'Unsorted','Units','MUA'};
dataType = {'Raw','Mean','Run Ind'};

% Pooling across runs requires vector of all stim, extract the eventIDs field, generate a cell array of unique stimuli
allStimuliVec = struct2cell(structfun(@(x) x.eventIDs, spikeDataBank,'UniformOutput', 0));
[allStimuliVec, ~, ic] = unique(vertcat(allStimuliVec{:}));
allStimuliVecCounts = accumarray(ic,1);

objFrameFiringRatesTotal = cell([length(allStimuliVec),length(groupingType), length(objListPlot), length(dataType)]);
%objFrameFiringRatesTotal{stim_ind, group_ind, obj_ind, data_ind};

channelUnitMat = [];
for run_ind = 1:length(runList)
  frameFiringRates = spikeDataBank.(runList{run_ind}).frameFiringRates;
  attendedObjVect = spikeDataBank.(runList{run_ind}).attendedObjData.attendedObjVect;
  eventIDs = spikeDataBank.(runList{run_ind}).eventIDs;
  chanCount = 0;
  unitCount = 0;
  
  % Per Stimuli Plot
  for channel_ind = 1:length(frameFiringRates)
    chanCount = chanCount + 1;
    for unit_ind = 1:length(frameFiringRates{channel_ind})
      unitCount = unitCount + 1;
      if unit_ind == 1
        groupInd = 1;  %Unsorted
      elseif unit_ind == length(frameFiringRates{channel_ind}) % MUA
        groupInd = 3;  %MUA
      else
        groupInd = 2;  %Unit
      end
      for stim_ind = 1:length(attendedObjVect)
        stimStoreInd = strcmp(spikeDataBank.(runList{run_ind}).eventIDs{stim_ind}, allStimuliVec);
        stimAttendedObj = attendedObjVect{stim_ind};
        spikeStruct = frameFiringRates{channel_ind}{unit_ind}{stim_ind};
        %ObjFrameFiringRatesStim = initNestedCellArray([length(objListPlot), length(dataType)],'zeros');
        for obj_ind = 1:length(objList)
          objRates = spikeStruct(strcmp(stimAttendedObj, objList{obj_ind}));
          
          % Identify correct object group, when using broad labels.
          if params.broadLabels
            switch find(strcmp(objList,objList{obj_ind}))
              case {1,7};         objStoreInd = 1;
              case {2,8};         objStoreInd = 2;
              case {3,4,9,10};    objStoreInd = 3;
              case {5,11};        objStoreInd = 4;
              case {6,12};        objStoreInd = 5;
              case 13;            objStoreInd = 6;
            end
          else
            objStoreInd = obj_ind;
          end
          if ~isempty(objRates)
            %ObjFrameFiringRatesStim{objStoreInd}{1} = [ObjFrameFiringRatesStim{objStoreInd}{1}; objRates];
            objFrameFiringRatesTotal{stimStoreInd, groupInd, objStoreInd, 1} = [objFrameFiringRatesTotal{stimStoreInd, groupInd, objStoreInd, 1}; objRates];
            objFrameFiringRatesTotal{stimStoreInd, groupInd, objStoreInd, 2} = [objFrameFiringRatesTotal{stimStoreInd, groupInd, objStoreInd, 2}; mean(objRates)];
            objFrameFiringRatesTotal{stimStoreInd, groupInd, objStoreInd, 3} = [objFrameFiringRatesTotal{stimStoreInd, groupInd, objStoreInd, 3}; run_ind];
          end
        end
      end
    end
  end
  % Structure is ObjFrameFiringRates{stim}{obj}{channel}{unit}{dataType}
  channelUnitMat = [channelUnitMat; [run_ind, chanCount, unitCount]];
end
disp('Finished calculating objFrameFiringRatesTotal...')

% Generate the correct labels depending on whether broad labels are used
% for objects.
% Broad Label - Swaps labels of individual stimuli for broader catagories,
% as defined in stimParamFile and the parameters.
tmp = load(params.stimParamsFilename);
for event_i = 1:length(tmp.paramArray)
  totalEventIDs{event_i} = tmp.paramArray{event_i}{1}; %per the stimParamFile spec, this is the event ID
end
totalEventIDs = totalEventIDs';
[~, paramSortVec] = ismember(allStimuliVec, totalEventIDs);
paramArray = tmp.paramArray(paramSortVec);

%Iterate across stimuli and assign new labels.
allStimuliVecBroad = cell(length(allStimuliVec),1);
for label_ind = 1:length(allStimuliVec)
  paramStimSet = paramArray{label_ind};
  allStimuliVecBroad{label_ind} = paramStimSet{ismember(paramStimSet,params.broadLabelPool)};
end
[uniqueStimLabels, ~, stimIndex] = unique(allStimuliVecBroad);

% Find out how many of each of the newly labeled stimuli there are
% uniqueStimLabelCounts = cell2mat(arrayfun(@(x)length(find(stimIndex == x)), unique(stimIndex), 'Uniform', false));

% Step 5 - Generate figures
% Produces 72 images
for stim_ind = 1:length(uniqueStimLabels)
  for group_ind = 1:length(groupingType)
    figTitle = sprintf('%s, Per Object Rates, %s - %s', uniqueStimLabels{stim_ind}, dataType{2}, groupingType{group_ind});
    figFilePath = fullfile(params.outputDir, figTitle);
    if ~exist([figFilePath, '.fig'],'file')
      % Fig 1 - 'Chasing, Per Object Rates Means - MUA'
      stimMatInd = strcmp(uniqueStimLabels(stim_ind), allStimuliVecBroad);
      stimObjData = objFrameFiringRatesTotal(stimMatInd, group_ind, :,  2);
      plotArray = cell(length(objListPlot),2);
      removeInd = false(size(stimObjData,3),1);
      for obj_ind = 1:size(stimObjData,3)
        plotArray{obj_ind,1} = vertcat(stimObjData{:,:,obj_ind});
        removeInd(obj_ind) = isempty(plotArray{obj_ind,1});
        plotArray{obj_ind,2} = objListPlot{obj_ind};
      end
      plotArray(removeInd,:) = [];
      h = figure('Name', figTitle, 'NumberTitle', 'off');
      sgtitle(figTitle)
      hold on
      uniqueObjLabels = plotArray(:,2);
      [tmpValues, tmpLabels] = deal([]);
      for obj_ind = 1:length(uniqueObjLabels)
        subplot(2,length(uniqueObjLabels),obj_ind+length(uniqueObjLabels))
        singleObjMeans = plotArray{obj_ind,1};
        histogram(singleObjMeans, 20)
        title(sprintf('%s (\x03bc  = %s)',uniqueObjLabels{obj_ind}, num2str(round(mean(singleObjMeans), 2)) ));
        tmpValues = [tmpValues; singleObjMeans];
        tmpLabels = [tmpLabels; repmat(uniqueObjLabels(obj_ind) ,[length(singleObjMeans),1])];
      end
      subplot(2,length(uniqueObjLabels),[1:length(uniqueObjLabels)])
      boxplot(tmpValues, tmpLabels, 'Symbol', 'o')
      linkaxes(h.Children(2:length(h.Children)-1),'xy')
      savefig(fullfile(params.outputDir, figTitle))
      close(h)
    end 
  end
end

% Produces 36 Images
for obj_ind = 1:length(objListPlot)
  for group_ind = 1:length(groupingType)
    % Fig 2 - 'Mean Face Rates - MUA, All Stimuli'
    figTitle = sprintf('%s %s rates, %s, All Stimuli', dataType{2}, objListPlot{obj_ind}, groupingType{group_ind});
    figFilePath = fullfile(params.outputDir, figTitle);    
    if ~exist([figFilePath, '.fig'],'file')
      % Collect relevant data for plot
      stimObjData = objFrameFiringRatesTotal(:, group_ind, obj_ind,  2);
      plotArray = cell(length(uniqueStimLabels),2);
      removeInd = false(length(uniqueStimLabels),1);
      for stim_i = 1:length(uniqueStimLabels)
        stimMatInd = strcmp(uniqueStimLabels(stim_i), allStimuliVecBroad);
        plotArray{stim_i,1} = vertcat(stimObjData{stimMatInd});
        removeInd(stim_i) = isempty(plotArray{stim_i,1});
        plotArray{stim_i,2} = uniqueStimLabels{stim_i};
      end
      plotArray(removeInd,:) = [];
      
      if ~isempty(plotArray)
        % Set up grid of subplots correctly
        h = figure('Name', figTitle, 'NumberTitle', 'off', 'units','normalized','outerposition',[0 0 1 1]);
        sgtitle(figTitle)
        hold on
        stim2Plot = size(plotArray,1);
        [tmpValues, tmpLabels] = deal([]);
        for stim_ind = 1:stim2Plot
          subplot(2,stim2Plot,stim_ind+stim2Plot)
          stimObjMeans = plotArray{stim_ind,1};
          histogram(stimObjMeans, 20);
          title(sprintf('%s (\x03bc  = %s)',plotArray{stim_ind,2},num2str(round(mean(stimObjMeans), 2))));
          tmpValues = [tmpValues; stimObjMeans];
          tmpLabels = [tmpLabels; repmat(plotArray(stim_ind,2) ,[length(stimObjMeans),1])];
        end
        subplot(2,stim2Plot,1:stim2Plot)
        boxplot(tmpValues, tmpLabels, 'Symbol', 'o')
        linkaxes(h.Children(2:length(h.Children)-1),'xy')
        savefig(fullfile(params.outputDir, figTitle))
        close(h)
      end
    end
  end
end

% Convert to Arrays of run indices into tables
frameFiringStruct.channelUnitMat = channelUnitMat;
frameFiringStruct.objList = objListPlot;
frameFiringStruct.groupList = groupingType;
frameFiringStruct.stimList = uniqueStimLabels;
frameFiringStruct.dataList = dataType(dataInd2Plot);

end

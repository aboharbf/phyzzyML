function eyeDataStruct = calcEyeObjectTrace(eyeInByEvent, channelUnitNames, psthParams, eventIDs, taskData, eyeDataStruct, stimDir)
% Functions goal is to return a vector, of the same structure as
% "eyeInByEvent", replacing the {event}(1*Channels*Trials*Millseconds)
% structure with a {event}{trial}{frame*1}, where ever item is the
% name of what is being fixated on.

plotType = 'none'; %'trace'; %'area';%
PixelsPerDegree = taskData.eyeCal.PixelsPerDegree;
dvaOrigin =  taskData.eyeCal.origin;
frameMotionData = taskData.frameMotionData;
scaleFactor = taskData.scaleFactor;
frameMotionDataNames = {frameMotionData(:).stimVid}';

% psthParams.psthPre
% psthParams.psthImDur

shapeColors = {[1. 0. 0.];[0 .4 1.];[.1 .8 .1];[.1 .8 .1];[0 0 0];[1 .4 0]; ...
  [.7 0 0];[0 0 .7];[0 .5 0];[0 .5 0];[0 0 0];[1 .6 0]};

% Initialize vector of indicies which will be passed out.
[frameStartInd, frameEndInd, attendedObjVect] = deal(cell(length(eventIDs),1));

% Remove pre and post stimuli eye data
extractStimEye = @(x) x(:,:, psthParams.psthPre:psthParams.psthPre+psthParams.psthImDur);
eyeInByEvent = cellfun(extractStimEye, eyeInByEvent, 'UniformOutput', 0);
origInd = 1:psthParams.psthImDur+1;
stimPresSec = (psthParams.psthImDur/1000);

if all(contains(eventIDs, 'headTurnCon'))
  paradigm = 1;
elseif all(contains(eventIDs, 'headTurnIso'))
  paradigm = 2;
else
  paradigm = 3;
end


idVidDir = fullfile(stimDir, 'idVideos');
objArray = [{'Face1', 'HandR1', 'Body1'}; {'Face2', 'HandR2', 'Body2'}];

eyeInByEventDSPerEvent = cell(size(eventIDs));
stimFrameMotionData = struct;
for stim_i = 1:length(eventIDs)
  % find the correct frameMotionData
  stimFrameMotionData = frameMotionData(strcmp(frameMotionDataNames, eventIDs(stim_i)));
  
  % Grab the eye signal for an event, and down sample it to the frame
  % rate of the video, frames in video and frames in output don't match
  % because videos can be longer than they were presented for.
  frames  = round(stimFrameMotionData.fps)*stimPresSec;
  dsInd = linspace(origInd(1), origInd(end) - 1000/stimFrameMotionData.fps, frames+1);
  eyeInByEventDS = zeros(2, size(eyeInByEvent{stim_i}, 2), frames);
  sampFreq = 1/stimFrameMotionData.fps*1e3;
  
  frameStartInd{stim_i} = round((0:frames-1)*sampFreq)+1; % Create an index of when each frame starts
  frameEndInd{stim_i} = round((1:frames)*sampFreq);       % each index is the number of points averaged to make a frame.
  
  for trial_i = 1:size(eyeInByEvent{stim_i}, 2)
    for eye_i = 1:2
      tmp = interp1(origInd, squeeze(eyeInByEvent{stim_i}(eye_i, trial_i, :)), dsInd, 'linear');
      eyeInByEventDS(eye_i, trial_i, :) = tmp(1:end-1);
    end
  end
  
  % move data to the dVa Origin, then shift to pixel space
  pixelOrigin = [stimFrameMotionData.width/2 stimFrameMotionData.height/2];
  for eye_ind = 1:size(eyeInByEventDS,1)
    tmp = eyeInByEventDS(eye_ind, :, :) + dvaOrigin(eye_ind);
    eyeInByEventDS(eye_ind, :, :) = tmp * (PixelsPerDegree(eye_ind)/scaleFactor) + pixelOrigin(eye_ind);
  end
  
  eyeInByEventDSPerEvent{stim_i} = eyeInByEventDS;
  stimFrameMotionDataSorted(stim_i) = stimFrameMotionData;
end
  
for stim_i = 1:length(eventIDs)
  
  % If there is a match, grab that data and use it to downsample frames.
  stimFrameMotionData = stimFrameMotionDataSorted(stim_i);
  eyeInByEventDS = eyeInByEventDSPerEvent{stim_i};
  
  % Initialize the output cell array.
  attendedObjVect{stim_i} = cell(size(eyeInByEventDS, 2), size(eyeInByEventDS, 3));
  
  if strcmp(taskData.paradigm, 'naturalSocial') && ~contains(eventIDs(stim_i), 'monkey')
    
    % For natural videos without monkeys (Objects, Scenes, Scrambles)
    if ~isempty(stimFrameMotionData.objNames) % Landscapes/Scrambles
      
      % Make variables for each object. Not a clean way of doing this
      objects = stimFrameMotionData.objNames;
      for obj_i = 1:length(objects)
        eval(sprintf('%s = stimFrameMotionData.objLoc{%d};', objects{obj_i}, obj_i));
      end
      
      % Reshape the info to make it easy to compare
      objRads = stimFrameMotionData.objRadii;
      objLocStack = [stimFrameMotionData.objLoc{:}];
      
      for frame_i = 1:size(eyeInByEventDS, 3)
        %for each frame, pull and reshape the apporpriate shape coords and
        %compare to the eye signal across trials.
        objLoc = reshape(objLocStack(frame_i,:)', 2,12);
        eyeCoord = [eyeInByEventDS(1, :, frame_i); eyeInByEventDS(2, :, frame_i)];
        for trial_i = 1:length(eyeCoord)
          % Calculate distances
          trialEye = eyeCoord(:,trial_i);
          objDist = sqrt(sum((objLoc - trialEye).^2, 1));
          objInd = objDist < objRads;
          % Assign the correct object. None = bkg, 1 = object, 2 = closest.
          if sum(objInd) == 0
            attendedObjVect{stim_i}(trial_i, frame_i) = {'bkg'};
          elseif sum(objInd) == 1
            attendedObjVect{stim_i}(trial_i, frame_i) = stimFrameMotionData.objNames(objInd);
          elseif sum(objInd) > 1
            [~,I] = min(objDist(objInd));
            indArray = find(objInd);
            objInd = indArray(I);
            attendedObjVect{stim_i}(trial_i, frame_i) = stimFrameMotionData.objNames(objInd);
          end
        end
        
      end
    else
      attendedObjVect{stim_i}(:) = {'bkg'}; %Landscapes/Scrambles
    end
    
  else
    
    % For the headTurn paradigm, load the video
    if strcmp(taskData.paradigm, 'headTurnCon') || strcmp(taskData.paradigm, 'naturalSocial')
      %headTurnCon videos each have their own video (unless there are
      %no agents, in which case there are no labels.
      
      % find and load the ID video for the stimuli
      idVid = dir(fullfile(idVidDir, [extractBefore(eventIDs{stim_i}, '.avi') '*']));
      
      % If the idVid is empty, and the eventID contains the code for
      % scenes and objects, just skip this object.
      if isempty(idVid)
        idString = extractAfter(eventIDs{stim_i}, '15');
        intNum = str2num(idString(1));
        
        if intNum == 7 || intNum == 8
          attendedObjVect{stim_i}(:, :) = deal({'bkg'});
          continue
        else
          error('No ID Video found for %s', eventIDs{stim_i});
        end
        
      end
      
    elseif strcmp(taskData.paradigm, 'headTurnIso')
      %headTurnIso video ID generalize across skin
      idVidName = [extractBefore(eventIDs{stim_i}, '.avi') '_ID*'];
      %           idVidName(24) = num2str(1);
      
      idVid = dir(fullfile(idVidDir, idVidName));
      
    end
    
    % Load the video
    idVidObj = VideoReader(fullfile(idVid(1).folder, idVid(1).name));
    
    % run through and identify if the pixel where the eye signal was
    % looking at was black, red, green, or blue.
    frameSize = [idVidObj.Width idVidObj.Height];
    for frame_i = 1:size(attendedObjVect{stim_i}, 2)
      
      % Pull the corresponding frame of the ID Video
      frameData = read(idVidObj, frame_i);
      
      % Cycle through the trial and match whatever is being looked at.
      for trial_i = 1:size(attendedObjVect{stim_i}, 1)
        eyeDataInd = round(eyeInByEventDS(:, trial_i, frame_i))';
        
        % If the eyeData is out of the video (either negative or above
        % the frame size) then put in background.
        if any(eyeDataInd > frameSize) || any(eyeDataInd <= 0)
          attendedObjVect{stim_i}{trial_i, frame_i} = 'bkg';
          continue
        end
        pixelValue = squeeze(frameData(eyeDataInd(2), eyeDataInd(1), :))';
        
        % If one of the values is above the threshold, select its
        % object. Otherwise, remove it.
        if any(pixelValue > 100)
          
          % Identify object based on color, monkey based on shade.
          objInd = find(pixelValue > 100, 1);
          monkInd = (pixelValue(objInd) <= 150) + 1;
          attendedObjVect{stim_i}{trial_i, frame_i} =  objArray{monkInd, objInd};
          
        else
          
          % Otherwise, background
          attendedObjVect{stim_i}{trial_i, frame_i} = 'bkg';
          
        end
        
      end
    end
    
  end
  
end

% if some of the stimuli have a different frame count, resample the ones
% with fewer frames to bring them in line.
frameCounts = cellfun(@(x) size(x, 2), attendedObjVect);
if length(unique(frameCounts)) > 1
  maxFrame = max(frameCounts);
  updatedFramesInd = [];
  for stim_i = 1:length(attendedObjVect)
    if ~isempty(attendedObjVect{stim_i})
      if size(attendedObjVect{stim_i},2) < maxFrame
        updatedFramesInd = [updatedFramesInd stim_i];
        %for now, since I know this solution works for the videos I'm using,
        %I'll be doubling the vector, then chopping off any extra at the end.
        tmpArray = cell(size(attendedObjVect{stim_i},1), size(attendedObjVect{stim_i},2)*2);
        for frame_ind = 1:size(attendedObjVect{stim_i},2)
          tmpArray(:,frame_ind*2) = attendedObjVect{stim_i}(:,frame_ind);
          tmpArray(:,(frame_ind*2)-1) = attendedObjVect{stim_i}(:,frame_ind);
        end
        attendedObjVect{stim_i} = tmpArray(:,1:maxFrame);
      else
        frameIndSwapIndex = stim_i;
      end
    end
  end
  
  %Overwrite the frame indexes with previous written ones for a longer stimuli
  [frameStartInd{updatedFramesInd}] = deal(frameStartInd{frameIndSwapIndex});
  [frameEndInd{updatedFramesInd}] = deal(frameEndInd{frameIndSwapIndex});
  
end

% For each stimulus, return the number of frames spent looking at each
% label
simpleObjList = {'Face', 'Body', 'Hand', 'Object', 'bkg'};
[attendedPerStim, attendedPerStimRate] = deal(nan(length(eventIDs), length(simpleObjList)));

for stim_i = 1:length(attendedObjVect)
  stimFrames = length(attendedObjVect{stim_i}(:));
  for obj_i = 1:length(simpleObjList)
    attendedPerStim(stim_i, obj_i) = sum(sum(contains(attendedObjVect{stim_i}, simpleObjList{obj_i})));
    attendedPerStimRate(stim_i, obj_i) = attendedPerStim(stim_i, obj_i)/stimFrames;
  end
  
end

objAttendRates = struct();
objAttendRates.simpleObjList = simpleObjList;
objAttendRates.attendedPerStim = attendedPerStim;
objAttendRates.attendedPerStimRate = attendedPerStimRate;

%Plot the Result as an area plot, where X = frames
areaColors = [shapeColors{:} {[.5 .5 .5]}]'; %Shape colors + Background
%Hard coded for now, will need adjustment for other task sets.
objList = [{'Face1'},{'Body1'},{'HandL1'},{'HandR1'},{'GazeFollow1'},{'Object1'},{'Face2'},{'Body2'},{'HandL2'},{'HandR2'},{'GazeFollow2'},{'Object2'},{'bkg'}]';
tracePlotData = cell(1,length(eventIDs));

switch plotType
  case 'none'
    % Used in the case of wanting colorCodedRaster.
    tag2color = @(n) areaColors{strcmp(objList, n)}; %Returns appropriate color for each object
    for channel_i = 1:length(channelUnitNames)
      for event_i = 1:length(eventIDs)
        %For every event, grab the correct color
        attendedObjectEvent = attendedObjVect{event_i};
        if ~isempty(attendedObjectEvent)
          attendObjColorMat = cellfun(tag2color, attendedObjectEvent, 'UniformOutput', false); % a Trial * Frame array of the first letter of each object.
          tracePlotData{event_i} = zeros(size(attendObjColorMat, 1), size(attendObjColorMat, 2), 3); %Frames * Unique * colors objects
          % Give each frame its correct color
          for trial_i = 1:size(attendObjColorMat,1)
            for frame_i = 1:size(attendObjColorMat,2)
              tracePlotData{event_i}(trial_i, frame_i, :) = attendObjColorMat{trial_i, frame_i};
            end
          end
          
        end
      end
    end
  case 'area'
    for channel_i = 1:length(channelUnitNames)
      for event_i = 1:length(eventIDs)
        psthTitle = sprintf('Per Image Area eye plot Ch %d, %s', channel_i, eventIDs{event_i});
        figure('Name',psthTitle,'NumberTitle','off');
        subplot(2, 1, 1)
        %For every event, grab the relevant data and initialize the "counts"
        %matrix for each object.
        attendedObjectEvent = attendedObjVect{event_i};
        if ~isempty(attendedObjectEvent)
          
          areaPlotData = zeros(size(attendedObjectEvent, 2), length(objList)); %Frames * Unique objects
          for frame_ind = 1:size(attendedObjectEvent, 2)
            %For each frame, cycle through objects and see whats attended. Create a
            %total count of the 13 objects and the number of trials in which they
            %are attended to.
            frameObjAttended = attendedObjectEvent(:,frame_ind);
            [A, ~, C] = unique(frameObjAttended);
            a_counts = accumarray(C,1);
            for obj_ind = 1:length(A)
              dataInd = strcmp(objList, A{obj_ind});
              assert(sum(dataInd) == 1);
              areaPlotData(frame_ind, dataInd) = a_counts(obj_ind);
            end
          end
          areaPlotDataNorm = areaPlotData./size(attendedObjectEvent,1);
          areaPlotHandle = area(areaPlotDataNorm); %Plot normalized values
          %Make it presentable.
          title(eventIDs{event_i})
          xlim([1 size(attendedObjectEvent, 2)])
          ylim([0 1])
          for face_i = 1:length(areaPlotHandle)
            areaPlotHandle(face_i).FaceColor = areaColors{face_i};
          end
          legend([frameMotionData(1).objNames {'bkg'}])
          
          %Underneath, plot the activity of all identified units with a PSTH
          unitCount = length(channelUnitNames{channel_i});
          psthHandle = subplot(2, 1, 2);
          %Create a PSTH of responses to this event specifically
          eventPSTH = zeros(unitCount, stimEndInd-stimStartInd);
          %Labeling related code
          if unitCount > 2
            unitLabels = cell(unitCount,1);
            unitLabels{1} = 'Unsorted';
            unitLabels{end} = 'MUA';
            for unit_ind = 1:unitCount - 2
              unitLabels{unit_ind+1} = ['Unit ' num2str(unit_ind)];
            end
          else
            unitLabels = {'Unsorted', 'MUA'};
          end
          %Gather the correct data
          for unit_i = 1:unitCount
            eventPSTH(unit_i,:) = channelUnitNames{channel_i}{unit_i}(event_i,stimStartInd:stimEndInd-1);
          end
          psthParams.psthPre = 0;
          psthParams.psthPost = 0;
          plotPSTH(eventPSTH, [], psthHandle, psthParams, 'color', psthTitle, unitLabels);
        end
      end
    end
  case 'trace'
    handleArray = gobjects(length(channelUnitNames), length(eventIDs));
    tag2color = @(n) areaColors{strcmp(objList, n)}; %Returns appropriate color for each object
    for channel_i = 1:length(channelUnitNames)
      for event_i = 1:length(eventIDs)
        psthTitle = sprintf('Per Image Trace Ch %d, %s', channel_i, eventIDs{event_i});
        handleArray(channel_i, event_i) = figure('Name', psthTitle, 'NumberTitle','off','visible','off');
        subplot(2, 1, 1)
        %For every event, grab the correct color
        attendedObjectEvent = attendedObjVect{event_i};
        if ~isempty(attendedObjectEvent)
          attendObjColorMat = cellfun(tag2color, attendedObjectEvent, 'UniformOutput', false); % a Trial * Frame array of the first letter of each object.
          tracePlotData{event_i} = zeros(size(attendObjColorMat, 1), size(attendObjColorMat, 2), 3); %Frames * Unique * colors objects
          for trial_i = 1:size(attendObjColorMat,1)
            for frame_i = 1:size(attendObjColorMat,2)
              tracePlotData{event_i}(trial_i, frame_i, :) = attendObjColorMat{trial_i, frame_i};
            end
          end
          image(tracePlotData{event_i});
          %Make it presentable.
          title(eventIDs{event_i});
          legend([frameMotionData(1).objNames {'bkg'}]);
          %Underneath, plot the activity of all identified units with a PSTH
          unitCount = length(channelUnitNames{channel_i});
          psthHandle = subplot(2, 1, 2);
          %Create a PSTH of responses to this event specifically
          eventPSTH = zeros(unitCount, stimEndInd-stimStartInd);
          %Labeling related code
          if unitCount > 2
            unitLabels = cell(unitCount,1);
            unitLabels{1} = 'Unsorted';
            unitLabels{end} = 'MUA';
            for unit_ind = 1:unitCount - 2
              unitLabels{unit_ind+1} = ['Unit ' num2str(unit_ind)];
            end
          else
            unitLabels = {'Unsorted', 'MUA'};
          end
          %Gather the correct data
          for unit_i = 1:unitCount
            eventPSTH(unit_i,:) = channelUnitNames{channel_i}{unit_i}(event_i,stimStartInd:stimEndInd-1);
          end
          psthParams.psthPre = 0;
          psthParams.psthPost = 0;
          plotPSTH(eventPSTH, [], psthHandle, psthParams, 'color', psthTitle, unitLabels);
        end
      end
    end
    
    attendedObjData.handleArray = handleArray;
    
end

%Package outputs
attendedObjData.objAttendRates = objAttendRates;
attendedObjData.eventIDs = eventIDs;
attendedObjData.attendedObjVect = attendedObjVect;
attendedObjData.tracePlotData = tracePlotData;
attendedObjData.colorCode = areaColors';
attendedObjData.objList = objList;
attendedObjData.frameStartInd = frameStartInd;
attendedObjData.frameEndInd = frameEndInd;
eyeDataStruct.attendedObjData = attendedObjData;

end

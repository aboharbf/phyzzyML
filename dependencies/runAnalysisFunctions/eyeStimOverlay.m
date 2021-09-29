function [eyeInByEventDS, spikeEyeData] = eyeStimOverlay(eyeInByEvent, eyeDataStruct, spikesByEventBinned, eventIDs, taskData, psthParams, eyeStimOverlayParams)
%Function will visualize eye signal (and objects) on top of stimulus.

saccadeByStim = eyeDataStruct.saccadeByStim;
attendedObjData = eyeDataStruct.attendedObjData;

shapeOverlay = eyeStimOverlayParams.shapeOverlay;         % Switch for shape overlay.
trialNumberOverlay = eyeStimOverlayParams.trialNumberOverlay;   % Switch for trial number overlay on eye signal.
spikeOverlay = eyeStimOverlayParams.spikeOverlay;
colorCodedTrace = eyeStimOverlayParams.colorCodedTrace;      % Trace changes colors corresponding to what is being looked at. 1 = per target of gaze, 2 = Saccade vs fixation.
frameCountOverlay = eyeStimOverlayParams.frameCountOverlay;    % Places the frame number in the lower left.
eyeSig.shape = eyeStimOverlayParams.eyeSig.shape;  
eyeSig.color = eyeStimOverlayParams.eyeSig.color;
videoOutput = eyeStimOverlayParams.videoOutput;
stimDir = eyeStimOverlayParams.stimDir;
outDir = eyeStimOverlayParams.outDir;
computeEyeSpikeMat = eyeStimOverlayParams.computeEyeSpikeMat;

frameMotionData = taskData.frameMotionData;
frameMotionDataNames = {frameMotionData(:).stimVid}';

% Remove pre and post stimuli eye data, generate variables for eye sig
% parcing.
extractStimEye = @(x) x(:,:, psthParams.psthPre:psthParams.psthPre+psthParams.psthImDur);
extractSaccEye = @(x) x(:, psthParams.psthPre:psthParams.psthPre+psthParams.psthImDur);

eyeInByEvent = cellfun(extractStimEye, eyeInByEvent, 'UniformOutput', 0);
saccadeByStim = cellfun(extractSaccEye, saccadeByStim, 'UniformOutput', 0);

origInd = 1:psthParams.psthImDur+1;
stimPresSec = (psthParams.psthImDur/1000);
PixelsPerDegree = taskData.eyeCal.PixelsPerDegree;
dvaOrigin =  taskData.eyeCal.origin;
scaleFactor = taskData.scaleFactor;

% Initialize some structures
[eyeInByEventDS, eyeImgByEventDS] = deal(cell(size(eventIDs)));

% Generate smoothing kernal
Xsmooth = 10;         % Smoothing in the x dimensions (pixels)
Ysmooth = Xsmooth;    % Smoothing in the y dimensions (pixels)
Zsmooth = 1;          % Smoothing in the Z dimension (Frames)

smoothingFilterX = normpdf(-3*Xsmooth:3*Xsmooth, 0, Xsmooth);
smoothingFilterY = normpdf(-3*Ysmooth:3*Ysmooth, 0, Ysmooth);
smoothingFilterZ = normpdf(-3*Zsmooth:3*Zsmooth, 0, Zsmooth);

TDGauss = convn(smoothingFilterX, smoothingFilterY', 'full');
ZGauss = reshape(smoothingFilterZ, [1, 1, length(smoothingFilterZ)]);
ThreeDGauss = convn(TDGauss, ZGauss, 'full');

colormapStored = colormap;
chanCount = length(spikesByEventBinned{1});
unitCount = cellfun(@(x) length(x), spikesByEventBinned{1});
spikeEyeData = initNestedCellArray(spikesByEventBinned);

%Run each stimuli presented
for stim_i = 1:length(eventIDs)
  %find the correct frameMotionData
  frameMotionDataInd = strcmp(frameMotionDataNames, eventIDs(stim_i));
  stimFrameMotionData = frameMotionData(frameMotionDataInd);
  attendedObjStim = attendedObjData.attendedObjVect{stim_i};
  convert2Ind = @(x) find(strcmp(x, attendedObjData.objList));
  objIndex = cellfun(convert2Ind, attendedObjStim);
  
  stimVidStruct = dir([stimDir '/**/' eventIDs{stim_i}]);
  stimVidPath = [stimVidStruct(1).folder filesep stimVidStruct(1).name];
  stimVid = VideoReader(stimVidPath);
  
  % Grab the eye signal for an event, and down sample it to the frame rate of the video
  %eyeInByEventDS = downSampleSig(eyeInByEvent{stim_i}, psthParams, stimFrameMotionData);
  eyeImgByEvent = downSample1D(saccadeByStim{stim_i}, psthParams, stimFrameMotionData);
  
  frames  = min(round(stimFrameMotionData.fps)*stimPresSec, stimVid.NumFrames);
  framesIter = 1:frames;
  dsInd = linspace(origInd(1), origInd(end) - 1000/stimFrameMotionData.fps, frames+1);
  eyeInByEventDS{stim_i} = zeros(2, size(eyeInByEvent{stim_i},2), frames);
  eyeImgByEventDS{stim_i} = zeros(size(eyeInByEvent{stim_i},2), frames);
  
  % Down sample various signals.
  for trial_i = 1:size(eyeInByEvent{stim_i},2)
    % Downsample the eye labels
    for eye_i = 1:2
      tmp = interp1(origInd, squeeze(eyeInByEvent{stim_i}(eye_i, trial_i, :)), dsInd, 'linear');
      eyeInByEventDS{stim_i}(eye_i, trial_i, :) = tmp(1:end-1);
    end
    % Downsample saccade labels
    tmp = interp1(origInd, saccadeByStim{stim_i}(trial_i, :), dsInd, 'linear');
    eyeImgByEventDS{stim_i}(trial_i, :) = tmp(1:end-1);
    
    % Brief modification for visual clarity - in cases where saccades are a
    % single frame, expand.
    sacInd = find(eyeImgByEventDS{stim_i}(trial_i, :) == 2);
    for sac_i = 1:length(sacInd)
      if sacInd(sac_i) ~= 1 && sacInd(sac_i) ~= size(eyeImgByEventDS{stim_i},2)
        preFrame = eyeImgByEventDS{stim_i}(sacInd(sac_i)-1);
        postFrame = eyeImgByEventDS{stim_i}(sacInd(sac_i)+1);
        if preFrame == 1 && postFrame == 1
          [eyeImgByEventDS{stim_i}(trial_i, sacInd(sac_i)-1), eyeImgByEventDS{stim_i}(trial_i, sacInd(sac_i)+1)] = deal(2);
        end
      end
    end
    
  end
  
  % move data to the dVa Origin, then shift to pixel space
  pixelOrigin = [stimFrameMotionData.width/2 stimFrameMotionData.height/2];
  for eye_ind = 1:size(eyeInByEventDS{stim_i},1)
    tmp = eyeInByEventDS{stim_i}(eye_ind, :, :) + dvaOrigin(eye_ind);
    eyeInByEventDS{stim_i}(eye_ind, :, :) = tmp * (PixelsPerDegree(eye_ind)/scaleFactor) + pixelOrigin(eye_ind); % scaleFactor makes the videos larger, meaning functionally, fewer pixels of video per degree.
  end
    
  %Shape related info
  if shapeOverlay
    objects = stimFrameMotionData.objNames;
    objectShapes = stimFrameMotionData.objShapes;
    objRad = stimFrameMotionData.objRadii;
    for obj_i = 1:length(objects)
      eval(sprintf('%s = stimFrameMotionData.objLoc{%d};', objects{obj_i}, obj_i));
    end
  end
  
  % Select the matrix to use for color
  if colorCodedTrace == 1
    colorIndMat = objIndex;%(shapeColors{obj_i}.*255);
  elseif colorCodedTrace == 2
    colorIndMat = eyeImgByEvent;
  else
    warning('colorCodedTrace choice wrong, using 1');
    colorIndMat = ones(size(objIndex));
  end
  
  % In cases where things were up sampled for storage, down sample back to
  % the appropriate frame count.
  if computeEyeSpikeMat
    if ~(size(eyeInByEventDS{stim_i},3) == size(colorIndMat,2))
      nInd = 1:size(eyeInByEventDS{stim_i},3);
      oInd = 1:size(colorIndMat,2);
      tmp = interp1(oInd, colorIndMat', nInd ,'nearest', 'extrap');
      colorIndMat = round(tmp');
    end
    
    % Create a 3D Matrix combining binned spike information and eye signal
    framesShown = ceil((psthParams.psthImDur/1000) * stimVid.FrameRate);
    spikeEyeMatrixTemplate = zeros(stimVid.Width, stimVid.Height, framesShown);
    
    activityShift = psthParams.psthPre + psthParams.movingWin(1)/2;     % Shift to stimulus activity.
    frameStarts = attendedObjData.frameStartInd{stim_i} + activityShift;
    frameEnds = attendedObjData.frameEndInd{stim_i} + activityShift;
    spikeBinnedStim = spikesByEventBinned{stim_i};
    
    % Cycle through chans and units.
    for chan_i = 1:chanCount
      for unit_i = 1:unitCount(chan_i)
        spikeEyeMat = spikeEyeMatrixTemplate;
        
        % If there are no spikes, go on to next unit.
        if all(all(spikeBinnedStim{chan_i}{unit_i} == 0))
          continue
        end
        
        % Cycle through frames
        for frame_i = framesIter
          % Bin the spikes which occur during this frame, accorinding to previous
          % defined times.
          spikesPerBinPerTrial = sum(spikeBinnedStim{chan_i}{unit_i}(:, frameStarts(frame_i):frameEnds(frame_i)), 2);
          
          % Place the spikes on the spikeEyeMat
          for trial_i = find(spikesPerBinPerTrial)'
            coords = round(eyeInByEventDS{stim_i}(:, trial_i, frame_i));
            if any(coords' > size(spikeEyeMat, 1:2)) || any(coords' < 0)
              continue
            else
              spikeEyeMat(coords(1), coords(2), frame_i) = spikeEyeMat(coords(1), coords(2), frame_i) + spikesPerBinPerTrial(trial_i);
            end
          end
          
        end
        
        % Convolve the stack with a 3D gaussian.
        spikeEyeMat = convn(spikeEyeMat, ThreeDGauss, 'same');
        
        % store the eye data, as a sparse matrix.
        spikeEyeMatR = reshape(spikeEyeMat, [stimVid.Width, stimVid.Height * framesShown]);
        spikeEyeMatRS = sparse(spikeEyeMatR);
        
        spikeEyeData{stim_i}{chan_i}{unit_i} = spikeEyeMatRS;
      end
    end
  else
    spikeEyeData = [];
  end
  
  % If Desired, create video of output
  if videoOutput
    
    if spikeOverlay
      outputVideo = cell(length(spikeEyeData), 1);
    else
      %Open a new video to save the results (video gets opened at end for
      %spikeOverlay).
      vidNameNoSpike = [outDir 'onlay_' stimVidStruct(1).name];
      outputVideo{1} = VideoWriter(vidNameNoSpike);
      outputVideo{1}.FrameRate = stimVid.FrameRate;
      open(outputVideo{1});
    end
    
    % Initialize a blend object, create indexes for better handling
    % sparse indices
    if spikeOverlay
      blend = vision.AlphaBlender('Operation', 'blend');
      sparseIndEnds = stimVid.Height:stimVid.Height:size(spikeEyeMatRS, 2);
      sparseIndStarts = [1 sparseIndEnds(1:end)+1];
    end
    
    for frame_i = framesIter
      
      % Pull a frame from the video
      img1 = read(stimVid, frame_i);
      
      % Add in Shapes (If the data exists, and switch is set)
      if shapeOverlay && ~isempty(stimFrameMotionData.objNames) %Landscapes/Scrambles
        for obj_i = 1:length(objects)
          coords = eval(eval(sprintf('objects{%d}',obj_i)));
          coords = round(coords(frame_i,:));
          if any(isnan(coords)) %NaN means the object isn't present.
            continue
          end
          switch objectShapes{obj_i}
            case 'Circle'
              img1 = insertShape(img1, objectShapes{obj_i}, [coords(1) coords(2) objRad(obj_i)], 'LineWidth', 2, 'color', (eyeSig.color{1}{obj_i}.*255));
            case 'Line'
              %the 2 lines define gazes - the first one comes from Face1, the
              %2nd from Face2.
              if obj_i == find(strcmp(objectShapes, 'Line'), 1, 'first')
                GazeSource = Face1;
              else
                GazeSource = Face2;
              end
              img1 = insertShape(img1,objectShapes{obj_i},[GazeSource(frame_i,1) GazeSource(frame_i,2) coords(1) coords(2)],'LineWidth',2);
          end
        end
      end
      
      % Add the eye signal, color coded according to the switch
      
      for trial_i = 1:size(eyeInByEventDS{stim_i}, 2)
        coords = round(eyeInByEventDS{stim_i}(:, trial_i, frame_i));
        % Place the shapes corresponding to gaze for every trial
        %       if colorCodedTrace ~= 2
        img1 = insertShape(img1, eyeSig.shape, [coords(1) coords(2) 1],'LineWidth', 10, 'Color', eyeSig.color{colorCodedTrace}{colorIndMat(trial_i, frame_i)}.*255);
        %img1 = insertShape(img1, eyeSig.shape, [coords(1) coords(2) 1],'LineWidth', 10, 'Color', 'r');
        %       else
        %         img1 = insertShape(img1, eyeSig.shape,[coords(1) coords(2) 1],'LineWidth',5,'Color',[eyeSig.color{colorCodedTrace}{colorIndMat(trial_i, frame_i)}]);
        %       end
        % If desired, put number on circle
        if trialNumberOverlay
          img1 = insertText(img1, [coords(1)-5 coords(2)-7], num2str(trial_i), 'BoxOpacity', 0, 'FontSize', 8);
        end
        
      end
      
      % Add the Frame to the lower left of the stimuli
      if frameCountOverlay
        img1 = insertText(img1,[10 255], num2str(frame_i), 'BoxOpacity' , 0,'FontSize', 20, 'TextColor', 'yellow');
      end
      
      if spikeOverlay
        % Extract spikeDataMat, Create an image with spikeData overlaid.
        for chan_i = 1:length(spikeEyeData{stim_i})
          for unit_i = length(spikeEyeData{stim_i}{chan_i})
            
            % Return the spikeEyeData to full, indexing to find the correct frame.
            spikeEyeMatRS = spikeEyeData{stim_i}{chan_i}{unit_i};
            spikeFrameInfo = full(spikeEyeMatRS(:, sparseIndStarts(frame_i):sparseIndEnds(frame_i))');
            
            if ~all(all(spikeFrameInfo == 0))
              % If there are spikes, create an image for them and blend.
              colorScaler = max(max(spikeFrameInfo))/256;
              img2 = uint8(ind2rgb(uint8(spikeFrameInfo/colorScaler), colormapStored) * 255);
              
              AlphaDataScale = spikeFrameInfo/max(max(spikeFrameInfo));
              % Combine the images
              blend.Opacity = AlphaDataScale;
              imgBlend = blend(img1, img2);
            else
              imgBlend = img1;
            end
            
            %Open a new video to save the results
            if frame_i == 1
              vidPath = fullfile(outDir, sprintf('onlay_Ch%d_U%d_%s', chan_i, unit_i, stimVidStruct(1).name));
              outputVideo{chan_i} = VideoWriter(vidPath);
              outputVideo{chan_i}.FrameRate = stimVid.FrameRate;
              open(outputVideo{chan_i});
            end
            
            % Write the frame
            writeVideo(outputVideo{chan_i}, imgBlend);  % Add the new Frame
          end
        end
      else
        % Write the frame
        writeVideo(outputVideo{1}, img1);  % Add the new Frame
        
        % To get example frames, have them saved in the same dir
        if 1
          imgName = sprintf('%s_Frame%d.png', extractBefore(vidNameNoSpike, '.avi'), frame_i);
          imwrite(img1, imgName);
        end
        
      end
      
    end
    
    for ii = 1:length(outputVideo)
      close(outputVideo{ii});
    end
  end
  
end

end
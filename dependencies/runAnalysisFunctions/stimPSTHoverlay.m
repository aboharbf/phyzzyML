function stimPSTHoverlay(psthByImage, sortMask, inclusionMask, stimDir, psthParams, lfpPaddedBy, taskEventList, outDir)
%Rearrange PSTH due to sorting which takes place w/ signifiance bars
%Creates a copy of the video of the stimulus with the PSTH traced below.
psthPre = psthParams.psthPre;
psthImDur = psthParams.psthImDur;
psthPost = psthParams.psthPost;
times = -psthPre:psthImDur+psthPost;
stimStartInd = psthPre+lfpPaddedBy + 1;
stimEndInd = stimStartInd + psthImDur - 1;
groupingType = {'Unsorted','Unit','MUA'};

for channel_i = 1:length(psthByImage)
  for unit_i = 1:length(psthByImage{channel_i})
    unitPSTH = psthByImage{channel_i}{unit_i};
    yMax = max(max(unitPSTH));
    if yMax == 0
      yMax = 1; % 0 spike unit protection.
    end
    yMin = min(min(unitPSTH));
    
    if unit_i == 1
      unitLabel = groupingType{1};
    elseif unit_i == length(psthByImage{channel_i})
      unitLabel = groupingType{3};
    else
      unitLabel = ['U' num2str(unit_i - 1)];
    end
    
    for epoch_i = 1:length(inclusionMask)
      
      runMask = (inclusionMask{epoch_i}{channel_i}{unit_i} < 0.05);   % Find out which stimuli produced the significant activity
      if sum(runMask) ~= 0
        
        sortVec = sortMask{epoch_i}{channel_i}{unit_i};
        stimOfInterest = sortVec(runMask);
        PSTHtoPlot = unitPSTH(stimOfInterest,:);       % Recover vectors of interest
        stimtoPlot = taskEventList(stimOfInterest);    % Sort the eventList
        
        for ii = 1:length(stimtoPlot)
          
          %Isolate the stimuli, its name, and path.
          stimInfo = dir(strcat(stimDir, '/**/', stimtoPlot{ii})); %use the taskEventList to find the file
          stimPath = [stimInfo(1).folder filesep stimInfo(1).name]; %create its path.
          
          %Open the Video, Get some Info on it
          stimVid = VideoReader(stimPath);
          
          %Draw the PSTH video to be added.
          psthTrace = PSTHtoPlot(ii, psthPre:psthPre+psthImDur);
          timeInd = 1:psthImDur;
          PSTHVideoPath =  [outDir 'tmpPSTH_' stimInfo(1).name];
          outputVideoPath = fullfile(outDir, sprintf('PSTH_Ch%d_%s_%s', channel_i, unitLabel, stimInfo(1).name));
          
          % new video
          PSTHVideo = VideoWriter(PSTHVideoPath);
          PSTHVideo.FrameRate = stimVid.FrameRate;
          framesPerSecond = round(1000/stimVid.FrameRate);
          open(PSTHVideo);
          
          %Initialize the line, make the figure plot the right size for later
          %adjoining to the stimulus.
          currFig = figure();
          an = animatedline('color',[1 0 0],'LineWidth',2);
          currFig.Children.Color = [0 0 0];
          currFig.Children.Units = 'pixels';
          currFig.Children.Position = [0 0 stimVid.Width stimVid.Height/5];
          currFig.Position(3:4) = [stimVid.Width stimVid.Height/5];
          axis([timeInd(1), timeInd(end), yMin, yMax])
          
          for time_ind = 1:timeInd(end)
            if mod(time_ind,framesPerSecond) == 0
              addpoints(an, timeInd(time_ind), psthTrace(time_ind));
              frame = getframe(gcf);
              writeVideo(PSTHVideo, frame.cdata);
            end
          end
          
          close(currFig);
          close(PSTHVideo);
          
          %Open the stimulus video, play it frame by frame with the PSTH, and save
          %as a new movie.
          
          PSTHVideo = VideoReader(PSTHVideoPath);
          
          %videoPlayer = vision.VideoPlayer;
          outputVideo = VideoWriter(outputVideoPath);
          outputVideo.FrameRate = stimVid.FrameRate;
          open(outputVideo);
          
          while hasFrame(stimVid) && hasFrame(PSTHVideo)
            img1 = readFrame(stimVid);
            img2 = readFrame(PSTHVideo);
            imgt = vertcat(img1, img2);
            % play video
            %step(videoPlayer, imgt);
            % record new video
            writeVideo(outputVideo, imgt);
          end
          
          %release(videoPlayer);
          close(outputVideo);
          clear PSTHVideo
          delete(PSTHVideoPath);
        end
      end
      
    end
  end
end


end

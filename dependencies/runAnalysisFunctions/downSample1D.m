function [downSampledEyeInByEvent,frameStartInd, frameEndInd] = downSample1D(eyeInByEvent, psthParams, stimFrameMotionData)
%Function which down samples a signal (assumed to be 3D) present in the 3rd
%dimension of the input vector (2 eyes * trials * samples for eye signal)
%by averaging across samples present between points of interest.

frames  = round(stimFrameMotionData.fps*(psthParams.psthImDur/1000));
sampFreq = psthParams.psthImDur/frames;
clear frameCountArray frameIndArray
frameCountArray = diff(round((1:frames)*sampFreq)); %each index is the number of points averaged to make a frame.
frameCountArray = [frameCountArray frameCountArray(end)]; %Diff chops off the end, so repeat it
frameIndArray = ones(frameCountArray(1), 1);
for avg_ind = 2:length(frameCountArray)
  frameIndArray = [frameIndArray; ones(frameCountArray(avg_ind), 1)*avg_ind];
end

% Downsample the signal
eyeInAvg = zeros(size(eyeInByEvent,1), frames);
for frame_i = 1:frames
  eyeInAvg(:,frame_i) = round(mean(eyeInByEvent(:,frameIndArray == frame_i),2));
end
downSampledEyeInByEvent = eyeInAvg;

%Frames to ms, for other functions to easily compare to 1kS/sec data.
frameStartInd = round((0:frames-1)*sampFreq)+1; %Create an index of when each frame starts
frameEndInd = round((1:frames)*sampFreq); %each index is the number of points averaged to make a frame.
end

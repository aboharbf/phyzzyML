function [newPacketData, newPacketTimes] = replacePacketData(monkeyLogicData, TrialRecord, oldPacketData, oldPacketTimes)
% A function intended to replace the packetData and packetTimes in the
% event that they don't line up with MonkeyLogic.

% Written to deal with issue between MKL and BKL machine. Applies to files
% record on and between 6/21/2021 and 6/25/2021. Resolved after by switch
% monkeylogic to 'send and clear' setting, which does not rely on strobe.

% Stack all the monkeyLogic times and data
trialStructs = [monkeyLogicData.BehavioralCodes];
allTrialCodes = vertcat(trialStructs.CodeNumbers);

% Shift all monkeylogic times to be relative to run start.
absCodeTimes = [];
for trial_i = 1:length(trialStructs)
  newTimes = trialStructs(trial_i).CodeTimes + monkeyLogicData(trial_i).AbsoluteTrialStartTime;
  absCodeTimes = [absCodeTimes ; newTimes];
end

% See that the trial starts are all present, and use these to find the time offset.
% Depending on the date, the strategy for finding the off times is
% different. on and before 0625, it is just finding 9s. After (if this code
% triggered), the 9's are 8's.

trialStartTimesMkl = absCodeTimes(allTrialCodes == 9);

[~, fileName, ~] = fileparts(TrialRecord.DataFile);
if contains(fileName, {'0621', '0622', '0623', '0624', '0625', '0627'})
  % 9's are preserved, rest of code is off.
  trialStartTimesBkl = oldPacketTimes(oldPacketData == 9);
elseif contains(fileName, {'0628', '0629', '0630'})
  % After fixing the first issue, it seems the 1 bit started having an
  % issue the day after. 
  trialStartTimesBkl = oldPacketTimes(oldPacketData == 8);
end

% A couple days when this happens, the normal alignment doesn't cut it. In
% those cases, use a subset of the trial starts to find the proper offset.
switch fileName
  case {'20210625Sam001', '20210625Sam004', '20210628Sam001', '20210628Sam002', '20210628Sam004', '20210630Sam006', '20210629Sam001', '20210630Sam001'}
    % '20210625Sam001' - Blackrock missed the first and some marker in the
    % middle
    % '20210625Sam004' - Blackrock missed the last 9.
    % '20210630Sam006' - 2 8's sent for some reason, messes up Algo.
    % Both procedures work for this purpose.
    % Vectors below make it clear whats the alignment issue in each case.
    mklDiff = round(diff(trialStartTimesMkl), -1);
    bklDiff = round(diff(trialStartTimesBkl), -1);
    lag = finddelay(mklDiff, bklDiff);
    if lag > 0
      % Bkl has extra
      trialStartTimesMkl = trialStartTimesMkl(1:50);
      trialStartTimesBkl = trialStartTimesBkl(2:51);
    elseif lag < 0
      % Mkl has extra
      trialStartTimesMkl = trialStartTimesMkl(2:51);
      trialStartTimesBkl = trialStartTimesBkl(1:50);
    elseif lag == 0
      % there is an extra bit somewhere.
      trialStartTimesMkl = trialStartTimesMkl(1:50);
      trialStartTimesBkl = trialStartTimesBkl(1:50);
    end

    % Check it worked
    comparisonVec = [round(diff(trialStartTimesMkl), -1), round(diff(trialStartTimesMkl), -1)];
    assert(all(~diff(comparisonVec, 1, 2)));
    
  otherwise
    
    % Make sure they're both the same date.
    assert(length(trialStartTimesMkl) == length(trialStartTimesBkl), 'Trial start counts dont match in MKL and BKL Data. find and verify alignment technique.')
    
    % A double check on the lines where the numbers match. This may trigger
    % needlessly in cases where strobe is sent right at the edges.
    differenceVector = [round(diff(trialStartTimesMkl), -1), round(diff(trialStartTimesMkl), -1)];
    assert(all(~diff(differenceVector, 1, 2)))

end

% Use a linear model to find the offset between the monkeyLogic collected
% times and blackrock's. Shift the monkeyLogic times to make them blackrock times
logVsBlkModel = fitlm(trialStartTimesMkl, trialStartTimesBkl);  % See how X (the mkl times) can be shifted to y (bkl times).
taskEventStartTimesFit = predict(logVsBlkModel, absCodeTimes); % Use the model to shift monkeyLogic times into blackrock space.

% Return new vectors
newPacketData = allTrialCodes;
newPacketTimes = taskEventStartTimesFit;

end
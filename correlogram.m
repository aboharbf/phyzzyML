function [ Cgram, C, shifts, confCgram, confC ] = correlogram( data1, data2, varargin )
%correlogram returns the tf correlogram of data1 and data2 
%   inputs: 
%     - data1, data2: nTrials x nSamples arrays
%     - varargin (for default values, see code):
%       - maxShift: bound for C computation
%       - matchTimeRanges: if 1, use the same time range for each entry in C
%       - useJacknife: if 0, use theoretical errorbars assuming independent samples
%       - normalize: if 1, return correlation coefficients; else, return raw pointwise product
%       - timeDifferenceBound: [minDiff, maxDiff], only compute Cgram entries in range [currently not implemented]
%       - timeAverageRange: 
%       - jacknifeDraws: default is true jacknife
%       - jacknifeTrialsPerDraw: defualt is true jacknife
%       - jacknifeParallelWorkers: if > 0, parallelize the jacknife loop
%       - confidenceThreshold: float in [0 100], % of distribution that +/- errorbars enclose
%       - usePercentileConf: bool, if 1, get confidence bound from empirical jacknife distribution 
%
%   outputs:
%     - Cgram: nSamples x nSamples array
%     - C: time-averaged autocorrelation
%     - shifts: shift times for C. Shift refers to data2, so positive
%       correlation at negative shift means data2 leads data1
%     - CramConf: + and - errorbars for Cgram
%     - confC: + and - errorbars for C
%   notes:
%     - error calculations for autocorrelation are off
if ~isempty(varargin)
  params = varargin{1};
else
  params.dummy = 1;
end
if isfield(params,'maxShift') && ~isempty(params.maxShift)
  maxShift = params.maxShift;
else
  maxShift = ceil(size(data1,2)/4);  %this means we'll include the middle half of the data, if matchTimeRanges == 1
end
if isfield(params,'matchTimeRanges') && ~isempty(params.matchTimeRanges)
  matchTimeRangeAllShifts = params.matchTimeRanges;
else
  matchTimeRangeAllShifts = 1;
end
if isfield(params, 'useJacknife') && ~isempty(params.useJacknife)
  useJacknife = params.useJacknife;
else
  useJacknife = 0;
end
if isfield(params,'normalize') && ~isempty(params.normalize)
  normalize = params.normalize;
else
  normalize = 1;
end
if isfield(params,'timeDifferenceBound') && ~isempty(params.timeDifferenceBound)
  timeDifferenceBound = params.timeDifferenceBound;
else
  timeDifferenceBound = [0,size(data1,2)];
end
if isfield(params,'timeAverageRange') && ~isempty(params.timeAverageRange)
  timeAverageStartInd = params.timeAverageRange(1);
  timeAverageEndInd = params.timeAverageRange(2);
else
  timeAverageStartInd = 1;
  timeAverageEndInd = size(data1,2);
end
if isfield(params, 'jacknifeDraws') && ~isempty(params.jacknifeDraws)
  jacknifeDraws = min(params.jacknifeDraws, size(data1,1));
else
  jacknifeDraws = size(data1,1);
end
if isfield(params, 'jacknifeTrialsPerDraw') && ~isempty(params.jacknifeTrialsPerDraw)
  jacknifeTrialsPerDraw = params.jacknifeTrialsPerDraw;
else
  jacknifeTrialsPerDraw = size(data1,1)-1;
end
if isfield(params, 'jacknifeParallelWorkers') && ~isempty(params.jacknifeParallelWorkers)
  jacknifeParallelWorkers = params.jacknifeParallelWorkers;
else
  jacknifeParallelWorkers = 0; %note that in matlab parallel worker counting, 0 means serial
end
if isfield(params, 'confidenceThreshold') && ~isempty(params.confidenceThreshold)
  tcrit = tinv(params.confidenceThreshold,1);
else
  tcrit = 1; %default error bar is +/- 1 std
  confidenceThreshold = 0.75; %tcdf(1,1) = 0.75
end
if isfield(params, 'usePercentileConf') && ~isempty(params.usePercentileConf)
  usePercentileConf = params.usePercentileConf;
else
  usePercentileConf = 0;
end
if isfield(params, 'useSmoothing') && ~isempty(params.useSmoothing)
  useSmoothing = params.useSmoothing;
else
  useSmoothing = [0, 0];
end
if isfield(params, 'smoothingFilter') && ~isempty(params.smoothingFilter)
  smoothingFilter = params.smoothingFilter;
else
  usePercentileConf = 0;
end
assert(all(size(data1) == size(data2)),'input arrays must be the same size');
if useSmoothing(1)
  data1 = convn(data1,smoothingFilter,'same');
end
if useSmoothing(2)
  data2 = convn(data2,smoothingFilter,'same');
end
Cgram = zeros(size(data1,2));
CgramErr = zeros(size(data1,2));
confCgram = zeros(size(data1,2),size(data1,2),2);
if normalize
  data1 = data1 - mean(data1,1);
  data2 = data2 - mean(data2,1);
end
numTrials = size(data1,1);
if ~normalize
  if ~useJacknife
    Output.DEBUG('using non-jacknife, non-normalized computation');
    for t1 = 1:size(data1,2)
      Cgram(t1,:) = mean(data1(:,t1).*data2,1);
    end
    CgramErr = CgramErr*((3/2)*sqrt(2)/numTrials);
  else
    Output.DEBUG('using jacknife, non-normalized computation');
    for t1 = 1:size(data1,2)
      for t2 = 1:size(data2,2)
        [CgramEntry, CgramStd, confCgramEntry] = jacknifeHelper([data1(:,t1),data2(:,t2)],@(d) mean(d(:,1).*d(:,2)),'draws',jacknifeDraws,...
          'trialsPerDraw',jacknifeTrialsPerDraw,'confPercentile',usePercentileConf*confidenceThreshold, 'parallelWorkers',jacknifeParallelWorkers,'separableF',1);
        Cgram(t1,t2) = CgramEntry;
        CgramErr(t1,t2) = CgramStd;
        confCgram(t1,t2,:) = confCgramEntry;
      end
    end
  end
else
  data1stds = std(data1,[],1);
  data2stds = std(data2,[],1);
  if ~useJacknife
    Output.DEBUG('using non-jacknife, normalized computation');
    for t1 = 1:size(data1,2)
      Cgram(t1,:) = mean(data1(:,t1).*data2,1)./(data1stds(t1)*data2stds);
    end
    CgramErr = CgramErr*((3/2)*sqrt(2)/numTrials);
  else
    Output.DEBUG('using jacknife, normalized computation');
    for t1 = 1:size(data1,2)
      for t2 = 1:size(data2,2)
        [CgramEntry, CgramStd, confCgramEntry] = jacknifeHelper([data1(:,t1),data2(:,t2)],@(d) mean(d(:,1).*d(:,2))/(std(d(:,1))*std(d(:,2))),'draws',jacknifeDraws,...
          'trialsPerDraw',jacknifeTrialsPerDraw,'confPercentile',usePercentileConf*confidenceThreshold, 'parallelWorkers',jacknifeParallelWorkers,'separableF',0);
        Cgram(t1,t2) = CgramEntry;
        CgramErr(t1,t2) = CgramStd;
        confCgram(t1,t2,:) = confCgramEntry;
      end
    end
  end
end


C = zeros(2*maxShift+1,1);
confC = zeros(2*maxShift+1,2);
shifts = -maxShift:1:maxShift;
for shift_i = 1:length(shifts)
  shift = shifts(shift_i);
  correlsAllTimes = diag(Cgram,shift);
  if ~useJacknife || ~usePercentileConf
    correlErrsAllTimes = diag(CgramErr,shift);
  else
    confUpperAllTimes = diag(squeeze(confCgram(:,:,1)),shift);
    confLowerAllTimes = diag(squeeze(confCgram(:,:,2)),shift);
  end
  if matchTimeRangeAllShifts && timeAverageStartInd <= 1 && timeAverageEndInd >= size(data1,2)
    % slice the diagonal so that we use the same time range for all shifts
    C(shift_i) = mean(correlsAllTimes(maxShift-abs(shift)+1:end-(maxShift-abs(shift))));
    if ~useJacknife || ~usePercentileConf
      Cstd = norm(correlErrsAllTimes(maxShift-abs(shift)+1:end-(maxShift-abs(shift))).^2)/sqrt(length(correlsAllTimes(maxShift-abs(shift)+1:end-(maxShift-abs(shift)))));
      confC(shift_i,1) = tcrit*Cstd;
      confC(shift_i,2) = tcrit*Cstd;
    else
      confC(shift_i,1) = norm(confUpperAllTimes(maxShift-abs(shift)+1:end-(maxShift-abs(shift))))/size(data1,1)/sqrt(length(correlsAllTimes(maxShift-abs(shift)+1:end-(maxShift-abs(shift)))));
      confC(shift_i,2) = norm(confLowerAllTimes(maxShift-abs(shift)+1:end-(maxShift-abs(shift))))/size(data1,1)/sqrt(length(correlsAllTimes(maxShift-abs(shift)+1:end-(maxShift-abs(shift)))));
    end
  else
    if timeAverageStartInd > 1 || timeAverageEndInd < size(data1,2)
      C(shift_i) = mean(correlsAllTimes(timeAverageStartInd:timeAverageEndInd));
      if ~useJacknife || ~usePercentileConf
        Cstd = norm(correlErrsAllTimes(timeAverageStartInd, timeAverageEndInd).^2)/sqrt(timeAverageEndInd-timeAverageStartInd+1);
        confC(shift_i,1) = tcrit*Cstd;
        confC(shift_i,2) = tcrit*Cstd;
      else
        confC(shift_i,1) = norm(confUpperAllTimes(timeAverageStartInd:timeAverageEndInd))/sqrt(timeAverageEndInd-timeAverageStartInd+1);
        confC(shift_i,2) = norm(confLowerAllTimes(timeAverageStartInd:timeAverageEndInd))/sqrt(timeAverageEndInd-timeAverageStartInd+1);
      end
    else
      C(shift_i) = mean(correlsAllTimes);
      if ~useJacknife || ~usePercentileConf
        Cstd = norm(correlErrsAllTimes.^2)/sqrt(length(correlErrsAllTimes));
        confC(shift_i,1) = tcrit*Cstd;
        confC(shift_i,2) = tcrit*Cstd;
      else
        confC(shift_i,1) = norm(confUpperAllTimes)/sqrt(length(correlErrsAllTimes));
        confC(shift_i,2) = norm(confLowerAllTimes)/sqrt(length(correlErrsAllTimes));
      end
    end
  end
end
end


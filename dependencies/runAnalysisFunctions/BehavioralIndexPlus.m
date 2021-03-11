function [behaviortime] = BehavioralIndexPlus(blinkInds, vecLength)
%function turns indexes into times by parsing at breaks in continuity -

behavVec = zeros(vecLength,1);
behavVec(blinkInds) = deal(1);
behavMat = {find(behavVec)', find(~behavVec)'};
behaviortimeArray = cell(length(behavMat),1);

for ii = 1:length(behavMat)
  behavind = behavMat{ii};
  %function turns indexes into times by parsing at breaks in continuity
  %function turns indexes into times by parsing at breaks in continuity -
  % taken from ClusterFix
  dind = diff(behavind);
  gaps =find(dind > 1);
  behaveind = zeros(length(gaps),50);
  if ~isempty(gaps)
    for gapind = 1:length(gaps)+1
      if gapind == 1
        temp = behavind(1:gaps(gapind));
      elseif gapind == length(gaps)+1
        temp = behavind(gaps(gapind-1)+1:end);
      else
        temp = behavind(gaps(gapind-1)+1:gaps(gapind));
      end
      behaveind(gapind,1:length(temp)) = temp;
    end
  else
    behaveind =  behavind;
  end
  behaviortime = zeros(2,size(behaveind,1));
  for index=1:size(behaveind,1)
    rowfixind = behaveind(index,:);
    rowfixind(rowfixind == 0) = [];
    behaviortime(:,index) = [rowfixind(1);rowfixind(end)];
  end
  
  behaviortimeArray{ii} = behaviortime;
  
end


% If intervals are too short, remove them, fuse events.
IntervalMin = 15;
nonBlinkTimes = behaviortimeArray{2};
InterBlinkInterval = diff(nonBlinkTimes,1);
interval2Remove = InterBlinkInterval < IntervalMin;

if any(interval2Remove)
  resortInd = [];
  for ii = 1:length(interval2Remove)
    if interval2Remove(ii)
      resortInd = [resortInd, nonBlinkTimes(1,ii):nonBlinkTimes(2,ii)];
    end
  end
  
  behavMat{1} = sort([behavMat{1}, resortInd]);
  
  for ii = 1:length(behavMat)
    behavind2 = behavMat{ii};
    % Identify times for both
    dind = diff(behavind2);
    gaps =find(dind > 1);
    behaveind = zeros(length(gaps),50);
    if ~isempty(gaps)
      for gapind = 1:length(gaps)+1
        if gapind == 1
          temp = behavind2(1:gaps(gapind));
        elseif gapind == length(gaps)+1
          temp = behavind2(gaps(gapind-1)+1:end);
        else
          temp = behavind2(gaps(gapind-1)+1:gaps(gapind));
        end
        behaveind(gapind,1:length(temp)) = temp;
      end
    else
      behaveind =  behavind2;
    end
    behaviortime = zeros(2,size(behaveind,1));
    for index=1:size(behaveind,1)
      rowfixind = behaveind(index,:);
      rowfixind(rowfixind == 0) = [];
      behaviortime(:,index) = [rowfixind(1);rowfixind(end)];
    end
    
    behaviortimeArray{ii} = behaviortime;
  end  
end

behaviortime = behaviortimeArray{1};

end

function [labelsTrim, resortInd] = resortAndTrimLabels(labels)
% A function which takes in a label vector of stimuli (specifically the
% natural social stimulus set ones) and returns a resorted and abbreviated list, 
% along with the index of that resorting

% Resort
if ~any(contains(labels, 'Reward')) && ~any(contains(labels, 'headTurn'))
  gdInd = find(contains(labels, 'GoalDir'));
  idleInd = find(contains(labels, 'Idle'));
  objInd = find(contains(labels, 'objects'));
  landInd = find(contains(labels, 'landscape'));
  socInd = find(contains(labels, 'monkey') & ~contains(labels, {'Idle', 'GoalDir'}));
  resortInd = [socInd; gdInd; idleInd; objInd; landInd];
  labels = labels(resortInd);
else
  labelsTrim = labels;
  resortInd = 1:length(labels);
  return
end

% Trim
labelsTrim = labels;

if contains(labels, '.avi')
  labelsTrim = extractBefore(labelsTrim, '.avi');
end

if any(contains(labels, 'monkey'))
  labelsTrim = strrep(labelsTrim, 'monkey', '');
end

if all(contains(labels, '_'))
  labelsTrim = strcat(extractBefore(labelsTrim, '_'), cellfun(@(x) x(end), labelsTrim));
end

end
function [labelsTrim, resortInd] = resortAndTrimLabels(labels)
% A function which takes in a label vector of stimuli (specifically the
% natural social stimulus set ones) and returns a resorted and abbreviated list, 
% along with the index of that resorting

% Resort
if ~any(contains(labels, 'Reward')) && ~any(contains(labels, 'headTurn ')) % space after headTurnto catch event aligned analyses
  
  if any(contains(labels, 'headTurn'))
    labels = eventIDMerge(labels);
  end
  
  gdInd = find(contains(labels, 'GoalDir'));
  idleInd = find(contains(labels, 'Idle'));
  objInd = find(contains(labels, 'objects'));
  landInd = find(contains(labels, 'landscape'));
  ChasInd = find(contains(labels, 'Chas') & ~contains(labels, {'Idle', 'GoalDir'}));
  FigInd = find(contains(labels, 'Fig') & ~contains(labels, {'Idle', 'GoalDir'}));
  MountInd = find(contains(labels, 'Mou') & ~contains(labels, {'Idle', 'GoalDir'}));
  GroomInd = find(contains(labels, 'Groom') & ~contains(labels, {'Idle', 'GoalDir'}));
  resortInd = [ChasInd; FigInd; MountInd; GroomInd; gdInd; idleInd; objInd; landInd];
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

if all(contains(labels, '_')) && any(contains(labels, 'monkey'))
  labelsTrim = strcat(extractBefore(labelsTrim, '_'), cellfun(@(x) x(end), labelsTrim));
else
  labelsTrim = strrep(labelsTrim, '_', ' ');
end

end
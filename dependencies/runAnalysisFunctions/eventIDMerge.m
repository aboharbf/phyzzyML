function eventIDsMerged = eventIDMerge(eventIDs)
% A function for creating a merged ID list for the animated paradigms.
% Based off of code in buildStimParamFileSocialVids_Auto.

% Used in preprocessLogFileMonkeyLogic.

% Identify the paradigm
if any(contains(eventIDs, 'headTurnCon'))
  paradigm = 'HTC';
elseif any(contains(eventIDs, 'headTurnIso'))
  paradigm = 'HTI';
elseif any(contains(eventIDs, 'Barney'))
  eventIDsMerged = eventIDs;
  return
elseif any(contains(eventIDs, 'monkey'))
  eventIDsMerged = eventIDs;
  return
end

% Identify the paradigm
stimuliNames = eventIDs;

% Identify the elements of the stimulus name which point out 'unique
% stimuli' as opposed to variants of the same stimulus.

isolatedCodes = extractBetween(stimuliNames, '_1', '.');

switch paradigm
  case 'HTI' %'headTurnIso'
    % code from the jsonGenerator used to make these stim and name them.- sprintf('headTurnIso_14%d%d%d_C%dM%dS%d', prefab_i, iter_i, turn_i, cam_i, mesh_i, skin_i);
    % prefab_i, turn_i, cam_i, mesh_i matter for unique stim, iter_i and
    % skin_i make variants of these stim.
    uniqueStimID = [2 4 7 9 12];
    
    prefabToSet = {'idle','bioMotion'};
    prefabTypes = {'_leftFull', '_leftHalf', '_core', '_rightHalf', '_rightFull'};
    cameras = {'_leftCam', '_frontCam', '_rightCam'};
    meshes = {'_Normal', '_LowRes', '_Dots'};
    timeTurn = {'_Early', '_Late'};
    
    nameElements = {prefabToSet; prefabTypes; cameras; meshes; timeTurn};
  case 'HTC' %'headTurnCon'
    % Head turn con code for naming stim sprintf('headTurnCon_15%d%dT%d_E%dC%d', intCode, intNum, turnCode, env_i, cam_i);
    % intCode, intNum, turnCode mmake unique stim, env_i, cam_i make
    % variants.
    uniqueStimID = [2, 3, 5];
    
    intArray = {'chasing', 'grooming', 'mounting', 'fighting', 'idle', 'goalDirected', 'objects', 'scene'};
    intNum = {'1', '2', '3', '4', '5', '6'};
    turnArray = {'_noTurn', '_Turn'};
    
    nameElements = {intArray; intNum; turnArray};
end

% Extract unique stimuli based on code, saving the membership to each
% label.
stimCode = cellfun(@(x) x(1, uniqueStimID), isolatedCodes, 'UniformOutput', false);
[A, ~, stim2Label] = unique(stimCode);
mergeTaskEvent = cell(length(A),1);

% Need to add 1 to the headTurn number to make it usable as an index.
if strcmp(paradigm, 'HTC')
  tmp = str2double(A) + 1;
  A = num2str(tmp);
  A = mat2cell(A, ones(size(A,1),1), size(A,2));
end

% Iterate across the codes, using the previously indicated vectors to turn
% them into names.
for ii = 1:length(A)
  wholeCode = A{ii};
  eventName = cell(size(wholeCode));
  
  for jj = 1:length(uniqueStimID)
    vect = nameElements{jj};
    vectInd = str2num(wholeCode(jj));
    if isempty(vectInd)
      switch wholeCode(jj)
        case 'E'
          vectInd = 1;
        case 'L'
          vectInd = 2;
      end
    end
    eventName(jj) = vect(vectInd);
  end
  
  % Join the elements
  mergeTaskEvent{ii} = strjoin(eventName, '');
end

% Add the mergeTaskEvent label to each stimuli's catagory labels.
eventIDsMerged = cell(size(stimCode));
for ii = 1:length(mergeTaskEvent)
  labelInd = stim2Label == ii;
  eventIDsMerged(labelInd) = deal(mergeTaskEvent(ii));
end

end
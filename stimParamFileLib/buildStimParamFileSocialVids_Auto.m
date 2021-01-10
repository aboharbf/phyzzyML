function [ ] = buildStimParamFileSocialVids_Auto()
%Creates a .m file containing cell arrays detailing the parameters of a
%particular stimulus set.

%File produced 2 Cell Arrays. The first is a ParamArray, with the name of
%the video. {{VideoID1; {label1}; {label2}...}{VideoID2; {label1};
%{label2}...}. The second is a TriggerLabels, containing all possible labels. 


%Find the folder with what you want
StimFolder = {'D:\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories', ...
    'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Stimuli and Code\headTurningStimuli\Final Videos Produced'};
stimType = '.avi';

%Find every file in this folder with the extension desired.
%Add other files
tmpStimFilesStack = [];
for ii = 1:length(StimFolder)
  tmpStimFiles = dir(fullfile(StimFolder{ii}, '**', ['*' stimType]));
  tmpStimFilesStack = [tmpStimFilesStack; tmpStimFiles];
end

% Remove duplicates
stimNames = {tmpStimFilesStack.name}';
[~, uniInd] = unique(stimNames);
stimList = {tmpStimFilesStack(uniInd).name}';

%% Load the eventTable
mergeHeadTurning = true;
mergeTurning = true;

eventDataPath = fullfile(StimFolder{1}, 'eventData.mat');
tmp = load(eventDataPath);
eventDataCell = table2cell(tmp.eventData);
eventDataStim = tmp.eventData.Properties.RowNames;
eventTitles = tmp.eventData.Properties.VariableNames;

eventPresMat = false(length(stimList), size(eventDataCell,2));
eventDataLog = cellfun(@(x) ~isempty(x.startFrame), eventDataCell);

for ii = 1:length(eventDataStim)
  stimInd = strcmp(stimList, eventDataStim{ii});
  eventPresMat(stimInd,:) = eventDataLog(ii,:);
end

if mergeHeadTurning
  % Identify events w/ any head turning, and merge those columns.
  mergeRows = logical(strcmp(eventTitles, {'headTurn_right'}) + strcmp(eventTitles, {'headTurn_left'}));
  headTurn = logical(sum(eventPresMat(:, mergeRows),2));
  eventPresMat = [headTurn, eventPresMat(:, ~mergeRows)];
  eventTitles = ['headTurn', eventTitles(~mergeRows)];
end

if mergeTurning
  mergeRows = logical(strcmp(eventTitles, {'headTurn'}) + strcmp(eventTitles, {'bodyTurn'}));
  allTurn = logical(sum(eventPresMat(:, mergeRows),2));
  eventPresMat = [allTurn, eventPresMat];
  eventTitles = ['allTurn', eventTitles];  
end

%% Turn that list into the appropriate file for phyzzy..
pictureLabels = cell(size(stimList));

for ii = 1:length(stimList)
  stim = stimList{ii};
  modWord = [];
  label = [];
  stimParts = split(stim,["_","."]);
  stimParts = stimParts(1:end-1); %Remove the tag
  code = stimParts{2};
  
  %First Label
  firstLabelArray = {'agents','objects','scramble','scene'};
  stimLabels = firstLabelArray(str2num(code(1)));
  
  % 2nd Label
  secondLabelArray = {'interaction', 'NONE', 'idle', 'headTurn_isolated', 'headTurn_context'};
  contextInd = str2double(code(2));
  if contextInd ~= 0
    stimLabels = [stimLabels, secondLabelArray(contextInd)];
  end
  
  if strcmp(code(1:3), '110')
    stimLabels = horzcat(stimLabels ,'goalDirected');
  end
  
  % Make sure to not catch controls for animated stimuli
  if (((length(stimParts) > 2) && ~strncmp(stimParts(3),'C',1)) || (length(stimParts) == 2)) && ~strcmp(stimParts{1}, 'scramble') && ~code(2) == 4
    socCatagoryArray = {'chasing', 'fighting', 'mounting', 'grooming', 'holding', 'following', 'observing', 'foraging', 'sitting'};
    stimLabels = [stimLabels, socCatagoryArray(str2num(code(3))), 'socialInteraction'];

  elseif ((length(stimParts) > 2) && strncmp(stimParts(3),'C',1))
    stimLabels = horzcat(stimLabels ,'animControl');
  end
  
  if code(2) == '4'     % Code in the name of isolated head turning paradigm stim.
    % code = sprintf('headTurnIso_14%d%d%d_C%dM%dS%d', prefab_i, iter_i, turn_i, cam_i, mesh_i, skin_i);
    codeAdd = stimParts{3};
    
    motionType = {'headIso', 'bioMotion'};
    turnDirection = {'leftFull', 'leftHalf', 'noTurn', 'rightHalf', 'rightFull'};   % Direction of Turn
    cameras = {'leftCam', 'frontCam', 'rightCam'};
    mesh = {'fullModel', 'smoothModel', 'dotModel', 'sphereModel'};                 % Mesh used.
    
    condInd = str2double(code(3));   % {'Idle' vs 'bioMotion'}
    turnDirInd = str2double(code(5));
    meshInd = str2double(codeAdd(4));
    camInd = str2double(codeAdd(2));
    
    stimLabels = horzcat(stimLabels , motionType{condInd});
    stimLabels = horzcat(stimLabels , turnDirection{turnDirInd});
    stimLabels = horzcat(stimLabels , mesh{meshInd});
    
    % Use a combination of the turn direction and camera to determine
    % either 'turnAway' or 'turnToward'.
    if turnDirInd ~= 3
      stimLabels = horzcat(stimLabels , 'headTurn');
      dynLabels = {'turnToward', 'turnAway'};
      if (camInd == 3 && turnDirInd > 3) || (camInd == 1 && turnDirInd < 3)
        turnConInd = 1;
      else
        turnConInd = 2;
      end
      stimLabels = horzcat(stimLabels , dynLabels{turnConInd});     
    else
      staticLabels = {'leftProfile', 'frontal', 'rightProfile'};
      stimLabels = horzcat(stimLabels , staticLabels{camInd});
    end
    
  end
  
  if code(2) == '5'     % Code in the name of isolated head turning context stim.
    % code = sprintf('headTurnCon_15%d%dT%d_E%dC%d', intCode, intNum, turnCode, env_i, cam_i);
    scenes = {'chasing', 'grooming', 'mating', 'fighting', 'idle', 'goalDirected', 'objects', 'scene'};    % Observed Scenes
    
    % Entirely too many cameras.
    cameraNum = {'CamRing_Chase1_1', 'CamRing_Chase1_2', 'CamRing_Chase1_3', 'CamRing_Chase1_4','CamRing_Chase2_1', 'CamRing_Chase2_2',...
      'CamRing_Chase2_3', 'CamRing_Chase2_4','CamRing_Chase3_1', 'CamRing_Chase3_2', 'CamRing_Chase3_3', 'CamRing_Chase3_4', ...
      'CamRing_Groom1','CamRing_Groom2','CamRing_Groom3','CamRing_Groom4','CamRing_Mate1','CamRing_Mate2','CamRing_Mate3','CamRing_Mate4', ...
      'CamRing_Fight1','CamRing_Fight2','CamRing_Fight3','CamRing_Fit4','CamRing_Control1', 'CamRing_Control2', 'CamRing_Control3', 'CamRing_Control4',...
      'CamRing_Scene1', 'CamRing_Scene2', 'CamRing_Scene3', 'CamRing_Scene4'};
  
    turnsMat = {'noTurn', 'Turn'};
    
    sceneInd = str2double(code(3));
    turnCode = logical(str2double(code(6)));
    
    stimLabels = horzcat(stimLabels ,scenes{sceneInd});
    if turnCode
      stimLabels = horzcat(stimLabels ,'headTurn');
    end
        
  end
  
  % Animation Processing
  if strcmp(stimParts{2}(end),'A')
    stimLabels = horzcat(stimLabels ,'animated');
    %A change needed for the naming convention on animated stimuli controls
    if length(stimParts) > 2 && strncmp(stimParts(3),'C',1)
      stimLabels = horzcat(stimLabels ,'animControl');
    elseif length(stimParts) == 2
      stimLabels = horzcat(stimLabels ,'animSocialInteraction');
    end
  end
  
  % Decomposed Stimuli processing
  if any(strcmp(stimParts, 'Dephased'))
    stimLabels = {'scramble'};
    label = 'Scrambled';
  end
  
  if (length(stimParts) > 2)
    if strncmp(stimParts(3),'Decomp',5)
      modWord = 'no';
      stimLabels = horzcat(stimLabels, 'Decomp');
      switch stimParts{3}
        case 'DecompA'
          switch stimParts{4}
            case 'B';   stimLabels = horzcat(stimLabels, 'faces','background'); label = 'Bodies';
            case 'F';   stimLabels = horzcat(stimLabels, 'bodies', 'hands','background'); label = 'Faces';
            case 'H';   stimLabels = horzcat(stimLabels, 'bodies', 'faces','background'); label = 'Hands';
            case 'FB';  stimLabels = horzcat(stimLabels, 'background'); label = 'FacesBodies';
            case 'FH';  stimLabels = horzcat(stimLabels, 'bodies','background'); label = 'FaceHands';
          end
        case 'DecompB'
          switch stimParts{4}
            case 'B';   stimLabels = horzcat(stimLabels, 'bodies','hands'); label = 'FacesScene';
            case 'F';   stimLabels = horzcat(stimLabels, 'faces'); label = 'BodiesScene';
            case 'H';   stimLabels = horzcat(stimLabels, 'hands'); label = 'FacesBodiesScene';
            case 'FB';  stimLabels = horzcat(stimLabels, 'faces', 'bodies'); label = 'Scene';
            case 'FH';  stimLabels = horzcat(stimLabels, 'faces', 'hands'); label = 'BodiesScene';
          end
      end
    elseif strncmp(stimParts(3),'C',1)
      label = strjoin(stimParts(3:end));
    end
  end
  
  % Event based labels, derived from eventData.mat in stim folder.
  if any(eventPresMat(ii, :))
    stimLabels = horzcat(stimLabels, eventTitles(eventPresMat(ii, :)));
  end
  
  % Accounts for stimuli with no Labels
  stimLabels = horzcat(stimLabels, 'allStim'); %#ok<*AGROW>
  
  % Generate final structure
  stimList{ii} = horzcat(stimList{ii}, stimLabels);
  
  % Make the picture Label
  if str2double(code(1)) == 1 % Agents
    if contextInd == 1 || contextInd == 3
      eventTag = stimParts{1};
      cutStartInd = find(isstrprop(eventTag, 'upper'));
      speciesTag = eventTag(1);
      eventTag = eventTag(cutStartInd: end);
      if strcmp(eventTag(end-2:end), 'ing')
        eventTag = eventTag(1:end-3);
        if strcmp(eventTag, 'Chas')
          eventTag = 'Chase';
        end
      end
      pictureLabels{ii} = [eventTag '_' speciesTag '_' stimParts{2}(4:end) modWord label];
    elseif contextInd == 4
      pictureLabels{ii} = [motionType{condInd} '_' turnDirection{turnDirInd} '_' cameras{camInd}];
    elseif contextInd == 5
      turns = {'noTurn', 'Turn'};
      turnInd = str2num(code(6)) + 1;
      pictureLabels{ii} = [scenes{sceneInd} code(4) '_' turns{turnInd}];
    else
      pictureLabels{ii} = [stimParts{1} '_' stimParts{2}(end) modWord label];
    end
  else
    pictureLabels{ii} = [stimParts{1} '_' stimParts{2}(end) modWord label];
  end
  
end

%% Adding merged stimuli

[stimList, pictureLabels] = mergeStimuliParam(stimList, pictureLabels);

%% Adding subEvents
% events within Phyzzy may be stimuli or subEvents within the stimuli. to
% facilitate the latter, the code below looks for an 'eventData.mat', and
% add each event within as a 'subEvent'.

eventDataPath = dir([StimFolder{1} '\**\eventData.mat']);
load(fullfile(eventDataPath(1).folder, eventDataPath(1).name), 'eventData');
subEvents2Add = eventData.Properties.VariableNames;
[stimAdd, labelAdd] = deal(cell(length(subEvents2Add),1));
for event_i = 1:length(subEvents2Add)
  event = subEvents2Add(event_i);
  stimLabels = ['subEvents', strsplit(event{1},'_')];
  stimAdd{event_i} = [event, stimLabels];
  % Generate Label
  stimParts = strsplit(event{1},'_');
  if length(stimParts) > 1
    labelAdd{event_i} = [stimParts{1}, '_' upper(stimParts{2}(1))];
  else
    labelAdd{event_i} = stimParts{1};
  end
end
stimList = [stimList; stimAdd];
pictureLabels = [pictureLabels; labelAdd];

%% Package Outputs
%pictureLabels = pictureLabels(:,1);

stimLabels = cellfun(@(x) x(2:end), stimList, 'UniformOutput', false);
catLabels = unique([stimLabels{:}]);

% categoryLabels = {'agents', 'objects', 'scramble', 'scene', 'interaction', 'nonInteraction', 'idle', 'bioMotion', 'socialInteraction','nonSocialInteraction'...
%   'goalDirected', 'chasing', 'fighting', 'mounting', 'grooming', 'holding', 'following', 'observing','animated', 'animControl',...
%   'foraging','sitting','faces','bodies','hands','background', 'subEvents', 'headTurn', 'bodyTurn', 'allTurn', 'eyeContact', 'allStim'};

categoryLabels = catLabels;

paramArray = stimList;

%% Remove replicates
eventIDs = cell(size(paramArray,1),1);
for event_i = 1:length(paramArray)
  eventIDs{event_i} = paramArray{event_i}{1}; %per the stimParamFile spec, this is the event ID
end

[~, uInd, ~] = unique(eventIDs,'stable');

paramArray = paramArray(uInd);
pictureLabels = pictureLabels(uInd);

save('StimParamFileSocialVids_Full.mat','paramArray','categoryLabels','pictureLabels')
end

  function [taskEventIDs, pictureLabels] = mergeStimuliParam(taskEventIDs, pictureLabels)
% Takes in the taskEventID string, relabels variants of a stimuli to the
% same name.

% Extract stimListName
stimListNames = cellfun(@(x) x(1), taskEventIDs);

% Identify the paradigm
paradigm = extractBefore(stimListNames, '_');

% Remove those which aren't headTurn paradigm
headTurnIsoInd = strcmp(paradigm, 'headTurnIso');
headTurnConInd = strcmp(paradigm, 'headTurnCon');

headTurnInds = [headTurnIsoInd, headTurnConInd];

for jj = 1:2
  
  % Knowing that headTurnCon and headTurnIso follow a specific code,
  % identify ones which are the same.
  stimListSub = stimListNames(headTurnInds(:,jj));
  paradigmCode = extractBefore(stimListSub{1}, '_');
  isolatedCodes = extractBetween(stimListSub, '_1', '.');
  
  switch paradigmCode
    case 'headTurnIso'
      % code from the jsonGenerator used to make these stim and name them.- sprintf('headTurnIso_14%d%d%d_C%dM%dS%d', prefab_i, iter_i, turn_i, cam_i, mesh_i, skin_i);
      % prefab_i, turn_i, cam_i, mesh_i matter for unique stim, iter_i and
      % skin_i make variants of these stim.
      stimCode = cellfun(@(x) x(1, [2 4 7 9]), isolatedCodes, 'UniformOutput', false);
      [A, labels2Add, uniqueStimGroup] = unique(stimCode);
      
      % insert these names to json file via "MonkeyPrefabs":"name". use this
      % index directly into name position 3.
      prefabToSet = {'idle','bioMotion'};
      prefabCounts = 3;
      prefabTypes = {'leftFull', 'leftHalf', 'core', 'rightHalf', 'rightFull'};
      prefabID = {'','_Red', '_Ashe', '_White'}; %The first prefab, 'Beige', doesn't have the color appended.
      
      % insert these names to json file via "CameraName" field. this index goes
      % into the name in position 4.
      cameras = {'leftCam', 'frontCam', 'rightCam'};
      
      % insert the index - 1 of the desired mesh below, and set
      % 'OverrideMeshType' to 'true'. Use the index in the name as code position
      % 5.
      % meshes = {'Normal', 'LowRes', 'Dots', 'Spheres'};
      meshes = {'Normal', 'LowRes', 'Dots'};
      
      mergeTaskEvent = cell(length(A),1);
      for ii = 1:length(A)
        wholeCode = A{ii};
        intInd = str2num(wholeCode(1));   % Extract int type.
        turnInd = str2num(wholeCode(2));  % Extract head Turn
        camInd = str2num(wholeCode(3));
        meshInd = str2num(wholeCode(4));
        mergeTaskEvent{ii} = strjoin([prefabToSet(intInd), prefabTypes(turnInd), cameras(camInd), meshes(meshInd)], '_');
      end
    case 'headTurnCon'
      % Head turn con code for naming stim sprintf('headTurnCon_15%d%dT%d_E%dC%d', intCode, intNum, turnCode, env_i, cam_i);
      % intCode, intNum, turnCode mmake unique stim, env_i, cam_i make
      % variants.
      stimCode = cellfun(@(x) x(1, [2, 3, 5]), isolatedCodes, 'UniformOutput', false);
      [A, labels2Add, uniqueStimGroup] = unique(stimCode);
      
      intArray = {'chasing', 'grooming', 'mating', 'fighting', 'idle', 'goalDirected', 'objects', 'scene'};
      turnArray = {'_noTurn', '_Turn'};
      
      mergeTaskEvent = cell(length(A),1);
      for ii = 1:length(A)
        indCode = A{ii};
        intInd = str2num(indCode(1));
        intNum = indCode(2);
        turnInd = str2num(indCode(3))+1;
        mergeTaskEvent{ii} = [intArray{intInd} intNum turnArray{turnInd}];
      end
  end
  
  % Generate vector of stimuli w/ the appropriate labels following to add
  % to the larger stimList
  stimList2Add = mergeTaskEvent;
  largerGroupInd = find(headTurnInds(:,jj));
  labels = cellfun(@(x) x(2:end), taskEventIDs(largerGroupInd(labels2Add)), 'UniformOutput', false);
  
  % Combine the stimList w/ the appropriate Labels, and add to the larger
  % list.
  stimListAddition = arrayfun(@(x) [stimList2Add(x), labels{x}], 1:length(stimList2Add), 'UniformOutput', false)';
  taskEventIDs = [taskEventIDs; stimListAddition];
  pictureLabels = [pictureLabels; stimList2Add]; % Just use normal names as labels.
  
end

end
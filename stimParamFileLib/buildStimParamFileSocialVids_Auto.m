function [ ] = buildStimParamFileSocialVids_Auto()
%Creates a .m file containing cell arrays detailing the parameters of a
%particular stimulus set.

%File produced 2 Cell Arrays. The first is a ParamArray, with the name of
%the video. {{VideoID1; {label1}; {label2}...}{VideoID2; {label1};
%{label2}...}. The second is a TriggerLabels, containing all possible labels.

% Combined catagories - which are either an '&' or '|' of other catagories.
comboCat.comboCategoryLabels = {'fullModel_headTurn', 'fullModel_noTurn', 'smoothModel_headTurn', 'smoothModel_noTurn', 'dotModel_headTurn', 'dotModel_noTurn'};
comboCat.comboCategoryMembers = {{'&', 'fullModel', 'headTurn'},  {'&', 'fullModel', 'noTurn'}, ...
                      {'&', 'smoothModel', 'headTurn'}, {'&', 'smoothModel', 'noTurn'},...
                      {'&', 'dotModel', 'headTurn'},    {'&', 'dotModel', 'noTurn'}};

%Find the folder with what you want
StimFolder = {'C:\Users\aboha\Onedrive\Lab\ESIN_Ephys_Files\Stimuli and Code\SocialCategories', ...
    'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Stimuli and Code\headTurningStimuli\Final Videos Produced',...
    'C:\Users\aboha\OneDrive\Lab\ESIN_Ephys_Files\Stimuli and Code\FamiliarFace Paradigm'};
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
emptyTable = eventDataCell(1);
replaceInd = cellfun(@(x) size(x, 2) == 8, eventDataCell);
eventDataCell(~replaceInd) = emptyTable;

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

% Variables for assigning category labels
socCatagoryArray = {'chasing', 'fighting', 'mounting', 'grooming', 'holding', 'following', 'observing', 'foraging', 'sitting'};
nonSocCatagoryArray = {'goalDirected', 'Idle', 'objects', 'landscape', 'scramble'};
for lab_i = 1:length(socCatagoryArray)
  socCatagory2Find{lab_i} = [upper(socCatagoryArray{lab_i}(1)) socCatagoryArray{lab_i}(2:end)];
end

% reference structure

for ii = 1:length(stimList)
  stim = stimList{ii};
  modWord = [];
  label = [];
  stimParts = split(stim,["_","."]);
  stimParts = stimParts(1:end-1); %Remove the tag
  code = stimParts{2};
  stimLabels = {};
  refStruct(ii).stimName = stim;
      
  if ~isempty(str2num(code(1)))
    
    % 2nd Label
    secondLabelArray = {'interaction', 'NONE', 'idle', 'headTurn_isolated', 'headTurn_context'};
    contextInd = str2double(code(2));
    if contextInd ~= 0
      stimLabels = [stimLabels, secondLabelArray(contextInd)];
    end
    
    %First Label
    firstLabelArray = {'agents','objects','scramble','scene'};
    setInd = str2double(code(1));
    %headTurn isn't labeled correctly here.
    if contextInd < 4
      stimLabels = [stimLabels, firstLabelArray(setInd)];
    end
    
    % Little labeling discrepency, where goalDirected is missing a unique
    % ID
    if strcmp(code(1:3), '110')
      stimLabels = horzcat(stimLabels ,'goalDirected');
      refStruct(ii).category = 'goalDirected';
      refStruct(ii).socInt = 0;
      refStruct(ii).paradigm = 'NaturalSocial';
    end
    
    % Adding 'social interaction' tag to Natural Stim paradigm
    socIntInd = ~cellfun('isempty', regexp(stimParts{1}, socCatagory2Find));
    if any(socIntInd) % Are any of the typical action strings in the title?
      if strcmp(stimParts{2}(end), 'A') % Check for A, signalling animation. 
        if length(stimParts) > 2 && (strncmp(stimParts(3),'C',1) || strncmp(stimParts(3),'P',1)) % Check if its an animated control
          stimLabels = horzcat(stimLabels ,'animControl');
        else
          stimLabels = horzcat(stimLabels ,'animSocialInteraction');
        end
      else
        stimLabels = [stimLabels, socCatagoryArray(socIntInd), 'socialInteraction']; % socCatagoryArray = {'chasing', 'fighting', 'mounting', 'grooming', 'holding', 'following', 'observing', 'foraging', 'sitting'};
        refStruct(ii).category = socCatagoryArray{socIntInd};
        refStruct(ii).socInt = 1;
        refStruct(ii).paradigm = 'NaturalSocial';
      end
    end
    
    % Adding 'nonSocialInteraction' tag ot Natural Stim paradigm
    if any(~cellfun('isempty', regexp(stimParts{1}, nonSocCatagoryArray))) % Are any of the typical action strings in the title?
      stimLabels = [stimLabels, 'nonSocialInteraction'];
    end
    
    if code(2) == '4'     % Code in the name of isolated head turning paradigm stim.
      refStruct(ii).paradigm = 'headTurnIso';
      
      % code = sprintf('headTurnIso_14%d%d%d_C%dM%dS%d', prefab_i, iter_i, turn_i, cam_i, mesh_i, skin_i);
      codeAdd = stimParts{3};
      
      motionType = {'headIso', 'bioMotion'};
      turnDirection = {'leftFull', 'leftHalf', 'noTurn', 'rightHalf', 'rightFull'};   % Direction of Turn
      cameras = {'leftCam', 'frontCam', 'rightCam'};
      mesh = {'fullModel', 'smoothModel', 'dotModel', 'sphereModel'};                 % Mesh used.
      
      condInd = str2double(code(3));   % {'Idle' vs 'bioMotion'}
      if condInd == 2
        % Mistake where bioMotion wasn't labeled w/ 3, instead 1.
        turnDirInd = 3;
      else
        turnDirInd = str2double(code(5));
      end
      meshInd = str2double(codeAdd(4));
      camInd = str2double(codeAdd(2));
      
      stimLabels = horzcat(stimLabels , motionType{condInd});
      stimLabels = horzcat(stimLabels , turnDirection{turnDirInd});
      stimLabels = horzcat(stimLabels , mesh{meshInd});
      
      refStruct(ii).motionType = motionType{condInd};
      refStruct(ii).turnDirection = turnDirection{turnDirInd};
      refStruct(ii).mesh = mesh{meshInd};
      
      if condInd == 1
        % Use a combination of the turn direction and camera to determine
        % either 'turnAway' or 'turnToward'.
        if turnDirInd ~= 3
          refStruct(ii).turn = 1;
          stimLabels = horzcat(stimLabels , 'headTurn');
          dynLabels = {'turnToward', 'turnAway'};
          if (camInd == 3 && turnDirInd > 3) || (camInd == 1 && turnDirInd < 3)
            turnConInd = 1;
          else
            turnConInd = 2;
          end
          stimLabels = horzcat(stimLabels , dynLabels{turnConInd});
          refStruct(ii).turnDir = dynLabels{turnConInd};
        else
          staticLabels = {'leftProfile', 'frontal', 'rightProfile'};
          stimLabels = horzcat(stimLabels , staticLabels{camInd});
          refStruct(ii).staticLabels = staticLabels{camInd};
        end
      end
      
    end
    
    if code(2) == '5'     % Code in the name of context head turning paradigm stim.
      % code = sprintf('headTurnCon_15%d%dT%d_E%dC%d', intCode, intNum, turnCode, env_i, cam_i);
      scenes = {'chasing', 'grooming', 'mounting', 'fighting', 'idle', 'goalDirected', 'objects', 'scene'};    % Observed Scenes
      
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
      refStruct(ii).category = scenes{sceneInd};
      
      % Agent label for decoding
      if sceneInd == 7 || sceneInd == 8
        refStruct(ii).agent = 0;
      else
        refStruct(ii).agent = 1;
        stimLabels = horzcat(stimLabels ,'agents');
      end

      % Social interaction label
      if sceneInd < 5
        stimLabels = horzcat(stimLabels ,'socialInteraction');
        refStruct(ii).socInt = 1;
      else
        stimLabels = horzcat(stimLabels ,'nonSocialInteraction');
        refStruct(ii).socInt = 0;
      end
      
      if turnCode
        stimLabels = horzcat(stimLabels ,'headTurn');
        refStruct(ii).turn = 1;
      else
        stimLabels = horzcat(stimLabels , 'noTurn');
        refStruct(ii).turn = 0;
      end
      
    end
    refStruct(ii).paradigm = 'headTurnCon';

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
    
  else
    refStruct(ii).paradigm = 'FamiliarFace';
    
    % Replace the Entering part w/ Movement.
    if strcmp('Entering', stimParts{2})
      stimParts{2} = 'Movement';
    end
    
    % Codes for the FamiliarFace Paradigm
    if ~isempty(str2num(stimParts{1}(1)))
      % Room number, preceeds empty
      stimParts{1} = ['R' stimParts{1}];
    else
      % Monkey name, use to ID Familiar
      unFamiliarMonkeys = {'Barney', 'Calvin', 'Hobbes', 'Luca'};
      if any(strcmp(unFamiliarMonkeys, stimParts{1}))
        stimLabels = horzcat(stimLabels , 'unFamiliar');
      else
        stimLabels = horzcat(stimLabels , 'familiar');
      end
      
      % Add the picture label to the catagory list
      stimLabels = horzcat(stimLabels ,[stimParts{1} '_' stimParts{2}]);
      
    end
    
    stimLabels = horzcat(stimLabels , stimParts{1});
    stimLabels = horzcat(stimLabels , stimParts{2});

    pictureLabels{ii} = [stimParts{1} '_' stimParts{2}];
    
    stimList{ii} = horzcat(stimList{ii}, stimLabels);
    
  end
  
end

%% Adding combined catagories
stimList = AddComboCategoryLabels(stimList, comboCat);

%% Adding merged stimuli
% Function which attaches a label to all of the stimuli variants in the
% headTurn paradigms with the 'core name'.
[stimList, refStruct, pictureLabels] = mergeStimuliParam(stimList, refStruct, pictureLabels);

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
% 
stimLabels = cellfun(@(x) x(2:end), stimList, 'UniformOutput', false);
catLabels = unique([stimLabels{:}]);

% categoryLabels = {'agents', 'objects', 'scramble', 'scene', 'interaction', 'nonInteraction', 'idle', 'bioMotion', 'socialInteraction','nonSocialInteraction'...
%   'goalDirected', 'chasing', 'fighting', 'mounting', 'grooming', 'holding', 'following', 'observing','animated', 'animControl',...
%   'foraging','sitting','faces','bodies','hands','background', 'subEvents', 'headTurn', 'bodyTurn', 'allTurn', 'eyeContact', 'allStim'};

categoryLabels = catLabels';

paramArray = stimList;

%% Remove replicates
eventIDs = cell(size(paramArray,1),1);
for event_i = 1:length(paramArray)
  eventIDs{event_i} = paramArray{event_i}{1}; %per the stimParamFile spec, this is the event ID
end

[~, uInd, ~] = unique(eventIDs,'stable');

paramArray = paramArray(uInd);
pictureLabels = pictureLabels(uInd);

save('StimParamFileSocialVids_Full.mat','paramArray','categoryLabels','pictureLabels', 'refStruct')
end

function [taskEventIDs, refStruct, pictureLabels] = mergeStimuliParam(taskEventIDs, refStruct, pictureLabels)

% Identify the paradigm
stimuliNames = cellfun(@(x) x(1), taskEventIDs);

% Identify the headTurning paradigm entries
headTurnLog = strncmp(stimuliNames, 'headTurn', 8);
headTurnIsoLog = strncmp(stimuliNames, 'headTurnIso', 11);
headTurnConLog = xor(headTurnLog, headTurnIsoLog);

headTurnInds = [headTurnIsoLog, headTurnConLog];

% Knowing that headTurnCon and headTurnIso follow a specific code,
% identify ones which are the same.

for paradigm_i = 1:2
  % Identify the elements of the stimulus name which point out 'unique
  % stimuli' as opposed to variants of the same stimulus.
  
  isolatedCodes = extractBetween(stimuliNames(headTurnInds(:,paradigm_i)), '_1', '.');
  
  switch paradigm_i
    case 1 %'headTurnIso'
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
    case 2 %'headTurnCon'
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
  if paradigm_i == 2
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
  catagoryLabels2Add = cell(size(stimCode));
  for ii = 1:length(mergeTaskEvent)
    labelInd = stim2Label == ii;
    catagoryLabels2Add(labelInd) = deal(mergeTaskEvent(ii));
  end
  
  headTurnInds2Add = find(headTurnInds(:,paradigm_i));
  
  for ii = 1:length(headTurnInds2Add)
    taskEventIDs{headTurnInds2Add(ii)} = [taskEventIDs{headTurnInds2Add(ii)}, catagoryLabels2Add(ii)];
    refStruct(headTurnInds2Add(ii)).stimNameMerge = catagoryLabels2Add{ii};
  end
  
end


end

function [stimList] = AddComboCategoryLabels(stimList, comboCat)
% Adds a combined label to each stimulus per the comboCat struct.

comboCategoryLabels = comboCat.comboCategoryLabels;
comboCategoryMembers = comboCat.comboCategoryMembers;

for ii = 1:length(comboCategoryLabels)
  catMembers = comboCategoryMembers{ii}(2:end);
  catLogic = comboCategoryMembers{ii}{1};
  comboCatMembers = false(length(stimList), length(catMembers));
  % Identify all the members
  for jj = 1:length(catMembers)
    comboCatMembers(:,jj) = cellfun(@(x) any(strcmp(x, catMembers{jj})), stimList);
  end
  
  % Identify the correct logical combo
  switch catLogic
    case '&'
      comboCatInds = arrayfun(@(rowInd) all(comboCatMembers(rowInd,:)), 1:size(comboCatMembers,1))';
    case '|'
      comboCatInds = arrayfun(@(rowInd) any(comboCatMembers(rowInd,:)), 1:size(comboCatMembers,1))';
  end
  
  % Use the vector to add the label
  comboCatInds2Add = find(comboCatInds);
  for jj = comboCatInds2Add'
    stimList{jj} = [stimList{jj}, comboCategoryLabels(ii)];
  end
  
end

end
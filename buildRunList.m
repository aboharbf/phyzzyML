function [dirList] = buildRunList(dataDir, tag, paradigms2Analyze, monkeyName)
% Constructs the runList variable for the processRunBatch function based on
% the selection of a directory and a tag, which both denotes the directory
% contents of interest but also changes output formating.

filesOfInterest = dir([dataDir filesep '**' filesep '*.nev']);
paradigmListFile = fullfile(dataDir, 'paradigmList.mat');

if ~exist(paradigmListFile, 'file') && ~isempty(paradigms2Analyze)
  MLfilesOfInterest = fullfile({filesOfInterest.folder}, {filesOfInterest.name})';
  MLfilesOfInterest = extractBefore(MLfilesOfInterest, '.');
  MLfilesOfInterest = strcat(MLfilesOfInterest, '.bhv2');
  
  % Identify the monkeyLogic output files
  for ii = 1:length(MLfilesOfInterest)
    if ~exist(MLfilesOfInterest{ii}, 'file')
      MLfilesOfInterest{ii} = [extractBefore(MLfilesOfInterest{ii}, '.bhv2'), '.mat'];
    end
  end
  
  % Identify the paradigm
  paradigmList = cell(length(MLfilesOfInterest),1);
  parfor ii = 1:length(MLfilesOfInterest)
    [A, ~, ~] = mlread(MLfilesOfInterest{ii});
    
    if isfield(A(1).VariableChanges, 'stimFolder')
      
      paradigmFolder = A(1).VariableChanges.stimFolder{end};
      if contains(paradigmFolder, 'headTurnCon')
        paradigmList{ii} = 'headTurnCon';
      elseif contains(paradigmFolder, 'headTurnIso')
        paradigmList{ii} = 'headTurnIso';
      elseif contains(paradigmFolder, 'SocVNonSoc')
        paradigmList{ii} = 'naturalSocial';
      elseif contains(paradigmFolder, 'Familiar')
        paradigmList{ii} = 'familiarFace';
      end
      
    else
      
      [~, fileName] = fileparts(MLfilesOfInterest{ii});
      switch fileName(1:4)
        case '2018'
          paradigmList{ii} = 'naturalSocial_old';
        case '2019'
          
          stimName = {A(1).TaskObject.Attribute.Name}; %Pull the string containing the stimulus name.
          stimInd = find(~strcmp(stimName, 'Fixation Point'));
          stimInRun = cell(length(A),1);
          for stim_i = 1:length(A)
            stimInRun{stim_i} = A(stim_i).TaskObject.Attribute(stimInd).Name;
          end
          
          stimInRun = extractBefore(stimInRun, '.avi');
          
          % Tag for old animated stimuli
          if any(strcmp(string(cellfun(@(x) x(end), stimInRun)), "A"))
            paradigmList{ii} = 'animatedSocial_old';
          elseif all(contains(stimInRun, 'scramble'))
            paradigmList{ii} = 'scrambles';
          elseif any(contains(stimInRun, 'monkey'))
            paradigmList{ii} = 'naturalSocial_old';
          end
          
        case '2020'
          % Before the addition of the paradigm folder
          
          stimName = {A(1).TaskObject.Attribute.Name}; %Pull the string containing the stimulus name.
          stimInd = find(~strcmp(stimName, 'Fixation Point'));
          stimInRun = cell(length(A),1);
          for stim_i = 1:length(A)
            stimInRun{stim_i} = A(stim_i).TaskObject.Attribute(stimInd).Name;
          end
          
          if any(contains(stimInRun, 'headTurnCon'))
            paradigmList{ii} = 'headTurnCon';
          elseif any(contains(stimInRun, 'headTurnIso'))
            paradigmList{ii} = 'headTurnIso';
          elseif any(contains(stimInRun, 'monkey'))
            paradigmList{ii} = 'naturalSocial';
          elseif any(contains(stimInRun, 'Barney'))
            paradigmList{ii} = 'familiarFace';
          end
          
      end
      
    end
  end
  
  save(paradigmListFile, 'paradigmList', 'filesOfInterest')
end

% If Paradigm, shift the cells
if strcmp(tag, 'paradigm')
  paradigmListStruct = load(paradigmListFile);
  
  if length(paradigmListStruct.filesOfInterest) ~= length(filesOfInterest)
    delete(paradigmListFile)
    error('Paradigm list deleted due to incompleteness, rerun')
  end
  
  paradigmList = paradigmListStruct.paradigmList;
  
  run2KeepInd = false(length(filesOfInterest),1);
  % Cycle through the paradigms
  for par_i = 1:length(paradigms2Analyze)
    run2KeepInd = run2KeepInd | strcmp(paradigmList, paradigms2Analyze{par_i});
  end
  filesOfInterest = filesOfInterest(run2KeepInd);
  
  % If Monkey name specified, only run that monkey's data
  if exist('monkeyName', 'var')
    fileNameList = {filesOfInterest.name}';
    filesOfInterest = filesOfInterest(contains(fileNameList, monkeyName));
  end
  
else
  
  %Allows for customization of list
  filesOfInterest = uiDirChoose(filesOfInterest);  
  
end

%Find the name in the filename.
[~,B,~] = fileparts(filesOfInterest(1).folder);
nameInd = regexp(B,'\D');
%Extract folder names and run names.
fileNames = {filesOfInterest.name}';
dirNames = extractBetween(fileNames,1,max(nameInd));
runNames = extractBetween(fileNames,max(nameInd)+1,'.');
tmpDirName = unique(dirNames);

%Cycle through those folders
dirList = cell(length(tmpDirName),1);
for ii = 1:length(tmpDirName)
  runInd = strcmp(dirNames, tmpDirName{ii});
  runListBlock = {tmpDirName{ii}, {runNames{runInd}}};
  dirList{ii} = runListBlock;
end

end

function filesOfInterestUpdated = uiDirChoose(filesOfInterest)
theList = {filesOfInterest.name}';
listCheckFig = figure();
listCheckFig.Position(3) = 200;

listCheck = uicontrol(listCheckFig,'style', 'listbox','String',theList,'Value',1);
listCheck.Units = 'Normalized';
listCheck.Position = [.055 .05 0.9 0.8];
listCheck.Max = length(theList);

selectButton = uicontrol(listCheckFig,'style','pushbutton','String','Select','Units','Normalized','Position',[.055 .9 0.3 0.1]);
selectButton.Callback = @(src,event)uiresume();
%Wait for the select button to be hit
uiwait(listCheckFig, 120)


filesOfInterestUpdated = filesOfInterest(listCheck.Value);
close(listCheckFig)


% fprintf('Using tmp files, hardcoded \n')
% load('filesOfInterestUpdatedTmp')

end
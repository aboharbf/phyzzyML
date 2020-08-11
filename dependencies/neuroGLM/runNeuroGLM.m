function  neuroGLMStruct = runNeuroGLM(spikesByEvent, lfpByEvent, taskData, trialIDsByEvent, catIndStruct, eyeDataStruct, eyeBehStatsByStim, eyeInByEvent, params)
% Function which implements neuroGLM package from jpillow github, Park
% 2014.
% Inputs spikesByEvent{event}{channel}{unit}(trial).times
% taskData struct with 
% .taskEventIDs - trial*1 cell array with stim name per trial.
% .taskEventStartTimes - trial*1 double array with trial starts.
% eyeDataStruct struct with
% .pupilByStim - stim*1 array,

disp('Running neuroGLM...');

% add the toolbox to path
addpath(genpath(params.neuroGLMPath))

%% Initialize neuroGLM expt, add variables
params.unitOfTime = 'ms';
params.binSize = 1;
params.uniqueID = [];
params.expParam = 'Random Note I want to store';

expt = buildGLM.initExperiment(params.unitOfTime, params.binSize, params.uniqueID, 'Random Data I want to store');

% Add all the expected variables

expt = buildGLM.registerTiming(expt, 'stimOn', 'Stimulus Onset');     % events that happen 0 or more times per trial (sparse)
expt = buildGLM.registerTiming(expt, 'stimOff', 'Stimulus Off');     % events that happen 0 or more times per trial (sparse)
expt = buildGLM.registerTiming(expt, 'reward', 'Reward');             % events that happen 0 or more times per trial (sparse)
expt = buildGLM.registerTiming(expt, 'sacc',   'Saccades');           % events that happen 0 or more times per trial (sparse
expt = buildGLM.registerTiming(expt, 'blink',   'Blinks');
expt = buildGLM.registerTiming(expt, 'fix',   'Fixations');
for event_i = 1:length(taskData.eventData.events)
  expt = buildGLM.registerTiming(expt, taskData.eventData.events{event_i}, taskData.eventData.events{event_i}); % events that happen 0 or more times per trial (sparse)
end

% Values: EventID - catIndStruct.eventIDs
% Values: catagoryLabel - catIndStruct.catagoryList and eventIDs. (social,
% Agent Containing)
% Values: Failed or succeed trial, possibly.
expt = buildGLM.registerValue(expt, 'stimID', 'stimulus ID'); % information on the trial, but not associated with time
expt = buildGLM.registerValue(expt, 'social', 'Social stim'); % information on the trial, but not associated with time
expt = buildGLM.registerValue(expt, 'agent', 'Agent containing stim'); % information on the trial, but not associated with time
% expt = buildGLM.registerValue(expt, 'catagory', 'Catagory stim'); % information on the trial, but not associated with time

% Continous: eyeX, eyeY, eye Velocity, eyeD, LFP
expt = buildGLM.registerContinuous(expt, 'eyePos', 'Eye Position', 2); % 2 dimensional observation
% expt = buildGLM.registerContinuous(expt, 'eyeVel', 'Eye Velocity', 1); % 1 dimensional observation
expt = buildGLM.registerContinuous(expt, 'pupilD', 'Pupil diameter', 1); % 1 dimensional observation
expt = buildGLM.registerContinuous(expt, 'LFP', 'Local Field Potential', 1); % 1 dimensional observation

% Loop through every unit and add each spike trian to the experiment
unitNames = initNestedCellArray(spikesByEvent{1}, 'cell');
for chan_i = 1:length(unitNames)
  unitCount = length(unitNames{chan_i});
  for unit_i = 1:unitCount
    % Add each unit w/ the appropriate name to the trial
    stringName = convertUnitToName(unit_i, unitCount, 1);
    [varName, unitNames{chan_i}{unit_i}] = deal(sprintf('Ch%d%s', chan_i, stringName));
    % General name for unit
    expt = buildGLM.registerSpikeTrain(expt, varName, varName');
  end
end
expt = buildGLM.registerSpikeTrain(expt, 'spikeTrain', 'Our Neuron'); % Spike train
% expt = buildGLM.registerSpikeTrain(expt, 'spTrain2', 'Neighbor Neuron');

% Cycle through each of the variables and retrieve their structures.
trialDur = size(eyeInByEvent{1}, 3);

%% Add Trials 
% A template trial, with all vars initialized.
trialTmp = struct(); 
vars2Add = fields(expt.type);
for var_i = 1:length(vars2Add)
  trialTmp.(vars2Add{var_i}) = [];
end

trialTmp.expt = expt;
trialTmp.duration = trialDur;
trialTmp.stimOn = params.psthPre;
trialTmp.stimOff = params.psthPre + params.psthImDur;

for stim_i = 1:length(catIndStruct.eventIDs)
  
  trialTmp2 = trialTmp; % An event specific template, for values which don't change between presentations of an event.
  
  % 1. Value variables
  trialTmp2.stimID = catIndStruct.eventIDs{stim_i};
  trialTmp2.social = catIndStruct.catIndMat(stim_i, strcmp(catIndStruct.categoryList, 'socialInteraction'));
  trialTmp2.agent = catIndStruct.catIndMat(stim_i, strcmp(catIndStruct.categoryList, 'agents'));
  
  % Timing variables that don't change per stim
  for event_i = 1:length(taskData.eventData.events)
    specEventData = taskData.eventData.eventDataStarts{catIndStruct.eventIDs{stim_i}, taskData.eventData.events{event_i}}{1};
    if ~isempty(specEventData)
      trialTmp2.(taskData.eventData.events{event_i}) = specEventData';
    end
  end
  
  eventEyeBeh = eyeBehStatsByStim{stim_i};
%   trialTmp2.agent = catIndStruct.catIndMat(event_i, strcmp(catIndStruct.categoryList, 'agents'));
  % 'catagory' - figure out some sort of parsimonious way to get to
  % 'chasing', 'fighting' etc.
  
  for trial_i = 1:size(eyeInByEvent{stim_i},2)
    
    trial = trialTmp2; % Copy the template
    trialEyeBeh = eventEyeBeh{trial_i};
    
    % 2. spikeTrains
    for chan_i = 1:length(spikesByEvent{stim_i})
      for unit_i = 1:length(spikesByEvent{stim_i}{chan_i})
        % Add each unit w/ the appropriate name to the trial

        % General name for unit
        spkTrain = spikesByEvent{stim_i}{chan_i}{unit_i}(trial_i).times' + params.psthPre;
        spkTrain = spkTrain(spkTrain > 0); % since spikeTimes come from lfpAlignParams, which is pstPre + buffer.
        trial.(unitNames{chan_i}{unit_i}) = spkTrain(spkTrain < trial.duration); % since spikeTimes come from lfpAlignParams, which is pstPre + buffer.
      end
    end
    
    % 3. Timing Variables
    % May cause problems if nan.
    rewardTimeTrial = taskData.rewardTimePerTrial(trialIDsByEvent{stim_i}(trial_i));
    if ~isnan(rewardTimeTrial)
      trial.reward = params.psthPre + rewardTimeTrial;
    end
    
    if ~isempty(trialEyeBeh.saccadetimes)
      trial.sacc = trialEyeBeh.saccadetimes(1,:) + params.psthPre;
    end
    
    if ~isempty(trialEyeBeh.fixationtimes)
      trial.fix = trialEyeBeh.fixationtimes(1,:) + params.psthPre;
    end
    
    if ~isempty(trialEyeBeh.blinktimes)
      trial.blink = trialEyeBeh.blinktimes(1,:) + params.psthPre;
    else
      
    end
    % trialStruct.events = eventData file may need to be put here, or output
    % from subEventAnalysis
    
    % 4. Continous variables
    trial.eyePos = squeeze(eyeInByEvent{stim_i}(:,trial_i,:))';
    % 'eyeVel', 
    trial.pupilD = eyeDataStruct.pupilByStim{stim_i}(trial_i,:)';
    trial.LFP = squeeze(lfpByEvent{stim_i}(:, :, trial_i, params.lfpBuffer:end-params.lfpBuffer));
    
    % Add to expt
    expt = buildGLM.addTrial(expt, trial, trialIDsByEvent{stim_i}(trial_i));

  end  
end

%% Build Design Matrix 
% which specifies how to generate the design matrix
% Each covariate to include in the model and analysis is specified.
dspec = buildGLM.initDesignSpec(expt);
binfun = expt.binfun;

% covariates to include
params.covar.unitHist = 0;
params.covar.stimOnOff = 0;
params.covar.saccade = 0;
params.covar.blink = 0;
params.covar.reward = 1;
params.covar.eyePos = 0;
params.covar.LFP = 0;
params.covar.subEvents = 1;

% Add filters for each variable.
if params.covar.unitHist
for chan_i = 1:length(unitNames)
  for unit_i = [1,3] %1:length(unitNames{chan_i})-1 % Last unit is MUA
    unitLab = unitNames{chan_i}{unit_i};
%     if length([expt.trial.(unitLab)]) > 50
      dspec = buildGLM.addCovariateSpiketrain(dspec, sprintf('%sCov', unitLab), unitLab, sprintf('%s history', unitLab));
%     end
  end
end
end

% boxcar for stimulus presentation period
if params.covar.stimOnOff
  dspec = buildGLM.addCovariateBoxcar(dspec, 'stimPres', 'stimOn', 'stimOff', 'stimuls presentation period');
end

% Smooth basis functions for saccades. Offset for possible pre-saccade activity.
if params.covar.saccade
  bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 400, 8, binfun);
  offset = -200;
  dspec = buildGLM.addCovariateTiming(dspec, 'sacc', 'sacc', [], bs, offset);
end
% Smooth basis functions for blinks.
if params.covar.blink
  bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 400, 8, binfun);
  offset = -200;
  dspec = buildGLM.addCovariateTiming(dspec, 'blink', 'blink', [], bs, offset);
end
% Reward activity
if params.covar.reward
  bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 200, 8, binfun);
  dspec = buildGLM.addCovariateTiming(dspec, 'reward', [], [], bs);
end
% Pupil and eye position, Raw
if params.covar.eyePos
  bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 40, 4, binfun);
  dspec = buildGLM.addCovariateRaw(dspec, 'eyePos', 'Eye position effect', bs);
end
if params.covar.LFP
  dspec = buildGLM.addCovariateRaw(dspec, 'LFP', 'LFP effect');
end

% dspec = buildGLM.addCovariateRaw(dspec, 'pupilD', 'Pupil dilation effect');

% subEvents, smooth.
if params.covar.subEvents
  bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 300, 8, binfun);
  for event_i = 1:length(taskData.eventData.events)
    eventLab = taskData.eventData.events{event_i};
    if length([expt.trial.(eventLab)]) > 65
      dspec = buildGLM.addCovariateTiming(dspec, eventLab, eventLab, [], bs);
    end
  end
end

%% Compile the Design matrix

socialTrials = find(~[dspec.expt.trial.social]);
% trialsToRegress = 1:length(dspec.expt.trial);
trialsToRegress = socialTrials;

dm = buildGLM.compileSparseDesignMatrix(dspec, trialsToRegress);

nanCount = full(sum(sum(isnan(dm.X))));
assert(nanCount == 0, sprintf('NaNs present - %d ', nanCount))

% Remove constant columns
dm = buildGLM.removeConstantCols(dm);
dm = buildGLM.zscoreDesignMatrix(dm);

% Visualize Design matrix for first 5 trials
if params.visualizeDesignMat
  trialIndices = randi(length(expt.trial),1,5);
  endTrialIndices = cumsum(binfun([expt.trial(trialIndices).duration]));
  X = dm.X(1:endTrialIndices(end),:);

  % normalize all the values to 1
  max_val = max(abs(X), [], 1);
  max_val(isnan(max_val)) = 1;
  X = bsxfun(@times, X, 1 ./ max_val);

  % Draw image
  figure(2); 
  clf; 
  imagesc(X);
  title('Design Matrix')
end

useGLMFit = 0;
dm = buildGLM.addBiasColumn(dm);

figure()
imagesc(dm.X' * dm.X)
condest(dm.X' * dm.X)


%% Perform Regression
% Time for the magic
for chan_i = 1:length(unitNames)
  for unit_i = 2%1:length(unitNames{chan_i})-1 % Last unit is MUA

    y = buildGLM.getBinnedSpikeTrain(expt, unitNames{chan_i}{unit_i}, trialsToRegress);
    
    if useGLMFit
      % Seems to generally not work
      [w, ~, ~] = glmfit(dm.X, y, 'poisson', 'link', 'log', 'Constant', 'off');
    else
      w = dm.X' * dm.X \ dm.X' * y;
    end
    assert(~full(any(isnan(w))), 'NaN city')

    ws = buildGLM.combineWeights(dm, w);
    
    % Plot the loadings onto the basis functions
    figure()
    title(sprintf('Weights Matrix for %s (%d spikes)', unitNames{chan_i}{unit_i}, full(sum(y))));
    hold on
    [cov2plot, covLegend] = deal(fields(ws));
    for cov_i = 1:length(cov2plot)
      covData = full(ws.(cov2plot{cov_i}).data);
      
      % If it a short cov, add the number to the label
      if length(covData) < 10
        covLegend{cov_i} = sprintf('%s (%d)', covLegend{cov_i}, covData);
      end
      
      % Plot the data
      plot(covData)
      
    end
    legend(covLegend)

    % Quick predict
    figure()
    subplot(3,1,1)
    plot(y)
    title('spikeTrain')
    
    subplot(3,1,2)
    plot(w)
    title('Weights for basis functions')
    
    subplot(3,1,3)
    pred = dm.X * w;
    plot(pred)
    title('Predicted spike rates')

  end
end

%% Use results of regression to look at prediction/fit

neuroGLMStruct = struct();

end

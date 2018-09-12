function [ ] = buildStimParamFileSocialVids()
% BuildStimParamFileSocialInt
%Creates a .m file containing cell arrays detailing the parameters of a
%particular stimulus set.

%File produced 2 Cell Arrays. The first is a ParamArray, with the name of
%the video. {{VideoID1; {label1}; {label2}...}{VideoID2; {label1};
%{label2}...}. The second is a TriggerLabels, containing all possible labels. 

% For experiment run 2018.06.27

categoryLabels = {'objects', 'twoMonkeys', 'twoHumans', 'interaction', 'socialInteraction', 'nonInteraction', ...
    'fighting', 'mounting', 'grooming', 'chasing','goalDirected','scene', 'scramble', 'walking', 'highMotion', 'lowMotion'};

paramArray = {...
{'Cut_obj_interact1_1.avi'; 'objects'; 'interaction'}; ...
{'Cut_obj_interact1_2.avi'; 'objects'; 'interaction'}; ...
{'Cut_Order1Interaction_1.avi'; 'twoMonkeys'; 'socialInteraction'; 'chasing';'highMotion'}; ...
{'Cut_Order1Interaction_2.avi'; 'twoMonkeys'; 'socialInteraction'; 'fighting';'highMotion'}; ...
{'Cut_Order1Interaction_3.avi'; 'twoMonkeys'; 'socialInteraction'; 'grooming';'lowMotion'}; ...
{'Cut_Order1Interaction_4.avi'; 'twoMonkeys'; 'socialInteraction'; 'grooming'}; ...
{'Cut_Order1Interaction_5.avi'; 'twoMonkeys'; 'socialInteraction'; 'mounting';'lowMotion'}; ...
{'Cut_Order2Interaction_1.avi'; 'twoMonkeys'; 'socialInteraction'; 'fighting';'highMotion'}; ...
{'Cut_Order2Interaction_2.avi'; 'twoMonkeys'; 'socialInteraction'; 'grooming'}; ...
{'Cut_Order2Interaction_3.avi'; 'twoMonkeys'; 'socialInteraction'; 'mounting'}; ...
{'Cut_Order2Interaction_4.avi'; 'twoMonkeys'; 'socialInteraction'; 'grooming'}; ...
{'Cut_Order2Interaction_5.avi'; 'twoMonkeys'; 'socialInteraction'; 'mounting';'lowMotion'}; ...
{'Cut_Order2LandscapesFull_1.avi'; 'scene'}; ...
{'Cut_Order2LandscapesFull_2.avi'; 'scene'}; ...
{'Order1GoalLeft_1_Order1GoalRight_1.avi'; 'twoMonkeys';'nonInteraction'; 'goalDirected'}; ...
{'Order1GoalLeft_2_Order1GoalRight_2.avi'; 'twoMonkeys';'nonInteraction'; 'goalDirected'}; ...
{'Order1GoalLeft_3_Order1GoalRight_3.avi'; 'twoMonkeys';'nonInteraction'; 'goalDirected'}; ...
{'Order1GoalLeft_4_Order1GoalRight_4.avi'; 'twoMonkeys';'nonInteraction'; 'goalDirected'}; ...
{'Order1GoalLeft_5_Order1GoalRight_5.avi'; 'twoMonkeys';'nonInteraction'; 'goalDirected'}; ...
{'Order1MonkLeft_1_Order1MonkRight_1.avi'; 'twoMonkeys';'nonInteraction'; 'idle'}; ...
{'Order1MonkLeft_2_Order1MonkRight_2.avi'; 'twoMonkeys';'nonInteraction'; 'idle'}; ...
{'Order1MonkLeft_3_Order1MonkRight_3.avi'; 'twoMonkeys';'nonInteraction'; 'idle'}; ...
{'Order1MonkLeft_4_Order1MonkRight_4.avi'; 'twoMonkeys';'nonInteraction'; 'idle'}; ...
{'Order1MonkLeft_5_Order1MonkRight_5.avi'; 'twoMonkeys';'nonInteraction'; 'idle'}; ...

{'Dephased_Order1Interaction_1.avi'; 'scramble'; 'highMotion'}; ...
{'Dephased_Order1Interaction_2.avi'; 'scramble'; 'highMotion'}; ...
{'Dephased_Order1Interaction_3.avi'; 'scramble'; 'lowMotion'}; ...
{'Dephased_Order1Interaction_4.avi'; 'scramble'; 'lowMotion'}; ...

{'Dephased_Order2Interaction_1.avi'; 'scramble'; 'highMotion'}; ...
{'Dephased_Order2Interaction_2.avi'; 'scramble'; 'highMotion'}; ...
{'Dephased_Order2Interaction_3.avi'; 'scramble'; 'highMotion'}; ...
{'Dephased_Order2Interaction_4.avi'; 'scramble'; 'highMotion'}; ...
{'Dephased_Order2Interaction_5.avi'; 'scramble'; 'highMotion'}; ...

{'Cut_human_interaction_1.avi'; 'twoHumans'; 'socialInteraction'; 'chasing'}; ...
{'Cut_human_interaction_2.avi'; 'twoHumans'; 'socialInteraction'; 'walking'}; ...
{'Cut_human_interaction_3.avi'; 'twoHumans'; 'socialInteraction'; 'grooming'}; ...
{'Cut_human_interaction_4.avi'; 'twoHumans'; 'socialInteraction'; 'fighting'}; ...
{'Cut_human_interaction_5.avi'; 'twoHumans'; 'socialInteraction'; 'mounting'}; ...

{'human_alone1_1_2_human_alone2_1_2.avi';'twoHumans';'nonInteraction'; 'idle'}; ...
{'human_alone1_1_3_human_alone2_1_3.avi';'twoHumans';'nonInteraction'; 'idle'}; ...
{'human_goal1_1_4_human_goal2_1_4.avi';'twoHumans';'nonInteraction'; 'goalDirected'}; ...
{'human_goal1_1_5_human_goal2_1_5.avi';'twoHumans';'nonInteraction';'goalDirected'}; ...

{'Cut_Order1Chasing1.avi'; 'twoMonkeys'; 'socialInteraction'; 'chasing';'highMotion'}; ...
{'Cut_Order1Chasing2.avi'; 'twoMonkeys'; 'socialInteraction'; 'chasing';'highMotion'}; ...
{'Cut_Order1Chasing3.avi'; 'twoMonkeys'; 'socialInteraction'; 'chasing';'highMotion'}; ...
{'Cut_Order1Chasing4.avi'; 'twoMonkeys'; 'socialInteraction'; 'chasing';'highMotion'}; ...
{'Cut_Order1Chasing5.avi'; 'twoMonkeys'; 'socialInteraction'; 'chasing';'highMotion'}; ...

{'Cut_Order1Fighting1.avi'; 'twoMonkeys'; 'socialInteraction'; 'fighting';'highMotion'}; ...
{'Cut_Order1Fighting2.avi'; 'twoMonkeys'; 'socialInteraction'; 'fighting';'highMotion'}; ...
{'Cut_Order1Fighting3.avi'; 'twoMonkeys'; 'socialInteraction'; 'fighting';'highMotion'}; ...
{'Cut_Order1Fighting4.avi'; 'twoMonkeys'; 'socialInteraction'; 'fighting';'highMotion'}; ...
{'Cut_Order1Fighting5.avi'; 'twoMonkeys'; 'socialInteraction'; 'fighting';'highMotion'}; ...

{'Cut_Order1Mounting1.avi'; 'twoMonkeys'; 'socialInteraction'; 'mounting';'highMotion'}; ...
{'Cut_Order1Mounting2.avi'; 'twoMonkeys'; 'socialInteraction'; 'mounting';'highMotion'}; ...
{'Cut_Order1Mounting3.avi'; 'twoMonkeys'; 'socialInteraction'; 'mounting';'highMotion'}; ...
{'Cut_Order1Mounting4.avi'; 'twoMonkeys'; 'socialInteraction'; 'mounting';'highMotion'}; ...
{'Cut_Order1Mounting5.avi'; 'twoMonkeys'; 'socialInteraction'; 'mounting';'highMotion'}; ...

{'Dephased_Order1Interactionchasing.avi'; 'scramble'; 'highMotion'}; ...
{'Dephased_Order1Interactionfighting.avi'; 'scramble'; 'highMotion'}; ...
};

pictureLabels = {...
    'Objects1'; ...
    'Objects2'; ...
    
    'monkeyChasing1'; ...
    'monkeyFighting1'; ...
    'monkeyGrooming1'; ...
    'monkeyGrooming2'; ...
    'monkeyMounting1'; ...
    'monkeyFighting2'; ...
    'monkeyGrooming3'; ...
    'monkeyMounting2'; ...
    'monkeyGrooming4'; ...
    'monkeyMounting3'; ...
    'monkeyLandscape1'; ...
    'monkeyLandscape2'; ...
    
    'monkeyGoalDirected1'; ...
    'monkeyGoalDirected2'; ...
    'monkeyGoalDirected3'; ...
    'monkeyGoalDirected4'; ...
    'monkeyGoalDirected5'; ... 
    
    'monkeyIdle1'; ...
    'monkeyIdle2'; ...
    'monkeyIdle3'; ...
    'monkeyIdle4'; ...
    'monkeyIdle5'; ...
    
    'scramble1';...
    'scramble2';...
    'scramble3';...
    'scramble4';...
    
    'scramble5';...
    'scramble6';...
    'scramble7';...
    'scramble8';...    
    'scramble9';...   
    
    'humanChasing1'; ...
    'humanWalking1'; ...
    'humanGrooming1'; ...
    'humanFighting1'; ...
    'humanMating1'; ...
    
    'humanIdle1'; ...
    'humanIdle2'; ...
    'humanGoalDirected1'; ...
    'humanGoalDirected2'; ...  
    
    'monkeyChasing2'; ...
    'monkeyChasing3'; ...
    'monkeyChasing1'; ...
    'monkeyChasing4'; ...
    'monkeyChasing5'; ...
    
    'monkeyFighting2'; ...
    'monkeyFighting3'; ...
    'monkeyFighting1'; ...
    'monkeyFighting4'; ...
    'monkeyFighting5'; ...

    'monkeyMounting1'; ...
    'monkeyMounting2'; ...
    'monkeyMounting3'; ...
    'monkeyMounting4'; ...
    'monkeyMounting5'; ...
    
    'scramble10';...
    'scramble11';...
    };

save('StimParamFileSocialVids.mat','paramArray','categoryLabels','pictureLabels')



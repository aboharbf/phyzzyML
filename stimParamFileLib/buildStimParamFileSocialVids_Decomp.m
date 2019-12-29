function [ ] = buildStimParamFileSocialVids_Decomp()
% BuildStimParamFileSocialInt
%Creates a .m file containing cell arrays detailing the parameters of a
%particular stimulus set.

%File produced 2 Cell Arrays. The first is a ParamArray, with the name of
%the video. {{VideoID1; {label1}; {label2}...}{VideoID2; {label1};
%{label2}...}. The second is a TriggerLabels, containing all possible labels. 

% For experiment run 2018.10.10

categoryLabels = {'chasing2', 'chasing4', 'scramble'};

paramArray = {...
{'scramble_3003.avi'; 'scramble'}; ...
{'scramble_3005.avi'; 'scramble'}; ...

{'Cut_Order1Chasing2.avi'; 'chasing2'}; ...
{'Decomposed1_circle_B_Order1Chasing2.avi'; 'chasing2'}; ...
{'Decomposed1_circle_FB_Order1Chasing2.avi'; 'chasing2'}; ...
{'Decomposed1_circle_FH_Order1Chasing2.avi'; 'chasing2'}; ...
{'Decomposed1_circle_F_Order1Chasing2.avi'; 'chasing2'}; ...
{'Decomposed1_circle_H_Order1Chasing2.avi'; 'chasing2'}; ...


{'Cut_Order1Chasing4.avi'; 'chasing4'}; ...
{'Decomposed1_circle_B_Order1Chasing4.avi'; 'chasing4'}; ...
{'Decomposed1_circle_FB_Order1Chasing4.avi'; 'chasing4'}; ...
{'Decomposed1_circle_FH_Order1Chasing4.avi'; 'chasing4'}; ...
{'Decomposed1_circle_F_Order1Chasing4.avi'; 'chasing4'}; ...
{'Decomposed1_circle_H_Order1Chasing4.avi'; 'chasing4'}; ...
{'Decomposed2_circle_B_Order1Chasing4.avi'; 'chasing4'}; ...
{'Decomposed2_circle_FB_Order1Chasing4.avi'; 'chasing4'}; ...
{'Decomposed2_circle_FH_Order1Chasing4.avi'; 'chasing4'}; ...
{'Decomposed2_circle_F_Order1Chasing4.avi'; 'chasing4'}; ...
{'Decomposed2_circle_H_Order1Chasing4.avi'; 'chasing4'}; ...
{'Dephased_Order1Chasing4.avi'; 'scramble'}; ...
};

pictureLabels = {...
'scramble3'; ...
'scramble5'; ...

'Chasing2'; ...
'Chasing2_Decom1_B'; ...
'Chasing2_Decom1_FB'; ...
'Chasing2_Decom1_FH'; ...
'Chasing2_Decom1_F'; ...
'Chasing2_Decom1_H'; ...

'Chasing4'; ...
'Chasing4_Decom1_B'; ...
'Chasing4_Decom1_FB'; ...
'Chasing4_Decom1_FH'; ...
'Chasing4_Decom1_F'; ...
'Chasing4_Decom1_H'; ...
'Chasing4_Decom2_B'; ...
'Chasing4_Decom2_FB'; ...
'Chasing4_Decom2_FH'; ...
'Chasing4_Decom2_F'; ...
'Chasing4_Decom2_H'; ...
'Chasing4Scramble'; ...

'scramble5'; ...

};


save('StimParamFileSocialVids_Decomp.mat','paramArray','categoryLabels','pictureLabels')
%% 1.1. Initialize Matlab

clear all
close all
clc

%% 1.2. You need to give the correct name of the alloy
alloy = '106_RT_R_x100_96';
save_folder = 'Results';
savename = [save_folder filesep alloy];

%% 2.1. Reconstruct grain : test

% run this part to find the correct angle and grainSize treshold to apply
% compare "Raw IPF_X" and "Grain Boundary" figures
% once happy with refining parameter, run 2.2 "reconstruct Grain"
setMTEXpref('xAxisDirection','east'); % orientation of the map
setMTEXpref('zAxisDirection','intoPlane');

angle_test = 3; % test value, play with it
grainSize_test = 0; % test value, play with it

fname = [mtexDataPath filesep 'EBSD' filesep alloy '.cpr'];
ebsd_test = EBSD.load(fname,'convertEuler2SpatialReferenceFrame','setting 4');
ebsd_test = ebsd_test('Cr'); % we only work on the FCC indexed phase named Cr

[grains_test,ebsd_test.grainId,ebsd_test.mis2mean] = calcGrains(ebsd_test,'angle',angle_test*degree);
ebsd_test(grains_test(grains_test.grainSize<=grainSize_test)) = [];
[grains_test,ebsd_test.grainId,ebsd_test.mis2mean] = calcGrains(ebsd_test,'angle',angle_test*degree);
oM_test = ipfHSVKey(ebsd_test); oM_test.inversePoleFigureDirection = xvector;
color_test = oM_test.orientation2color(ebsd_test.orientations);

figure; plot(ebsd_test,color_test);mtexTitle('raw IPFX');
figure; plot(grains_test);mtexTitle('Current Grain Boundary');

%% 2.2. Reconstruct grain : Main

% change the inputs with the chosen values after the test in 2.1.
% 3/20/3 works well for x4000 maps 3/200/3 for x500 maps
angle_thresh = 10; % max misorientation before its a GB
grain_thresh = 5; % min pixel count to be a grain
smooth_thresh = 3; % amount of GB smoothing

B_Grain_reconstruction;

save([savename '_ebsd.mat'],'ebsd') 
save([savename '_grains.mat'],'grains')

%% 2.2. Or load if already exist
ebsd = load([savename '_ebsd.mat']).ebsd;
grains = load([savename '_grains.mat']).grains;


%% 3.1. Test Grain Boundaries

% Step 0 : run twin_fraction_test one first time to see how the segmentation went
twin_thresh = 5; % allowed deviation from 60° to be considered as TwinBoundary
min_twin_frac = 0.01; % min fraction of detected twin boundary to consider the whole GB as TB

E_TwinFraction_test_1; 

% Step 1: identify uncorrect twins in red (too big or not straigth lines), use "Twin Identification maps that just popped 
% to select grains, view the IDs in the command window and enter the grain pairs in the following list

%% 3.2. Correct Twin fraction

% Step 2: corrected twins: a/ remove boundaries that are neither GB or TB in "Removed_Boundaries" (from maps of test_1)
% then run once test 2 and identify TB_to_GB and GB_to_TB with correct IDs
% b/ Change TB to GB in TB_to_GB; c/ Change GB to TB in GB_to_TB
% d/ inverse twin/parent regions for grains where twin area is bigger than parents in Inv_grains
% please start from top left to right bottom to keep track of grain with increasing ID

savename = [save_folder filesep alloy '_CorrectedBoundaries.mat'];

Removed_Boundaries = [ % use the purple grey map IDs
   
    ];


TB_to_GB = [ % use the red grey map IDs
    
    
    ];

GB_to_TB = [
    767,758;767,732;767,765;
    1161,1194;1171,1205;
    1228,1229;
    1167,1141;1141,1116;
    898,950;
    589,584;
    977,928;977,971;
    689,718;
    1307,1343;1343,1329;
    1480,1573;
    ];

TB_to_no = [
    
    ];

Inv_grains = [
    
    ];

params.angle_thresh = angle_thresh; params.grain_thresh = grain_thresh;
params.smooth_thresh = smooth_thresh; params.twin_thresh = twin_thresh;
params.min_twin_frac = min_twin_frac; 
save(savename, "Removed_Boundaries","TB_to_GB","GB_to_TB","TB_to_no","Inv_grains","params");



E_TwinFraction_test_2;

% look at the distribution graphs and open the excel, check the merged grains with largest twin percentages
% track their id into targetId and enter it in the variable below to see if the grain has been done correctly

%% 3.2. or if done previously

savename = [save_folder filesep alloy '_CorrectedBoundaries.mat'];

% if grains performed previously you can load the data of the commented lines below instead of filling the arrays
PrevDa=load(savename); Removed_Boundaries=PrevDa.Removed_Boundaries; 
TB_to_GB=PrevDa.TB_to_GB; GB_to_TB=PrevDa.GB_to_TB; 
Inv_grains=PrevDa.Inv_grains;
params=PrevDa.params;

E_TwinFraction_test_2;

%% 3.3. Run TwinFraction

grain_dist_thresh = 0.1; % you can remove grains that are too small in the distribution

E_TwinFraction_test_3;


%% 4.0. FCC Twin relation and rotations(not required anymore)
% C_Inputs;
% D_TwinThickness;

%% 5.1. Texture
F_Texture_FCC;

%% all together

% clear all
% close all
% clc
% 
% B_Grain_reconstruction;
% C_Inputs;
% D_TwinThickness;
% E_TwinFraction_FCC;
% F_Texture_FCC;

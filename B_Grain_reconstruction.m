%% Load and save EBSD locally 

% alloy = '33RT_31_4';

save_folder = 'Results';
savename = [save_folder filesep alloy];

fname = [mtexDataPath filesep 'EBSD' filesep alloy '.cpr'];
ebsd = EBSD.load(fname,'convertEuler2SpatialReferenceFrame','setting 4');
ebsd = ebsd('Cr'); % we only work on the FCC indexed phase named Cr



%% Reconstruct grains

%angle_thresh = 10;
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',angle_thresh*degree);
%grain_thresh = 10;
ebsd(grains(grains.grainSize<=grain_thresh)) = [];
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',angle_thresh*degree);
% smooth the grains
%smooth_thresh = 3;
grains = grains.smooth(smooth_thresh);
big_grains = grains(grains.grainSize >1500);

%% Test

setMTEXpref('xAxisDirection','east'); % orientation of the map
setMTEXpref('zAxisDirection','intoPlane');

% hg_ebsd = scaleBar(ebsd, 'um');
% hg_grain = scaleBar(grain, 'um');
oM = ipfHSVKey(ebsd('Cr')); % IPF of the map
oM_grains = ipfHSVKey(grains('Cr'));
oM.inversePoleFigureDirection = xvector; % define the direction of the ipf
oM_grains.inversePoleFigureDirection = xvector;
color = oM.orientation2color(ebsd('Cr').orientations);
color_grain = oM_grains.orientation2color(grains.meanOrientation);

figure; plot(ebsd('Cr'),color);mtexTitle('IPF X');
saveFigure([savename '_raw IPF X'])

figure; plot(grains, color_grain);mtexTitle('Child Grain Reconstruction');
txt = sprintf('min angle GB: %d     grain size min: %d    smoothing iterations: %d', angle_thresh, grain_thresh, smooth_thresh);
text(big_grains, int2str(big_grains.id))
saveFigure([savename '_reconstructed IPF X'])

%%
%figure; plot(ebsd, ebsd('Cr').prop.bs), mtexTitle('BandContrast');

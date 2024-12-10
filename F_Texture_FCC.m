%% input EBSD Data
save_folder = 'Results';
% alloy = '33RT_31_4';
% fname = [mtexDataPath filesep 'EBSD' filesep alloy '.cpr'];
% ebsd = EBSD.load(fname,'convertEuler2SpatialReferenceFrame','setting 4');
% ebsd = ebsd('Cr'); % we only work on the FCC indexed phase named Cr
% 
% % cleaning grains
% [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',5*degree);
% grain_thresh = 25;
% ebsd(grains(grains.grainSize<=grain_thresh)) = [];
% [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',5*degree);
% cs = grains.CS; %get the crystal system after cleaning

%% set reference orientations

oM = ipfHSVKey(grains('Cr'));
oM.inversePoleFigureDirection = xvector;
colors = oM.orientation2color(ebsd('Cr').orientations);
figure; plot(ebsd('Cr'),colors);
figure; plot(oM)

%% PF and IPF figures plots

% Define the sample symmetry (assuming X for tensile)
oM.inversePoleFigureDirection = xvector;

psi = calcKernel(grains('Cr').meanOrientation);
% Compute the Orientation Density Function (ODF)
odf = calcDensity(ebsd('Cr').orientations,'kernel',psi);
max_odf_value = max(odf);

% Plot the IPF and Distribution with Contour Levels
figure;
plotIPDF(odf, xvector);
%mtexColorMap hsvColorMap; % Set a color map suitable for pole figures
mtexColorbar;          % Add a colorbar to visualize the density
%caxis([0 max_odf_value]);

save_name = [save_folder filesep alloy '_IPF_X_distribution'];
saveFigure(save_name)

% Plot the Pole Figures (PDF) for Specific Miller Indices
h = [Miller(1,0,0,odf.CS), Miller(1,1,0,odf.CS), Miller(1,1,1,odf.CS)];
figure;
plotPDF(odf,h);
%mtexColorMap hsvColorMap;
mtexColorbar;

save_name = [save_folder filesep alloy '_PF_X_distribution'];
saveFigure(save_name)

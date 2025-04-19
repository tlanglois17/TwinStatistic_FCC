% %% input and cleaning EBSD Data

save_folder = 'Results';
savename = [save_folder filesep alloy];

%% identify twin boundaries

cs = grains.CS;
big_grains = grains(grains.grainSize > 4000);
gB = grains.boundary;
gB_fcc = gB('Cr', 'Cr');

% Filter out wrong boundaries that are neither GB or TB
RemovedBoundary = false(size(gB_fcc.grainId, 1), 1);
for i = 1:size(Removed_Boundaries, 1)
    grain1ID = Removed_Boundaries(i, 1);
    grain2ID = Removed_Boundaries(i, 2);
    RemovedBoundary = RemovedBoundary | ...
                 ((gB_fcc.grainId(:, 1) == grain1ID & gB_fcc.grainId(:, 2) == grain2ID) | ...
                  (gB_fcc.grainId(:, 1) == grain2ID & gB_fcc.grainId(:, 2) == grain1ID));
end

gB_corr_0 = gB_fcc(RemovedBoundary);
[grains_v1, ParentId_v1] = merge(grains, gB_corr_0,'calcMeanOrientation'); % merge removes the selected boundaries gB_0
%gB_corr_1 = gB_fcc(~RemovedBoundary); % equals to TB+BG (all black)
gB_corr_1 = grains_v1.boundary; % alternative definition but with different grains ids

%figure; plot(grains_v1.boundary,'lineColor','grey');text(grains_v1, int2str(grains_v1.id),'color', 'grey');
% figure; plot(grains_v1, 'FaceColor','white'); hold on; plot(gB_corr_1,'linecolor','black'); 
% text(1,1, 'wrong boundaries removed'); hold off
% figure; plot(grains_v1, 'FaceColor','white'); hold on; plot(gB_corr_1,'linecolor','black'); plot(gB_corr_0,'linecolor','green', 'lineWidth', 2); hold off
% text(1,1, 'wrong boundaries highlighted'); hold off

%%
% Angular criterion for twinning boundaries
angle_crit = (gB_corr_1.misorientation.angle > 59*degree) & (gB_corr_1.misorientation.angle < 60*degree); 
mori = gB_corr_1.misorientation(angle_crit);
mori_mean = mean(mori, 'robust');  % Determine the mean of the cluster
round2Miller(mori_mean);  % Execute to get special orientations for FCC systems

% Twinning boundary criterion
twinning = orientation.map(...
    Miller(1, 2, 2, cs), Miller(0, 1, 0, cs), ... 
    Miller(0, 1, -1, cs, 'uvw'), Miller(-1, 0, -1, cs, 'uvw'));
isTwinning = angle(gB_corr_1.misorientation, twinning) < twin_thresh*degree;
twinBoundary_v1 = gB_corr_1(isTwinning); 

% twinBoundaries: type grainBoundary (with 1 for TB and 0 for GB)
% identify all TB inside gB_correct_1 (all GB and TB) but some are wrongs...

% figure; plot(grains_v1, 'FaceColor','white'); hold on; 
% plot(gB_corr_1,'linecolor','black');plot(twinBoundary_v1,'linecolor','red','lineWidth', 2); hold off

%%
% Filter out wrong twins that should be GB
Wrong_Twins_in_GB = false(size(gB_corr_1.grainId, 1), 1);
for i = 1:size(TB_to_GB, 1)
    grain1ID = TB_to_GB(i, 1);
    grain2ID = TB_to_GB(i, 2);
    Wrong_Twins_in_GB = Wrong_Twins_in_GB | ...
                     ((gB_corr_1.grainId(:, 1) == grain1ID & gB_corr_1.grainId(:, 2) == grain2ID) | ...
                      (gB_corr_1.grainId(:, 1) == grain2ID & gB_corr_1.grainId(:, 2) == grain1ID));
end

Wrong_Twins_from_TB = false(size(twinBoundary_v1.grainId, 1), 1);
for i = 1:size(TB_to_GB, 1)
    grain1ID = TB_to_GB(i, 1);
    grain2ID = TB_to_GB(i, 2);
    Wrong_Twins_from_TB = Wrong_Twins_from_TB | ...
                     ((twinBoundary_v1.grainId(:, 1) == grain1ID & twinBoundary_v1.grainId(:, 2) == grain2ID) | ...
                      (twinBoundary_v1.grainId(:, 1) == grain2ID & twinBoundary_v1.grainId(:, 2) == grain1ID));
end

parentBoundary_v1 = [gB_corr_1(~isTwinning), gB_corr_1(Wrong_Twins_in_GB)]; %wrt grains.id

twinBoundary_v1 = gB_corr_1(isTwinning);
twinBoundary_v2 = twinBoundary_v1(~Wrong_Twins_from_TB);

% figure; 
% plot(grains_v1, 'FaceColor','white'); hold on; 
% plot(gB_corr_1(Wrong_Twins_in_GB),'linecolor','blue','lineWidth', 2);
% plot(parentBoundary_v1,'linecolor','black','lineWidth', 2);
% plot(twinBoundary_v2,'linecolor','red','lineWidth', 1);
% text(1,1,'red lines are twins, will be englobed by parent, blacks lines are parents')
% hold off


%%
% Now add Missing_TB_pairs into twinBoundary
Missing_TB = false(size(parentBoundary_v1.grainId, 1), 1);
for i = 1:size(GB_to_TB, 1)
    grain1ID = GB_to_TB(i, 1);
    grain2ID = GB_to_TB(i, 2);
    % Find boundaries between grain1ID and grain2ID in gB_filtered
    Missing_TB = Missing_TB | ...
            ((parentBoundary_v1.grainId(:, 1) == grain1ID & parentBoundary_v1.grainId(:, 2) == grain2ID) | ...
             (parentBoundary_v1.grainId(:, 1) == grain2ID & parentBoundary_v1.grainId(:, 2) == grain1ID));
end

parentBoundary_v2 = parentBoundary_v1(~Missing_TB);
twinBoundary_v3 = [twinBoundary_v2; parentBoundary_v1(Missing_TB)];

% figure; 
% plot(grains_v1, 'FaceColor','white'); hold on; 
% plot(parentBoundary_v2,'linecolor','black','lineWidth', 2);
% plot(twinBoundary_v3,'linecolor','red','lineWidth', 2);
% text(1,1,'final corrections')
% hold off

%% Refining twins so that deformed parent grain boundaries are not considered as twins - R.Vargas Copyrights   

[listTwin,indexList] = unique(twinBoundary_v3.grainId,'rows'); % listTwin returns a sorted list without repetition of grain pairs with common TB
perimeterGrains = grains_v1(twinBoundary_v3.grainId).perimeter; % return pair of GB length associated to a pair of grains 
twinPerimeter = full(twinBoundary_v3.perimeter(unique(listTwin))); % return sorted list with total TB length surrounding a grain with TB
twinPerimeter(end:end+length(grains_v1)-length(twinPerimeter)) = 0; % add elements to twinPerimeter list for grains w/o TB to zero value 
[~,sortedTwin] = sort(grains_v1.grainSize); % list sortedTwin is rearranged from TwinPerimeter with grain sorted by ascending TB length 

for i=1:length(twinPerimeter) %assigning TB to the correct grain
    I1 = find(sortedTwin(i)==listTwin(:,1)); % find the index of left grain of the TB of range i, if TB i is not null
    I2 = find(sortedTwin(i)==listTwin(:,2)); % find the index of righ grain of the TB of range i, if TB i is not null    
    if isempty(I1) && isempty(I2)
        continue
    else
        if ~isempty(I1)
            twinPerimeter(listTwin(I1,2)) = twinPerimeter(listTwin(I1,2))-twinPerimeter(listTwin(I1,1)); % grains sorted by ascending TB: remove TB lenght of bigger grains assuming, its the parent grain
        end
        if ~isempty(I2)
            twinPerimeter(listTwin(I2,1)) = twinPerimeter(listTwin(I2,1))-twinPerimeter(listTwin(I2,2)); % same
        end
        I = true(size(listTwin,1),1);
        I(I1) = false;
        I(I2) = false;
        listTwin = listTwin(I,:);
    end
end

perimeterTwin = twinPerimeter(twinBoundary_v3.grainId);
percentageTwin = max(perimeterTwin./grains_v1(twinBoundary_v3.grainId).area,[],2);
length_crit = percentageTwin > 0.01;
twinBoundary_new = twinBoundary_v3(length_crit); 

% figure; plot(gB_corr_1); hold on
% plot(twinBoundary_v1,twinBoundary_v1.misorientation.angle./degree,'linewidth',2,'displayName','\Sigma3')
% mtexColorbar('Special boundaries'); mtexTitle('Special Boundaries'); %visualise the twin
% plot(twinBoundary_new,'linecolor','r','linewidth',1,'linestyle', '-','displayName','twin boundary')
% hold off

%% Reconstruct ParentGrain structure

[mergedGrains,parentId] = merge(grains_v1, twinBoundary_new);
 
% figure; plot(twinBoundary_v3,'linecolor','r','linewidth',1,'displayName','twin boundary'); hold on
% plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-','displayName','merged grains');
% mtexTitle('Merging Parent grains');hold off

% figure; plot(grains_v1.boundary, 'lineColor', 'grey'); hold on
% text(grains_v1, int2str(grains_v1.id),'Color', 'grey'); 
% plot(mergedGrains.boundary, 'lineColor','cyan'); 
% text(mergedGrains, int2str(mergedGrains.id),'Color', 'cyan'); hold off

%% Overall Twin fraction (per surface)
gArea = grains_v1.area;
big_grains_v1 = grains_v1(grains_v1.grainSize > 500);

% loop over mergedGrains and determine children that are not twins
isTwin = true(grains_v1.length,1); % initiate the twin identification 

for i = 1:mergedGrains.length % check in every merged grain
   MG_Id = find(parentId==i); % get child ids
   [fId,center] = calcCluster(grains_v1.meanOrientation(MG_Id),'maxAngle',... % fid is an ID from 1 to max diff average orientation
       15*degree,'method','hierarchical','silent'); % cluster grains of similar orientations
   clusterArea = accumarray(fId,gArea(MG_Id)); % compute area of each cluster
   [~,fParent] = max(clusterArea); % label the grains of largest cluster as original grain
   isTwin(MG_Id(fId==fParent)) = false;
end

for i = 1:length(TB_to_no)
    % Find indices of grains associated with the current Parent ID
    isTwin(TB_to_no(i)) = ~isTwin(TB_to_no(i));
end

% visualize the result
figure; plot(grains_v1(~isTwin),'FaceColor','darkgray','displayName','not twin'); hold on
plot(grains_v1(isTwin),'FaceColor','red','displayName','twin'); hold on
plot(mergedGrains.boundary,'linecolor','k','linewidth',2,'linestyle','-','displayName','merged grains');
plot(grains_v1.boundary,'linecolor','grey');
text(big_grains_v1, int2str(big_grains_v1.id),'Color', 'yellow');
hold off;

% figure; plot(grains, color_grain); hold on
% plot(grains_v1(isTwin).boundary,'LineColor','yellow','linewidth',1,'displayName','twin'); hold on
% plot(mergedGrains.boundary,'linecolor','k','linewidth',1,'linestyle','-','displayName','merged grains');
% text(big_grains_v1, int2str(big_grains_v1.id),'Color', 'yellow');
% mtexTitle('twin identification manually'); 
% hold off;

%% manual corrected twin/parent inversion

big_mergedGrains = mergedGrains(mergedGrains.grainSize > 25);

% Loop through each parent ID in the inversion list and invert `isTwin` condition for its children
for i = 1:length(Inv_grains)
    % Find indices of grains associated with the current Parent ID
    InvGrainIndices = find(parentId == Inv_grains(i));
    
    % Invert `isTwin` condition for these grains
    isTwin(InvGrainIndices) = ~isTwin(InvGrainIndices);
end

%Visualize the updated twin identification with inverted regions
figure;
plot(grains_v1(~isTwin), 'FaceColor', 'darkgray', 'DisplayName', 'Not Twin'); hold on;
plot(grains_v1(isTwin), 'FaceColor', 'red', 'DisplayName', 'Twin'); hold on;
plot(mergedGrains.boundary, 'LineColor', 'k', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Merged Grains');
text(big_mergedGrains, int2str(big_mergedGrains.id),'Color', 'yellow');
mtexTitle('Twin Identification with Inversions');
hold off;

fig = figure('Visible', 'off');
plot(grains_v1(~isTwin), 'FaceColor', 'darkgray', 'DisplayName', 'Not Twin'); hold on;
plot(grains_v1(isTwin), 'FaceColor', 'red', 'DisplayName', 'Twin'); hold on;
plot(mergedGrains.boundary, 'LineColor', 'k', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Merged Grains');
%legend('Location', 'southeast');
savename = [save_folder filesep alloy '_TwinMap'];
saveFigure([savename])
close(fig); % Close the figure after saving to free up memory

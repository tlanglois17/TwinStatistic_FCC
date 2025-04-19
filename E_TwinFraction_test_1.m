oM = ipfHSVKey(ebsd('Iron fcc')); % IPF of the map
oM_grains = ipfHSVKey(grains('Iron fcc'));
oM.inversePoleFigureDirection = xvector; % define the direction of the ipf
oM_grains.inversePoleFigureDirection = xvector;
color = oM.orientation2color(ebsd('Iron fcc').orientations);
color_grain = oM_grains.orientation2color(grains.meanOrientation);

cs = grains.CS;
big_grains = grains(grains.grainSize >4000);
gB = grains.boundary;
gB_fcc = gB('Iron fcc','Iron fcc');
angle_crit = (gB_fcc.misorientation.angle > 59*degree) & (gB_fcc.misorientation.angle < 60*degree); %angular criterion for twinning boundaries
mori = gB_fcc.misorientation(angle_crit);

mori_mean = mean(mori,'robust'); % determine the mean of the cluster
round2Miller(mori_mean); % execute this line to get the special orientations and fill the next lines for fcc systems
twinning = orientation.map(...
    Miller(1, 2, 2, cs), Miller(0, 1, 0, cs), ... 
    Miller(0, 1, -1, cs, 'uvw'), Miller(-1, 0, -1, cs, 'uvw')); 
isTwinning = angle(gB_fcc.misorientation,twinning) < twin_thresh*degree;
twinBoundary = gB_fcc(isTwinning);

%% Refining twins so that deformed parent grain boundaries are not considered as twins - R.Vargas Copyrights   

[listTwin,indexList] = unique(twinBoundary.grainId,'rows'); % listTwin returns a sorted list without repetition of grain pairs with common TB
perimeterGrains = grains(twinBoundary.grainId).perimeter; % return pair of GB length associated to a pair of grains 
twinPerimeter = full(twinBoundary.perimeter(unique(listTwin))); % return sorted list with total TB length surrounding a grain with TB
twinPerimeter(end:end+length(grains)-length(twinPerimeter)) = 0; % add elements to twinPerimeter list for grains w/o TB to zero value 
[~,sortedTwin] = sort(grains.grainSize); % list sortedTwin is rearranged from TwinPerimeter with grain sorted by ascending TB length 

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

perimeterTwin = twinPerimeter(twinBoundary.grainId);
percentageTwin = max(perimeterTwin./grains(twinBoundary.grainId).area,[],2);
length_crit = percentageTwin > min_twin_frac;
twinBoundary_new = twinBoundary(length_crit); 

%% Reconstruct ParentGrain structure

[mergedGrains,parentId] = merge(grains,twinBoundary_new);
gArea = grains.area;

% loop over mergedGrains and determine children that are not twins
isTwin = true(grains.length,1);
for i = 1:mergedGrains.length
   childId = find(parentId==i); % get child ids
   [fId,center] = calcCluster(grains.meanOrientation(childId),'maxAngle',...
       15*degree,'method','hierarchical','silent'); % cluster grains of similar orientations
   clusterArea = accumarray(fId,gArea(childId)); % compute area of each cluster
   [~,fParent] = max(clusterArea); % label the grains of largest cluster as original grain
   isTwin(childId(fId==fParent)) = false;
end

% visualize the result
figure; plot(grains(~isTwin),'FaceColor','cyan','displayName','not twin'); hold on
plot(grains(isTwin),'FaceColor','purple','displayName','twin'); hold on
plot(mergedGrains.boundary,'linecolor','k','linewidth',2,'linestyle','-','displayName','merged grains');
plot(grains.boundary, 'lineColor','grey');
text(grains, int2str(grains.id));
mtexTitle('twin identification'); hold off

figure; plot(grains, color_grain); hold on
plot(mergedGrains.boundary,'linecolor','k','linewidth',2,'linestyle','-','displayName','merged grains');
plot(grains.boundary, 'lineColor','grey'); text(grains, int2str(grains.id));
oM.inversePoleFigureDirection = xvector; color = oM.orientation2color(ebsd.orientations); hold off
figure; plot(ebsd('Iron fcc'),color);mtexTitle('IPF X');
oM.inversePoleFigureDirection = yvector; color = oM.orientation2color(ebsd.orientations); hold off
figure; plot(ebsd('Iron fcc'),color);mtexTitle('IPF Y');
oM.inversePoleFigureDirection = zvector; color = oM.orientation2color(ebsd.orientations);
figure; plot(ebsd('Iron fcc'),color);mtexTitle('IPF Z');

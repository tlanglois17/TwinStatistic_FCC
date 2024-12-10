% %% input and cleaning EBSD Data

save_folder = 'Results';
savename = [save_folder filesep alloy];

%% identify twin boundaries

cs = grains.CS;
big_grains = grains(grains.grainSize >4000);
gB = grains.boundary;
gB_fcc = gB('Cr','Cr');

isInGrainPairs = false(size(gB_fcc.grainId, 1), 1); % will remove wrong twins
for i = 1:size(Not_TB_Pairs, 1)
    grain1ID = Not_TB_Pairs(i, 1);
    grain2ID = Not_TB_Pairs(i, 2);
    
    % Update the mask for boundaries between grain1ID and grain2ID
    isInGrainPairs = isInGrainPairs | ...
                     ((gB_fcc.grainId(:, 1) == grain1ID & gB_fcc.grainId(:, 2) == grain2ID) | ...
                      (gB_fcc.grainId(:, 1) == grain2ID & gB_fcc.grainId(:, 2) == grain1ID));
end

gB_filtered = gB_fcc(~isInGrainPairs); 
angle_crit = (gB_filtered.misorientation.angle > 59*degree) & (gB_filtered.misorientation.angle < 60*degree); %angular criterion for twinning boundaries
mori = gB_filtered.misorientation(angle_crit);
mori_mean = mean(mori,'robust'); % determine the mean of the cluster
round2Miller(mori_mean); % execute this line to get the special orientations and fill the next lines for fcc systems
twinning = orientation.map(...
    Miller(1, 2, 2, cs), Miller(0, 1, 0, cs), ... 
    Miller(0, 1, -1, cs, 'uvw'), Miller(-1, 0, -1, cs, 'uvw')); 
isTwinning = angle(gB_filtered.misorientation,twinning) < twin_thresh*degree;
twinBoundary = gB_filtered(isTwinning);

%%

cs = grains.CS;
big_grains = grains(grains.grainSize > 4000);
gB = grains.boundary;
gB_fcc = gB('Cr', 'Cr');

% Filter out wrong twins based on Not_TB_Pairs
isInGrainPairs = false(size(gB_fcc.grainId, 1), 1);
for i = 1:size(Not_TB_Pairs, 1)
    grain1ID = Not_TB_Pairs(i, 1);
    grain2ID = Not_TB_Pairs(i, 2);
    
    % Update the mask for boundaries between grain1ID and grain2ID
    isInGrainPairs = isInGrainPairs | ...
                     ((gB_fcc.grainId(:, 1) == grain1ID & gB_fcc.grainId(:, 2) == grain2ID) | ...
                      (gB_fcc.grainId(:, 1) == grain2ID & gB_fcc.grainId(:, 2) == grain1ID));
end

gB_filtered = gB_fcc(~isInGrainPairs); 

% Angular criterion for twinning boundaries
angle_crit = (gB_filtered.misorientation.angle > 59*degree) & (gB_filtered.misorientation.angle < 60*degree); 
mori = gB_filtered.misorientation(angle_crit);
mori_mean = mean(mori, 'robust');  % Determine the mean of the cluster
round2Miller(mori_mean);  % Execute to get special orientations for FCC systems

% Twinning boundary criterion
twinning = orientation.map(...
    Miller(1, 2, 2, cs), Miller(0, 1, 0, cs), ... 
    Miller(0, 1, -1, cs, 'uvw'), Miller(-1, 0, -1, cs, 'uvw'));
isTwinning = angle(gB_filtered.misorientation, twinning) < twin_thresh*degree;
twinBoundary = gB_filtered(isTwinning);

% Now add Missing_TB_pairs into twinBoundary
for i = 1:size(Missing_TB_Pairs, 1)
    grain1ID = Missing_TB_Pairs(i, 1);
    grain2ID = Missing_TB_Pairs(i, 2);
    
    % Find boundaries between grain1ID and grain2ID in gB_filtered
    isMissingTB = ((gB_filtered.grainId(:, 1) == grain1ID & gB_filtered.grainId(:, 2) == grain2ID) | ...
                   (gB_filtered.grainId(:, 1) == grain2ID & gB_filtered.grainId(:, 2) == grain1ID));
    
    % Add these missing twin boundary pairs to twinBoundary
    twinBoundary = [twinBoundary; gB_filtered(isMissingTB)];
end

% Optional: Display the final twin boundaries
plot(twinBoundary, 'linecolor', 'r'); % Example to visualize the twin boundaries

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
length_crit = percentageTwin > 0.01;
twinBoundary_new = twinBoundary(length_crit); 

figure; plot(gB); hold on
plot(twinBoundary,twinBoundary.misorientation.angle./degree,'linewidth',2,'displayName','\Sigma3')
mtexColorbar('Special boundaries'); mtexTitle('Corrected twins'); %visualise the twin
plot(twinBoundary_new,'linecolor','r','linewidth',1,'linestyle', '-','displayName','twin boundary')
hold off

%% Reconstruct ParentGrain structure

[mergedGrains,parentId] = merge(grains,twinBoundary_new);
% 
% figure; plot(twinBoundary,'linecolor','r','linewidth',1,'displayName','twin boundary'); hold on
% plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-','displayName','merged grains');
% mtexTitle('Merging Parent grains');hold off

%% Overall Twin fraction (per surface)
length_crit = percentageTwin > 0.01;
twinBoundary_new = twinBoundary(length_crit);
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
figure; plot(grains(~isTwin),'FaceColor','darkgray','displayName','not twin'); hold on
plot(grains(isTwin),'FaceColor','red','displayName','twin'); hold on
plot(mergedGrains.boundary,'linecolor','k','linewidth',2,'linestyle','-','displayName','merged grains');
mtexTitle('twin identification'); 
%text(mergedGrains, int2str(mergedGrains.id));
% Add text label only for the specific merged grain ID
if exist('targetId', 'var') && ~isempty(targetId)
    for i = 1:length(mergedGrains)
        if mergedGrains.id(i) == targetId
            text(mergedGrains(i), int2str(mergedGrains(i).id), 'Color', 'blue'); % Add text for the target grain
            break; % Exit the loop after labeling the specific grain
        end
    end
else
    disp('targetId is not defined or is empty. Skipping this section.');
end
mtexTitle('twin identification'); 
hold off;
save_name = [save_folder filesep alloy '_TwinMap'];
saveFigure([save_name])

%% Parent grains statistics

%name all the columns
columntitle = {'merged ID','merged parent area','M.P. major axis','M.P. aspect ratio','M.P. percent internal boundary length',...
'parent x','parent y','Orientation of parent at xy','GOS of main parent','Count Parent grains','Total Parent Area',...
'Count Twin grains','Total Twinned Area','Total Other Area','Percent of grain that twinned',...
'Schmid for Var1','Schmid for Var2','Schmid for Var3','Schmid for Var4','Schmid for Var5','Schmid for Var6',...
'Area Twinned for Var1','Area Twinned for Var2','Area Twinned for Var3','Area Twinned for Var4','Area Twinned for Var5','Area Twinned for Var6'...
'Rank for Var1','Rank for Var2','Rank for Var3','Rank for Var4','Rank for Var5','Rank for Var6',...
'Twin Count for Var1','Twin Count for Var2','Twin Count for Var3','Twin Count for Var4','Twin Count for Var5','Twin Count for Var6'...
'M.P. total boundary length','M.P. border boundary length','Parent grain','total secondary twins','total Area secondary twins'};

% info on the columns
% merge ID = ID of the parent grain
% MP = merged parent
% MP major axis = big diameter of the elipse (vs minor axis)
% percent internal boundary length: amount of boundary that was recored,
% edge grain have a internal boundary length , 100)
% parent x and y ; position of the barycenter of the grain 
% MP orientation: phi1, PHI, phi2 at xy
% GOS: grain orientation spread in MP
% Count: amount of twins inside the MP


minSize=8; %pixel size allowed for parent grains
%create an array with one line for each of the merged grains, and add the
TwAr = cell(length(mergedGrains), length(columntitle)); 
TwAr = [columntitle;TwAr]; %TwAr is the table where we are going to have all the statistics

%% Remaining code from Magnesium example

%% Set up a matrix with one line per grain

columntitle = {'grain ID','merged ID','grain type','grain area','major axis','minor axis',...
'percent grain boundary internal','deg dev from ideal','GOS','Twin Rank','Active Twin SF'};

GrowCount=length(parentId);
%must be total columns of GrAr-2 to allow for m and parentID
GrAr=cell(GrowCount,length(columntitle)-2);

%the row number is the grain ID so number sequentially
m = reshape(1:GrowCount, [1 GrowCount])';

%concatenate horizontally and add in the parentId
GrAr = [num2cell(m), num2cell(parentId), GrAr];

%label the columns by vertically concatenate with the titles 
GrAr = [columntitle; GrAr];

%% Loop to query parent xy, and retrieve info from grains


r=xvector;

% define hkl and miller indices for FCC matrix

% En1 as extension twin n direction
En1 = [Miller(1, 1, -2, cs, 'uvw'); Miller(-1, -1, 2, cs, 'uvw'); Miller(1, -1, 2, cs, 'uvw'); Miller(-1, 1, -2, cs, 'uvw')];

% Ek1 as extension twin k plane
Ek1 = [Miller(1, 1, 1, cs, 'hkl'); Miller(-1, -1, -1, cs, 'hkl'); Miller(1, -1, 1, cs, 'hkl'); Miller(-1, 1, -1, cs, 'hkl')];

% extension twin axis rotation
Eta=[Miller(-1,2,-1,0,cs,'uvw'); Miller(1,-2,1,0,cs,'uvw'); Miller(-1,-1,2,0,cs,'uvw'); Miller(1,1,-2,0,cs,'uvw'); Miller(2,-1,-1,0,cs,'uvw'); Miller(-2,1,1,0,cs,'uvw')];


% determine the pixel size 
dx = max(ebsd.unitCell(:,1))-min(ebsd.unitCell(:,1));
dy = max(ebsd.unitCell(:,2))-min(ebsd.unitCell(:,2));
pixSize=abs(dx*dy);

%loop through all the merged parent data, once per entry in mergedGrains
for TwinNum=1:length(mergedGrains)
%output the data from the merged grain

%store the merged grain ID for reference
TwAr{(TwinNum+1),1}=mergedGrains(TwinNum).id;
%store the merged area NOTE if you use grain_selected.area this
%incorporates the GB curvature.  Better to use the # of pixels x pixel
%spacing.
area=full(mergedGrains(TwinNum).grainSize)*pixSize;
TwAr{(TwinNum+1),2}=area;

%store the grain major and minor axis.  Since we can easily get the area
%and the aspect ratio, we calculate based on area=pi*minor axis *major axis
%and minor axis * aspect ratio = major axis, so a(major axis/2)=
%area*aspect ratio divided by pi
majAxis=2*sqrt(area*(mergedGrains(TwinNum).aspectRatio)/pi);
%store the major axis length of the parent
TwAr{(TwinNum+1),3}=majAxis;
%store the aspect ratio
TwAr{(TwinNum+1),4}=mergedGrains(TwinNum).aspectRatio;

%store the total boundary length
boundZero=sum(any(mergedGrains(TwinNum).boundary.grainId==0,2));
boundTotal=mergedGrains(TwinNum).boundary.length;
boundPct=100*(boundTotal-boundZero)/boundTotal;
TwAr{(TwinNum+1),5}=boundPct;
TwAr{(TwinNum+1),40}=boundTotal;
TwAr{(TwinNum+1),41}=boundZero;

%Now, determine how many child grains went into this parent
%first put list of rows (which equal id # from grains variable) into an array
[a]=find(parentId==mergedGrains(TwinNum).id);

% the length of this tells how many entries from grains went into the merged grain
%entry (how many times this loop must run)
gLoop=length(a);

%Determine the parent based on which grain referenced by a is the largest
%first, get the area of all grains in a
sizea=full(grains(a).grainSize)*pixSize;
%determine the largest row
[M,I]=max(sizea);
%get the index in corresponding row of a
index = a(I);

%this grain is the parent.
TwAr{(TwinNum+1),42}=index;

%determine the parent orientation by selecting a point of the grain with
 %lowest angle to the meanOrientation, 
%first get id# for all EBSD orientations within that default parent grain
indexAll=find(ebsd.grainId==grains(index).id);
%now get the actual orientations included within default parent grain
oriGindex=ebsd(indexAll).orientations;
%now get the angle between the meanorientation of the parent and all other
%orientations within the parent grain
angG2M=angle(oriGindex, grains(index).meanOrientation)/degree;
%find the ebsd orientation that is closest to the meanOrientation
check1=abs(angG2M);
[angM1,angI1]=min(check1);
%get the x and y of this data
x=ebsd(indexAll(angI1)).x;
y=ebsd(indexAll(angI1)).y;

%Save these values of parent to the TwAr
TwAr{(TwinNum+1),6}=x;
TwAr{(TwinNum+1),7}=y;
oriP=ebsd(x,y).orientations;
[alpha,beta,gamma] = Euler(oriP,'bunge');
TwAr{(TwinNum+1),8} = strcat(num2str(rad2deg(alpha)),',',num2str(rad2deg(beta)),',',num2str(rad2deg(gamma)));

%record the GOS of the parent grain
TwAr{(TwinNum+1),9}=grains(index).GOS./degree;

%prealocate SF
SF=zeros(1,4);
%calculate the parent schmid factor 
for j=1:4
%calculate the angle of the twin plane normal in the parent
ExtTwinNormal=oriP*Ek1(j);
%angle between force and normal
theta=cos(angle(r, ExtTwinNormal));
%calculate the orientation of the slip directions in the parent
ExtTwinDir=oriP*En1(j);
%note that the 'antipodal' command is not used because twinning is
%directional.  this will give negative values of cos (range 1 to -1) when angle >90 degrees
%calculate the angle between the slip and the force
lmda=cos(angle(r,ExtTwinDir));
SF(j)=theta.*lmda;
% Now write schmid factor to array
colWrite=j+15;
TwAr{(TwinNum+1),colWrite}=SF(j);
end
%end of schmid factor calculation loop

%determine the rank of all schmid factors.  sortSF holds the SF in
%ascending order, and inv gives the rank in ascending order (highest schmid
%=6)
[sortSF, extra, invrank] = unique(SF);
%switch rank to descending order, ie highest schmid =1
rank=7-invrank;

%calculate the 6 possible orientations generated by this parent and record
%schmid rank too
for j=1:4
% Write schmid rank to array
colWrite=j+27;
TwAr{(TwinNum+1),colWrite}=rank(j);

%calculate the orientations of the six twin variants
%find the vector3d which is the direction of the twin axis in the parent
vari(j)=oriP*Eta(j);

%set a rotation around that axis
rot(j)=rotation('axis',vari(j),'angle',60*degree);

% rotate the c axis around the twin axis, and save orientation
oriV(j)=rot(j)*oriP;
end

%reset the running totals for number of twins and count
%also count of parent segments and running total of areas which don't match anything (ie errors)
totalAreaT=0;
totalCountT=0;
totalAreaP=0;
totalCountP=0;
totalAreaO=0;
totalCountA=0;
totalAreaA=0;

%these two are arrays of six rows (one per variant)
varCountT=zeros(1,6);
varAreaT=zeros(1,6);

%start loop here through all grains which went into merged grain 
for i=1:gLoop

%get the id of the grain to analyse in this loop from array a
index=a(i);

%and start storing data in GrAr, one line per grain
%and get the number of pixels
N=full(grains(index).grainSize);
%and save it as the area
GrAr{(index+1),4}=N*pixSize;

%store the grain major and minor axis.  Since we can easily get the area
%and the aspect ratio, we calculate based on area=pi*minor axis *major axis
%and minor axis * aspect ratio = major axis, so a(major axis/2)=
%area*aspect ratio divided by pi
majAxis=2*sqrt((N*pixSize)*(grains(index).aspectRatio)/pi);
GrAr{(index+1),5}=majAxis;
GrAr{(index+1),6}=majAxis/grains(index).aspectRatio;

%store the boundary info in the grain array
boundZero=sum(any(grains(index).boundary.grainId==0,2));
boundTotal=grains(index).boundary.length;
boundPct=100*(boundTotal-boundZero)/boundTotal;
GrAr{(index+1),7}=boundPct;

%record the orientation spread
GrAr{(index+1),9}=grains(index).GOS./degree;

%determine if the grain is part of the parent, a twin, or something else
ang=angle(oriP, grains(index).meanOrientation)/degree;
%check if it's a parent
if ang<9.999
%if it is, update GrAr with grain type
GrAr{(index+1),3}='P';
%and record deviation from ideal parent
GrAr{(index+1),8}=ang;

%if the is large enough, add to counters for parent # and %area
if N>minSize
totalCountP=totalCountP+1;
totalAreaP=totalAreaP+N*pixSize;
end
%and exit this loop of grain checking
continue
end

%if the grain is not a parent, we must check for twin boundary
%Extract the GB information to a separate variable
gB = grains(index).boundary;

%To not include the grain-border in calculations, create a new variable for 
%only the mg-mg segments
gB_fcc = gB('Cr','Cr');

%Next we check for each boundary segment whether it is a twinning boundary, 
%i.e., whether boundary misorientation is close to the twinning. 
% restrict to twinnings with threshold 3 degree
isTwinning = angle(gB_fcc.misorientation,twinning) < 3*degree;
%however, this is not a sufficient condition on its own because both sides
%of the boundary will qualify.  The twin must also be within a certain
%range of the allowable twin variants (this prevents 'other' type grains
%from qualifying erroniously)

%compare the 6 variant orientations to the grain in question, 
for j=1:4
angV(j)=angle(oriV(j), grains(index).meanOrientation)/degree;
end
%is the closest variant within 10 degrees?  be careful of negatives
check=abs(angV);
[angM,angI]=min(check);
%if the closest variant is within ~10 degrees of the grain in question,
%this is a twin.  
%only if thre exist a minimum of 3 boundary segments that qualify as twinning and
 %the orientation within 10 degrees of parent will the grain qualify as
 %twin
if sum(isTwinning)>3 && angM<9.999
%if it is, update GrAr
GrAr{(index+1),3}='T';

%and record deviation from ideal twin orientation
GrAr{(index+1),8}=angM;

%and use the index of the lowest angle to retrieve the corresponding schmid
%and rank
GrAr{(index+1),10}=rank(angI);
GrAr{(index+1),11}=SF(angI);

%if the grain has more than minSize, add to counters for twin # and
%area
if N>minSize
%add the area to the running total for twin area
totalAreaT=totalAreaT+N*pixSize;
%increment the twin count.
totalCountT=totalCountT+1;

%if it is a twin, the closest variant is in angI,
%increment the count for that variant

varCountT(angI)=varCountT(angI)+1;
%increment the area for that variant
varAreaT(angI)=varAreaT(angI)+N*pixSize;
%and exit this loop
end
continue
end

%if you're still in the loop, then it's something odd.  Increment area of
%other grains and finish loop if it's big enough to count
GrAr{(index+1),3}='O';
GrAr{(index+1),8}=ang;

if N>minSize
%if the grain is big enough, count the area
totalAreaO=totalAreaO+N*pixSize;
end

end
%end of loop through all grains in a given merged grain

%Check for secondary twins if both other grains and twins (of significant size) exist in
%the mergedGrain
if totalAreaO>0 && totalAreaT>0
%get a list of all applicable other grains and loop through them.
%list all other grains
typea=false(length(a),1);
for i=1:length(a)
row=a(i)+1;
if GrAr{row,3}=='O'
typea(i)=1;
else
typea(i)=0;
end
end
%eliminate all non-others from typea
outO=a(typea);

%for each 'other' grain determine if there is an adjacent twin with which it has a
%misorientation relationship by first checking the boundary orientation and
%relationship.
for i=1:length(outO)
%put the index of the applicable 'other grain' in variable j
indexA=outO(i);
bound_other=grains(indexA).boundary;
%To not include the grain-border in calculations, create a new variable for only the mg-mg segments
gB_fcc = bound_other('Magnesium','Magnesium');
%Next we check for each boundary segment whether it is a twinning boundary, %i.e., whether boundary misorientation is close to the twinning. 
% restrict to twinnings with threshold 3 degree 

%creates a logical of which points are twins
isTwinning = angle(gB_fcc.misorientation,twinning) < 3*degree;
%reduces the data set of boundaries to only include twins
twinBoundary = gB_fcc(isTwinning);
%get the ID of all grains involved
CheckGrains=unique(twinBoundary.grainId);

%remove the 'O' grain from this list so we only have the grains on the
%opposite (outside) of the boundary.
CheckGrains(all(CheckGrains==indexA,2),:)=[];

%for each grain across a twin boundary, check if they are first within the grain and then
%classified as 'T'.  must be this order or may hit uncategorised grains on
%first run through of the original code.
for k=1:length(CheckGrains)
%first check if the grain is in a (which contains the list of all grains
%for this merged grain)
indexT=CheckGrains(k);
if any(a==indexT)
%now check if that grain is T
if GrAr{indexT+1,3}=='T'
%if it is, we've found a secondary twin.  Change the designation of the OTHER GRAIN to A
GrAr{indexA+1,3}='A';
%load both the primay twin (T) and secondary twin (A) orientations into variables.
oriT=grains(indexT).meanOrientation;
oriA=grains(indexA).meanOrientation;
%calculate the possible secondary twin variants 
%determine which one was created
%output the deviation from this orientation to GrAr column 8.

%prealocate SF and verify empty
SF=zeros(1,4);
angV=zeros(1,4);
clear vari
clear rot
clear oriV

%calculate the schmid factors of the secondary twin and the variant
%orientations
for j=1:4
%calculate the angle of the twin plane normal in the parent
ExtTwinNormal= oriT*Ek1(j);
%angle between force and normal
theta=cos(angle(r, ExtTwinNormal));
%calculate the orientation of the slip directions in the parent
ExtTwinDir=oriT*En1(j);
%note that the 'antipodal' command is not used because twinning is
%directional.  this will give negative values of cos (range 1 to -1) when angle >90 degrees
%calculate the angle between the slip and the force
lmda=cos(angle(r,ExtTwinDir));
SF(j)=theta.*lmda;

%calculate the orientations of the six twin variants
%find the vector3d which is the direction of the twin axis in the parent
vari(j)=oriT*Eta(j);

%set a rotation around that axis
rot(j)=rotation('axis',vari(j),'angle',86.3*degree);

% rotate the c axis around the twin axis, and save orientation of the
% variant
oriV(j)=rot(j)*oriT;
%save the deviation of the varinat from the other grain
angV(j)=angle(oriV(j), oriA)/degree;
end
%end of schmid factor calculation loop

check=abs(angV);
[angM,angI]=min(check);
%this is the deviation from an ideal secondary twin
GrAr{(indexA+1),8}=angM;

%determine the rank the active variant.  sortSF holds the SF in
%ascending order, and inv gives the rank in ascending order (highest schmid
%=6, opposite what we want)
[sortSF, extra, invrank] = unique(SF);
%switch rank to descending order, ie highest schmid =1
rank=7-invrank;
%and use the index of the lowest angle to retrieve the corresponding schmid
%and rank
GrAr{(indexA+1),10}=rank(angI);
GrAr{(indexA+1),11}=SF(angI);

%if the secondary twin is large enough add it to the counter and size running totals
N=full(grains(indexA).grainSize);
if N>minSize
totalCountA=totalCountA+1;
%calc area of other grain, and add to running total
totalAreaA=totalAreaA+N*pixSize;
end
%we can now proceed to the next 'other' grain on the list, no more twins
%need to be checked
break
else
continue
end
end
end
end
end
%Check for secondary twinning is now complete.

%write the running totals of perimiter and area to the
%matrix  We must remember to subtract the grains which were changed to
%secondary twins from the 'other' count
totalAreaO=totalAreaO-totalAreaA;
TwAr{(TwinNum+1),10}=totalCountP;
TwAr{(TwinNum+1),11}=totalAreaP;
TwAr{(TwinNum+1),12}=totalCountT;
TwAr{(TwinNum+1),13}=totalAreaT;
TwAr{(TwinNum+1),14}=totalAreaO;
TwAr{(TwinNum+1),43}=totalCountA;
TwAr{(TwinNum+1),44}=totalAreaA;
%only write other twin info if there are twins
if totalCountT>0
%calculate the percent of grain that twinned
TwAr{(TwinNum+1),15}=100*totalAreaT/(totalAreaP+totalAreaT);
%write twin areas by variant
TwAr{(TwinNum+1),22}=varAreaT(1);
TwAr{(TwinNum+1),23}=varAreaT(2);
TwAr{(TwinNum+1),24}=varAreaT(3);
TwAr{(TwinNum+1),25}=varAreaT(4);
TwAr{(TwinNum+1),26}=varAreaT(5);
TwAr{(TwinNum+1),27}=varAreaT(6);
%write twin counts by variant
TwAr{(TwinNum+1),34}=varCountT(1);
TwAr{(TwinNum+1),35}=varCountT(2);
TwAr{(TwinNum+1),36}=varCountT(3);
TwAr{(TwinNum+1),37}=varCountT(4);
TwAr{(TwinNum+1),38}=varCountT(5);
TwAr{(TwinNum+1),39}=varCountT(6);
end

end
%end loop through all merged grains

%%
% parent_area = TwAr.Total_Parent_Area;
% twinned_area = TwAr.Total_Twinned_Area;
% grain_size = TwAr.M_P_major_axis; % Assuming 'M.P. major axis' represents the grain size

%%
% Assuming TwAr is your cell array, with specific columns in it

% Extract the relevant columns from TwAr, skipping the header
merged_ID = TwAr(2:end, 1);           % 'merged ID'
merged_parent_area = cell2mat(TwAr(2:end, 2));   % 'merged parent area' (skip the header)
count_twin_grains = cell2mat(TwAr(2:end, 12));    % 'Count Twin grains' (skip the header)
total_twinned_area = cell2mat(TwAr(2:end, 13));   % 'Total Twinned Area' (skip the header)
total_other_area = cell2mat(TwAr(2:end, 14));     % 'Total Other Area' (skip the header)

% Calculate TwinPercentage
TwinPercentage = 100 * (total_twinned_area + total_other_area) ./ merged_parent_area;

% Initialize TwinPercentage for the first row
TwinPercentage = [TwinPercentage]; % Set first value to NaN or 0

% Create a table with the extracted data and the calculated TwinPercentage
T = table(merged_ID, merged_parent_area, count_twin_grains, total_twinned_area, total_other_area, TwinPercentage, ...
    'VariableNames', {'Merged_ID', 'Merged_Parent_Area', 'Count_Twin_Grains', 'Total_Twinned_Area', 'Total_Other_Area', 'TwinPercentage'});

% Specify the filename to export the data

save_folder = 'Results';
save_name = [save_folder filesep alloy '_TwinFraction.xlsx'];

% Write the table to an Excel file
writetable(T, save_name);

disp(['Data successfully exported to ', save_name]);


%%
% Assuming you have the previous code leading to the creation of T table

% Step 1: Identify the 35% smallest merged parent area grains and exclude them
num_grains = height(T);
num_smallest = floor(0.35 * num_grains); % Calculate 35% of total grains
sorted_grains = sortrows(T, 'Merged_Parent_Area'); % Sort by merged parent area
grains_to_keep = sorted_grains(num_smallest+1:end, :); % Keep only grains excluding the smallest 35%

% Prepare TwinPercentage for remaining grains
total_twinned_area_remaining = grains_to_keep.Total_Twinned_Area; % Total twinned area of remaining grains
total_other_area_remaining = grains_to_keep.Total_Other_Area; % Total other area of remaining grains
merged_parent_area_remaining = grains_to_keep.Merged_Parent_Area; % Merged parent area of remaining grains

% Calculate TwinPercentage for the remaining grains
TwinPercentage_remaining = 100 * (total_twinned_area_remaining + total_other_area_remaining) ./ merged_parent_area_remaining;

% Initialize TwinPercentage for the first row of the remaining grains (optional)
TwinPercentage_remaining = [NaN; TwinPercentage_remaining(2:end)]; % Set first value to NaN or 0 if desired

edges = 0.001:2.5:100;
midpoints = edges(1:end-1) + diff(edges) / 2;

% Step 2: Calculate the counts of TwinPercentage in each bin for remaining grains
[counts, ~] = histcounts(TwinPercentage_remaining, edges); % Calculate counts using only remaining TwinPercentage

% Prepare Y-data for three different versions using the counts from remaining grains
total_twin_grains_remaining = sum(grains_to_keep.Count_Twin_Grains);
total_grains_remaining = height(grains_to_keep); % Total number of grains remaining after excluding the 35%

% Prepare Y-data for three different versions
Y_data_a = counts; % (a) Current counts of twin grains
Y_data_b = counts / total_twin_grains_remaining; % (b) Counts of twin grains / Total amount of twin grains
Y_data_c = counts / total_grains_remaining; % (c) Counts of twin grains / Total number of grains

% Step 3: Create separate bar plots for distribution
% (a) Current counts of twin grains
figure; % Create a new figure
bar(midpoints, Y_data_a, 'FaceColor', 'b'); % Standard blue
xlabel('Twin Percentage Range (%)');
ylabel('Number of Twinned Grains');
title('Distribution of Twin Percentage (Number of Twinned Grains)');
grid on;
xlim([0 100]);
xticks(0:5:100);
save_name = [save_folder filesep alloy '_TwinnedGrainDistribution'];
saveFigure([save_name])

% (b) Counts of twin grains / Total amount of twin grains
figure; % Create a new figure
bar(midpoints, Y_data_b, 'FaceColor', [0.678, 0.847, 0.902]); % Sky blue
xlabel('Twin Percentage Range (%)');
ylabel('Twin Grains / Total Twin Grains');
title('Distribution of Twin Percentage (Twin Grains / Total Twin Grains)');
grid on;
xlim([0 100]);
xticks(0:5:100);
save_name = [save_folder filesep alloy '_NormalisedTwin_distribution'];
saveFigure([save_name])

% (c) Counts of twin grains / Total number of grains
figure; % Create a new figure
bar(midpoints, Y_data_c, 'FaceColor', [0.5, 0.5, 0.5]); % Grey
xlabel('Twin Percentage Range (%)');
ylabel('Twin Grains / Total Grains');
title('Distribution of Twin Percentage (Twin Grains / Total Grains)');
grid on;
xlim([0 100]);
xticks(0:5:100);
save_name = [save_folder filesep alloy '_TotalTwin_distribution'];
saveFigure([save_name])

% Save the figures if desired
% saveas(gcf, [save_folder filesep alloy '_TwinPercentageDistribution_a.png']);
% saveas(gcf, [save_folder filesep alloy '_TwinPercentageDistribution_b.png']);
% saveas(gcf, [save_folder filesep alloy '_TwinPercentageDistribution_c.png']);

%%
% Assuming TwAr is defined and contains the necessary data
merged_ID = TwAr(2:end, 1);           % 'merged ID'
merged_parent_area = cell2mat(TwAr(2:end, 2));   % 'merged parent area' (skip the header)
count_twin_grains = cell2mat(TwAr(2:end, 12));    % 'Count Twin grains' (skip the header)
total_twinned_area = cell2mat(TwAr(2:end, 13));   % 'Total Twinned Area' (skip the header)
total_other_area = cell2mat(TwAr(2:end, 14));     % 'Total Other Area' (skip the header)

% Calculate TwinPercentage
TwinPercentage = 100 * (total_twinned_area + total_other_area) ./ merged_parent_area;

% Initialize TwinPercentage for the first row
TwinPercentage = [TwinPercentage]; % Set first value to NaN or 0

% Create a table with the extracted data and the calculated TwinPercentage
T = table(merged_ID, merged_parent_area, count_twin_grains, total_twinned_area, total_other_area, TwinPercentage, ...
    'VariableNames', {'Merged_ID', 'Merged_Parent_Area', 'Count_Twin_Grains', 'Total_Twinned_Area', 'Total_Other_Area', 'TwinPercentage'});

% Determine the number of grains and the cutoff for the 65% largest grains
num_grains = height(T);
num_largest = floor(0.65 * num_grains); % 65% largest grains

% Sort the table based on the merged parent area and get the remaining 65%
[~, idx_sorted] = sort(T.Merged_Parent_Area, 'descend'); % Sort in descending order
T_remaining = T(idx_sorted(1:num_largest), :); % Keep only the 65% largest grains

% Calculate total twin area for the remaining grains
totalTwinAreaRemaining = T_remaining.Total_Twinned_Area + T_remaining.Total_Other_Area;

% Prepare data for plotting
remainingMergedParentArea = T_remaining.Merged_Parent_Area;

% Create the plot
figure;
hold on;

% Plot the total twin area vs merged parent area
scatter(remainingMergedParentArea, totalTwinAreaRemaining, 'filled', 'MarkerFaceColor', 'blue', 'DisplayName', 'Total Twin Area');
xlabel('grain area (\mu m^2)');
ylabel('twinned area (\mu m^2)');
grid on;

% Optional: add a line of best fit if desired
p = polyfit(remainingMergedParentArea, totalTwinAreaRemaining, 1); % Linear fit
yfit = polyval(p, remainingMergedParentArea);
plot(remainingMergedParentArea, yfit, 'r--', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

% Display the coefficients of the linear fit on the plot
slope = p(1);
text_position = [0.9 0.8]; % Position to place the text relative to the axis
str = sprintf('y = %.2fx ', slope); % Format the text
text('String', str, 'Position', text_position, 'Units', 'normalized', 'FontSize', 10, 'Color', 'r');
hold off;

save_name = [save_folder filesep alloy '_TwinArea_vs_GrainArea'];
saveFigure([save_name])




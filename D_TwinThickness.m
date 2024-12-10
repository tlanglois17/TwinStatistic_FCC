%% Load the inputs/data

save_folder = 'Results';
savename = [save_folder filesep alloy];


%% Build inital fragment graph 
% setMTEXpref('xAxisDirection','east'); % orientation of the map
% setMTEXpref('zAxisDirection','intoPlane');

[G,mGbGrains] = InitialGraph(ebsd,grains,opt);

%Plot full graph
labelNodes=true;labelEdges=false;plotG=true;legendOn=false;
fhandle = plotGraph(grains,[],G,...
    grains.meanOrientation,G.Nodes.Id,...
    labelNodes,labelEdges,legendOn,plotG,[]);
mtexTitle('Full graph')      


%Plot boundary based clusters
figure;
plot(grains,grains.meanOrientation,'noBoundary');hold on;
plot(mGbGrains.boundary,'lineWidth',2,'lineColor','w');
hold off
mtexTitle('Fig0. Boundary based merged grains')  

%figure; plot(mGbGrains,mGbGrains.meanOrientation); mtexTitle('show the GB after merging')


%% Build cluster graph by removing, adding, etc... edges

[G_clust,G,mGrains] = ClusterGraph(G,grains,opt);

%Plot mGrains
figure;plot(grains,grains.meanOrientation,'noBoundary');hold on;
plot(mGrains.boundary,'lineWidth',2,'lineColor','w');
text(mGrains,int2str(mGrains.id));hold off
mtexTitle('Fig1. Cluster graph merged grains')  

%% Edit merged clusters

%plot quantities such as FamilyID that can identify problem clusters
figure;plot(grains,G.Nodes.FamilyID,'noBoundary');hold on;
plot(mGrains.boundary,'lineWidth',2,'lineColor','k');
text(mGrains,int2str(mGrains.id));hold off;
mtexTitle('Fig 2. Cluster graph Families')  

% groups=unique(G_clust.Nodes.Group(G_clust.Nodes.FamilyID>5));
groups=[]
value=G.Nodes.FamilyID;
GraphEditor(groups,1,[],G_clust,G,grains,mGrains,value,0,1,0,1,2);

%% Build family graph for each cluster

computemGrainId=[];
computemGrainId=mGrains.id;
[G_Family,G_clust,G] = FamilyGraph(G_clust,G,grains,computemGrainId,opt);


%% Clean the family tree

cleanGroups=unique(G_Family.Nodes.Group);
[G_Family,G_clust,G,exflagGroup]= CleanFamilyTree(cleanGroups,...
    G_Family,G_clust,G,grains,mGrains,opt);   


%% Visualize the results

%Plot grain fragment twin type
figure; plot(grains,G.Nodes.type,'noBoundary');
hold on;plot(mGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-','displayName',...
    'merged grains');hold off
mtexColorbar;mtexTitle('Twin Type');
saveFigure([savename '_Twin type'])

%Plot the number of generations in each cluster
figure; plot(grains,G.Nodes.Generation,'noBoundary');
hold on;plot(mGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-','displayName',...
    'merged grains');hold off;
mtexColorbar;mtexTitle('Generation');
saveFigure([savename '_Twin Generation'])

figure; plot(grains,G.Nodes.EffSF,'noBoundary');
hold on;plot(mGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-','displayName',...
    'merged grains');hold off;
mtexColorbar;mtexTitle('EFF SF');

figure; plot(grains,G.Nodes.nSFAV,'noBoundary');
hold on;plot(mGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-','displayName',...
    'merged grains');hold off;
mtexColorbar;mtexTitle('SF active variant');

figure; plot(grains,G.Nodes.nSFAVR,'noBoundary');
hold on;plot(mGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-','displayName',...
    'merged grains');hold off;
mtexColorbar;mtexTitle('Active variant SF rank (lower is most stresses)');

%% Edit/Visualize the Family tree

%Use the exit flag to construct groups to manually address
%exflagGroup=-1 single frag
%exflagGroup=1 no twin relationships but multiple families
%exflagGroup=2 too many parents
%exflagGroup=3 no parent but twin relationships exist
%exflagGroup=4 max number of generations hit while assigning generation
%exflagGroup=5 not all fragments were related... something is wrong
%exflagGroup=6 entered fix circular relationship fnc too many times
groups=cleanGroups(exflagGroup>0)
% groups=unique(G_clust.Nodes.Group(find(G_clust.Nodes.computeFamily)))
groups=[]
value=G.Nodes.FamilyID;
[G_Family] = GraphEditor(groups,1,G_Family,G_clust,grains,mGrains,value,0,1,0,1,1);


%% Compute the twin fraction 

[mGrains,twinVF] = getTwinFractions(G,grains,mGrains,opt);


%% Twin thickness 

[G] = TwinThickness(G,grains,opt);
figure;
h = histogram(G.Nodes.twinThickness(G.Nodes.twinThickness>0), 'BinWidth', 0.01);
xlabel('Twin Fragment thickness (um)')
ylabel('Counts')


% Extract histogram data
binEdges = h.BinEdges; % Bin edges
binCounts = h.Values;   % Counts for each bin

% Calculate mean and median of twin thickness
twinThicknessValues = G.Nodes.twinThickness(G.Nodes.twinThickness > 0);
meanThickness = mean(twinThicknessValues);
medianThickness = median(twinThicknessValues);

% Plot mean and median lines
hold on; % Retain current graph
yLimits = ylim;
xLimits = xlim; % Get current y-axis limits
plot([meanThickness meanThickness], yLimits, 'r--', 'LineWidth', 2, 'DisplayName', 'Mean'); % Mean line
plot([medianThickness medianThickness], yLimits, 'g--', 'LineWidth', 2, 'DisplayName', 'Median'); % Median line


xCenter = mean(xLimits); % Center x position between mean and median
yCenter = mean(yLimits); % Center y position

% Annotate the mean and median values on the graph
text(xCenter, yCenter, ['Mean: ', num2str(meanThickness * 1000), ' nm'], ...
     'Color', 'r', 'HorizontalAlignment', 'center', 'FontSize', 10);
text(xCenter, 0.8*yCenter, ['Median: ', num2str(medianThickness * 1000), ' nm'], ...
     'Color', 'g', 'HorizontalAlignment', 'center', 'FontSize', 10);

hold off;

% Add a legend
legend('Histogram', 'Mean', 'Median');

saveFigure([savename '_Twin Fragment Thickness'])

% Prepare the data for export
dataToExport = table(binEdges(1:end-1)', binCounts', ...
                     'VariableNames', {'Bin_Edges', 'Counts'});

% Write data to an Excel file
writetable(dataToExport, [savename '_TwinFragmentThicknessData.xlsx']);
%% Twin count

[mGrains] = CountTwins(G,grains,mGrains,opt);

figure;histogram(mGrains.prop.twinCount,'BinMethod','integers');
xlabel('Number of twin fragments in cluster')
ylabel('Counts')
saveFigure('Twin Fragments in cluster')
figure;histogram(mGrains.prop.twinFamilyCount,'BinMethod','integers');
xticks([0,1,2]); 
xlabel('Unique twin fragments in cluster')
ylabel('Counts')
saveFigure([savename '_Unique twin fragments in cluster'])


%% Twin variant ranking

figure;histogram(G.Nodes.nSFAVR(G.Nodes.nSFAVR>0),'BinLimits',[0.5,6.5],'BinMethod','integers');
xlabel('Rank (1 is largest schmid factor)')
ylabel('Counts')
xticks([1,2,3,4,5,6]); 
saveFigure([savename '_Twin Variant Rank'])


%% Transfer results to grains 

[grains] = transferGtoGrains(G,grains);


%% Save important variable
%save('Segmented_data.mat','G','grains','mGrains','opt','twinVF');


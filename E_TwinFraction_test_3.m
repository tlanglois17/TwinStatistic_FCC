%% Table grain statistic

% A. Create the table
columntitle = {
    'MG ID','MG area',... % 1-2: ID and total area of merged grains
    'MP count','MP Area',... % 3-4: Amount of parent regions and areas
    'TW count','TW area','Twinning % in MG',... % 5-7
    'GB+TB length','GB length',... % 8-9
    'TB length', 'TB % length',... % 10-11
    'TW thickness',... % 12
    'count 2nd twins','Area 2nd twins' % 14-15
    };
TwAr = cell(length(mergedGrains), length(columntitle)); %create an array
TwAr = [columntitle;TwAr]; % add labels

% B. Initiate calculation 
dx = max(ebsd.unitCell(:,1))-min(ebsd.unitCell(:,1));
dy = max(ebsd.unitCell(:,2))-min(ebsd.unitCell(:,2));
pixSize=abs(dx*dy); % determine the pixel size 


for i=1:length(mergedGrains) % statistic for all the Merged grains
    
    % 1-2: MG Id and areas    
    TwAr{(i+1),1}=mergedGrains(i).id;
    area=full(mergedGrains(i).grainSize)*pixSize; % no GB into accounts is better
    TwAr{(i+1),2}=area; 
    
    % 3-4-5-6-7: MP and TB distinction, count and areas
    % gArea = grains_v1.area; % this definition incorporates grain boundaries to the area count
    gArea = full(grains_v1.grainSize)*pixSize; % no GB into accounts is better
    childId = find(parentId==i);
    mpId = childId(~isTwin(childId));
    twId = childId(isTwin(childId));
    mpCount = length(mpId);
    mpArea = sum(gArea(mpId));
    twCount = length(twId);
    twArea = sum(gArea(twId));
    twLength = 0;
    for k=1:length(twId)
        twLength = twLength + grains_v1(twId(k)).boundary.length;
    end
    TwAr{(i+1),3}=mpCount;
    TwAr{(i+1),4}=mpArea;
    TwAr{(i+1),5}=twCount;
    TwAr{(i+1),6}=twArea;
    TwAr{(i+1),7}=100*twArea/area;
    
    % 8-9-10-11: store boundary lengths
       
    boundGB=mergedGrains(i).boundary.length;
    boundTB=twLength;
    boundGBTB = boundGB + boundTB; %all TB are parts of two grains 
    boundPct=100*(boundGBTB-boundGB)/boundGBTB;
    TwAr{(i+1),8}=boundGBTB;
    TwAr{(i+1),9}=boundGB;
    TwAr{(i+1),10}=boundTB;
    TwAr{(i+1),11}=boundPct;
    
    % 12: twin thickness
    twinThicknessValues = calculateTwinThickness(grains_v1(twId));
    %twinThicknessValues = calculateTwinThickness(grains_v1(childId));
    TwAr{(i+1), 12} = twinThicknessValues; 
end

savename = [save_folder filesep alloy '_TwAr.mat'];

save(savename,"TwAr");

%% Plots

% 1. Twin Fraction distribution 
% Extract the relevant columns from TwAr, skipping the header
merged_ID = TwAr(2:end, 1);          
merged_grain_area = cell2mat(TwAr(2:end, 2));  
parent_area = cell2mat(TwAr(2:end, 4));
count_twins = cell2mat(TwAr(2:end, 5));     
twinned_area = cell2mat(TwAr(2:end, 6));   
twin_percentage = cell2mat(TwAr(2:end, 7));
twin_length = cell2mat(TwAr(2:end, 10));
twin_thickness = cell2mat(TwAr(2:end, 12));


T = table(merged_ID, merged_grain_area, parent_area, count_twins, twinned_area,twin_percentage, twin_length, ...
    'VariableNames', {'Merged_ID', 'Tot_Area','Parent_Area', 'Count_Twins', 'Twin_Area', 'Twin_%', 'Twin_Length'});

% savename = [save_folder filesep alloy '_TwinFraction.xlsx'];
% writetable(T, savename);
% disp(['Data successfully exported to ', savename]);

num_grains = height(T);
num_smallest = floor(grain_dist_thresh * num_grains); % Calculate 35% of total grains
sorted_grains = sortrows(T, 'Tot_Area'); % Sort by merged area
grains_to_keep = sorted_grains(num_smallest+1:end, :); % Keep only grains excluding the smallest 35%

% Prepare stats for remaining grains
total_twinned_area_remaining = grains_to_keep.Twin_Area; % Total twinned area of remaining grains
grain_area_remaining = grains_to_keep.Tot_Area; % Merged parent area of remaining grains
TwinPercentage_remaining = grains_to_keep.("Twin_%");
TwinLength_remaining = grains_to_keep.Twin_Length;

% Initialize TwinPercentage for the first row of the remaining grains (optional)
TwinPercentage_remaining = [NaN; TwinPercentage_remaining(2:end)]; % Set first value to NaN or 0 if desired

edges = 0.001:2.5:100;      
midpoints = edges(1:end-1) + diff(edges) / 2;

% Step 2: Calculate the counts of TwinPercentage in each bin for remaining grains
[counts, ~] = histcounts(TwinPercentage_remaining, edges); % Calculate counts using only remaining TwinPercentage

% Prepare Y-data for three different versions using the counts from remaining grains
total_twin_grains_remaining = sum(grains_to_keep.Count_Twins);
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



% 2. Twin Thickness

twinThicknessVals = cell2mat(TwAr(2:end, 12));  % Extract values from the table, excluding header
figure;
h = histogram(twinThicknessVals(twinThicknessVals > 0), 'BinWidth', 0.01);
xlabel('Twin Thickness'); ylabel('Count'); title('Distribution of Twin Thickness');
hold on;

meanThickness = mean(twinThicknessVals);
medianThickness = median(twinThicknessVals);
yMax = max(h.Values);
ylim([0, yMax + 1]);
yLimits = ylim;
xLimits = xlim;
xCenter = mean(xLimits); % Center x position between mean and median
yCenter = mean(yLimits); % Center y position

plot([meanThickness meanThickness], yLimits, 'r--', 'LineWidth', 2, 'DisplayName', 'Mean'); % Mean line
plot([medianThickness medianThickness], yLimits, 'g--', 'LineWidth', 2, 'DisplayName', 'Median'); % Median line

text(xCenter, yCenter, ['Mean: ', num2str(meanThickness * 1000), ' nm'], ...
     'Color', 'r', 'HorizontalAlignment', 'center', 'FontSize', 10);
text(xCenter, 0.8*yCenter, ['Median: ', num2str(medianThickness * 1000), ' nm'], ...
     'Color', 'g', 'HorizontalAlignment', 'center', 'FontSize', 10);
legend('Histogram', 'Mean', 'Median');
hold off;
saveFigure([savename '_Twin Fragment Thickness'])


% 3. Twin Length
merged_grain_area = cell2mat(TwAr(2:end, 2));
twin_length = cell2mat(TwAr(2:end, 10));
figure;
hold on
scatter(merged_grain_area, twin_length, 'filled', 'MarkerFaceColor', 'blue', 'DisplayName', 'TwinBoundaries length per Grain Area');
xlabel('grain area');
ylabel('sum twin length');
grid on;
p = polyfit(merged_grain_area, twin_length, 1); % Linear fit
yfit = polyval(p, merged_grain_area);
plot(merged_grain_area, yfit, 'r--', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

% Display the coefficients of the linear fit on the plot
slope = p(1);
text_position = [0.9 0.8]; % Position to place the text relative to the axis
str = sprintf('y = %.2fx ', slope); % Format the text
text('String', str, 'Position', text_position, 'Units', 'normalized', 'FontSize', 10, 'Color', 'r');
hold off;

save_name = [save_folder filesep alloy '_TwinLength_vs_GrainArea'];
saveFigure([save_name])



function twinThickness = calculateTwinThickness(grains)
    % Get minimum axis lengths using fitEllipse (assumed function)
    [~, ~, b_min] = fitEllipse(grains);  % Assuming fitEllipse is defined and returns b_min
    twinThickness = b_min;  % Use the minimum axis length directly as the twin thickness
end





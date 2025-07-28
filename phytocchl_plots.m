% new plotting script for new phytocchl table 7/24/25 - EM

% assigning each row a season designation
m = datetime(cchl100.date_sampled, 'InputFormat','yyyy-MM-dd HH:mm:ss+00:00');

seasons = table;
seasons.month = month(m, 'shortname');
seasons.season = strings(size(seasons.month, 1), 1);

for i = 1:height(seasons.month)
    if ismember(seasons.month(i), {'Dec','Jan', 'Feb'}) == 1
        seasons.season(i) = 'Win';
    end
    if ismember(seasons.month(i), {'Mar', 'Apr', 'May'}) == 1
        seasons.season(i) = 'Spr';
    end
    if ismember(seasons.month(i), {'Jun', 'Jul', 'Aug'}) == 1
        seasons.season(i) = 'Sum';
    end
    if ismember(seasons.month(i), {'Sep', 'Oct', 'Nov'}) == 1
        seasons.season(i) = 'Fal';
    end
end
clear m

%% boxplot all data, including outliers using breakyaxis
%7/24/25

boxplot(phytocchl.C_chl_ratio_total, seasons.month)
ylim([-20 2700])

[q, ~, uniqueid] = unique(seasons.month(~isnan(phytocchl.C_chl_ratio_total)));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -20, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end

xlabel('Month')
ylabel('C:chl a')
title('C:chl a for all samples')
breakyaxis([900 2500]) % this is from an imported package

%clear q counts uniqueid qsorted labels xt

%% boxplot all data, including outliers, with max y value for outliers
% 7/24/25
boxplot(phytocchl.C_chl_ratio_total, seasons.month)
ymax = 500;
ylim([-20 ymax])
hold on

% Get unique group labels
labels = unique(seasons.month, 'stable');  % 'stable' preserves original order

% Loop through and plot capped outliers
for i = 1:length(phytocchl.C_chl_ratio_total)
    val = phytocchl.C_chl_ratio_total(i);
    if val > ymax
        % Find x-position by matching label
        xpos = find(strcmp(seasons.month{i}, labels));
         % Apply small horizontal jitter
        jitterAmount = 0.25;  % adjust for more/less spread
        xjitter = xpos + (rand - 0.5) * 2 * jitterAmount;
        
        % Plot the capped outlier
        plot(xjitter, ymax, '^', ...
            'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', 6)
    end
end

[q, ~, uniqueid] = unique(seasons.month(~isnan(phytocchl.C_chl_ratio_total)));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -15, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end

xlabel('Month')
ylabel('C:chl a')
title('C:chl a for all samples')

%% scatter by month
seasons.monthcat = categorical(seasons.month, {'Jan','Feb', 'May','Jul', 'Aug', 'Sep', 'Oct', 'Nov'}, 'Ordinal', 1);

scatter(seasons.monthcat, phytocchl.C_chl_ratio_total, 40, phytocchl.depth_m, 'filled',  'o', 'jitter', 'on', 'jitteramount', 0.15) % can scatter actual data points on

colormap(flipud(parula))
c = colorbar;
c.Direction = 'reverse';
xlabel('Month')
ylabel('C:chl a')
set(gca, 'FontSize', 14, 'LineWidth', 1)
title("C:chl a total fraction by month", 'FontWeight', 'bold')

ymax = 500;
ylim([-20 ymax])
hold on

% Get unique group labels
labels = unique(seasons.month, 'stable');  % 'stable' preserves original order

% Loop through and plot capped outliers
for i = 1:length(phytocchl.C_chl_ratio_total)
    val = phytocchl.C_chl_ratio_total(i);
    if val > ymax
        % Find x-position by matching label
        xpos = find(strcmp(seasons.month{i}, labels));
         % Apply small horizontal jitter
        jitterAmount = 0.25;  % adjust for more/less spread
        xjitter = xpos + (rand - 0.5) * 2 * jitterAmount;
        
        % Plot the capped outlier - can't figure out how to get these to
        % plot with the same color scheme as colorbar
        scatter(xjitter, ymax, 50, phytocchl.depth_m(i), 'filled')
    end
end

%% plotting boxplot and scatter on top of each other - 7/24 finally works!
figure(3)
ymax = 500;

[~, labels] = findgroups(seasons.monthcat);  % group index for each row

ax1 = axes();
scatter(ax1, seasons.monthcat, phytocchl.C_chl_ratio_total, 40, phytocchl.depth_m, "filled", 'o', 'jitter', 'on', 'jitteramount', 0.1)
set(ax1, 'YLim', [-20 ymax])
pos = get(ax1, 'Position');
colormap(flipud(parula))
c = colorbar;
c.Label.String = 'Depth (m)';
c.Location = 'eastoutside';
c.Direction = 'reverse';
set(ax1, 'Position', pos)
hold on

% Find capped outliers
isCapped = phytocchl.C_chl_ratio_total > ymax;

% Convert categorical groups to numeric positions
[~, xGroupNums] = ismember(seasons.monthcat, categories(seasons.monthcat));

% Extract x positions and apply jitter
jitterAmt = 0.25;
xCapped = xGroupNums(isCapped) + (rand(sum(isCapped),1) - 0.5) * 2 * jitterAmt;

% All y-values get capped to ymax
yCapped = repmat(ymax, sum(isCapped), 1);

% Extract color variable
colorCapped = phytocchl.depth_m(isCapped);

% Plot capped outliers
scatter(ax1, xCapped, yCapped, 50, colorCapped, ...
    'filled', 'Marker', 'o', 'MarkerEdgeColor', 'k')

ax2 = axes('Position', get(ax1, 'Position'), 'YLim', get(ax1, 'YLim'),  'PlotBoxAspectRatio', get(ax1, 'PlotBoxAspectRatio'), 'Visible','off', 'Color', 'none');
hold(ax2, 'on')
boxplot(ax2, phytocchl.C_chl_ratio_total, seasons.monthcat, 'Positions', 1:length(labels))
set(ax2, 'Ylim', get(ax1, 'Ylim'))
set(ax2, 'Position', get(ax1, 'position'))

labels2 = unique(seasons.month, 'stable');  % 'stable' preserves original order

% Get unique group labels
% Loop through and plot capped outliers
for i = 1:length(phytocchl.C_chl_ratio_total)
    val = phytocchl.C_chl_ratio_total(i);
    if val > ymax
        % Find x-position by matching label
        xpos = find(strcmp(seasons.month{i}, labels2));
         % Apply small horizontal jitter
        jitterAmount = 0.25;  % adjust for more/less spread
        xjitter = xpos + (rand - 0.5) * 2 * jitterAmount;
        
        % Plot the capped outlier - can't figure out how to get these to
        % plot with the same color scheme as colorbar
        scatter(ax2, xjitter, ymax, 50, phytocchl.depth_m(i), 'filled')
    end
end

ax2.XLim = [min(double(seasons.monthcat)) - 0.5, max(double(seasons.monthcat)) + 0.5];
ax2.XTick = 1:numel(categories(seasons.monthcat));
ax2.XTickLabel = categories(seasons.monthcat);

[q, ~, uniqueid] = unique(seasons.month(~isnan(phytocchl.C_chl_ratio_total)));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -15, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end

xlabel(ax1, 'Month')
ylabel(ax1, 'C:chl a')
title(ax1, 'C:chl a by month, total size fraction')
set(ax1, 'Fontsize', 14)

clearvars -except phytocchl seasons

%% boxplot + scatter for only upper 100m of ocean
figure(3)
ymax = 500;
req = phytocchl.depth_m <= 100;
[~, labels] = findgroups(seasons.monthcat);  % group index for each row

ax1 = axes();
scatter(ax1, seasons.monthcat(req), phytocchl.C_chl_ratio_total(req), 40, phytocchl.depth_m(req), "filled", 'o', 'jitter', 'on', 'jitteramount', 0.1)
set(ax1, 'YLim', [-20 ymax])
pos = get(ax1, 'Position');
colormap(flipud(parula))
c = colorbar;
c.Label.String = 'Depth (m)';
c.Location = 'eastoutside';
c.Direction = 'reverse';
set(ax1, 'Position', pos)
hold on

% Find capped outliers
isCapped = phytocchl.C_chl_ratio_total(req) > ymax;

% Convert categorical groups to numeric positions
[~, xGroupNums] = ismember(seasons.monthcat(req), categories(seasons.monthcat(req)));

% Extract x positions and apply jitter
jitterAmt = 0.25;
xCapped = xGroupNums(isCapped) + (rand(sum(isCapped),1) - 0.5) * 2 * jitterAmt;

% All y-values get capped to ymax
yCapped = repmat(ymax, sum(isCapped), 1);

% Extract color variable
colorCapped = phytocchl.depth_m(isCapped);

% Plot capped outliers
scatter(ax1, xCapped, yCapped, 50, colorCapped, ...
    'filled', 'Marker', 'o', 'MarkerEdgeColor', 'k')

ax2 = axes('Position', get(ax1, 'Position'), 'YLim', get(ax1, 'YLim'),  'PlotBoxAspectRatio', get(ax1, 'PlotBoxAspectRatio'), 'Visible','off', 'Color', 'none');
hold(ax2, 'on')
boxplot(ax2, phytocchl.C_chl_ratio_total(req), seasons.monthcat(req), 'Positions', 1:length(labels))
set(ax2, 'Ylim', get(ax1, 'Ylim'))
set(ax2, 'Position', get(ax1, 'position'))

labels2 = unique(seasons.month(req), 'stable');  % 'stable' preserves original order

% Get unique group labels
% Loop through and plot capped outliers
for i = 1:length(phytocchl.C_chl_ratio_total(req))
    val = phytocchl.C_chl_ratio_total(req(i));
    if val > ymax
        % Find x-position by matching label
        xpos = find(strcmp(seasons.month{i}, labels2));
         % Apply small horizontal jitter
        jitterAmount = 0.25;  % adjust for more/less spread
        xjitter = xpos + (rand - 0.5) * 2 * jitterAmount;
        
        % Plot the capped outlier - can't figure out how to get these to
        % plot with the same color scheme as colorbar
        scatter(ax2, xjitter, ymax, 50, phytocchl.depth_m(i), 'filled')
    end
end

ax2.XLim = [min(double(seasons.monthcat(req))) - 0.5, max(double(seasons.monthcat(req))) + 0.5];
ax2.XTick = 1:numel(categories(seasons.monthcat(req)));
ax2.XTickLabel = categories(seasons.monthcat(req));

[q, ~, uniqueid] = unique(seasons.month(req & ~isnan(phytocchl.C_chl_ratio_total)));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -15, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end

xlabel(ax1, 'Month')
ylabel(ax1, 'C:chl a')
title(ax1, 'C:chl a by month, total size fraction, upper 100m')
set(ax1, 'Fontsize', 14)

clearvars -except phytocchl seasons
%% plotting by in/mid/offshore and 0-50, 50-100, 100+ depth

% first need to index where each sample is coming from (in/mid/offshore)

distfromshore = strings(size(cchl100.nearest_station, 1), 1);
for i = 1:height(distfromshore)
    if ismember(cchl100.nearest_station(i), {'MVCO','L1', 'L2'}) == 1
        distfromshore(i) = 'innershelf';
    end
    if ismember(cchl100.nearest_station(i), {'L3', 'L4','L5','L6', 'L7', 'L8', 'L9'}) == 1
        distfromshore(i) = 'midshelf';
    end
    if ismember(cchl100.nearest_station(i), {'L10', 'L11'}) == 1
        distfromshore(i) = 'outershelf';
    end
end

%month_num = month(phytocchl.date_sampled);
%monthOrder = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
%month_categories = categorical(month_num, 1:12, monthOrder);

%%
monthorder = {'Jan','Feb', 'May','Jul','Aug','Sep','Oct', 'Nov'};
% then plot as subplots

% 0-50 m depth (1-3)
subplot(3, 3, 1)

loc = 'innershelf'; % 'innershelf', 'midshelf', 'outershelf' are options
frac = 'total';
dotidx = append("C_chl_ratio_", frac);

req = (phytocchl.depth_m <= 50 &  ~isnan(phytocchl.(dotidx)) ...
    & ~isoutlier(phytocchl.(dotidx), "median") &  (distfromshore == loc));
unique(seasons.month(req))

boxplot(phytocchl.(dotidx)(req), seasons.month(req), 'GroupOrder', monthorder)

xlabel('Month')
ylabel('C:chl a')
ylim([-20 Inf])
title(append(frac, '-', loc, '-', 'upper50m'))
[q, ~, uniqueid] = unique(seasons.month(req));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [2 1 4 3 6 5]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');
%%
for i = 1:length(labels)
    text(xt(i), -15, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end

clear q counts uniqueid qsorted labels xt frac dotidx loc req monthorder

%% 50-100 m depth (plots 4-6)
subplot(3, 3, 6) % 4 will be empty bc MVCO-L2 don't have depth > 50 m
monthorder = {'Jan','Feb','Jul', 'Aug', 'Sep','Oct', 'Nov'};
loc = 'outershelf'; %'innershelf', 'midshelf', 'outershelf' are options
frac = 'greaterthan20';
dotidx = append("C_chl_ratio_", frac);

req = (phytocchl.depth_m > 50 & phytocchl.depth_m <= 100 & ~isnan(phytocchl.(dotidx)) ...
    & ~isoutlier(phytocchl.(dotidx), "median") & (distfromshore == loc));
unique(seasons.month(req))

boxplot(phytocchl.(dotidx)(req), seasons.monthcat(req), 'GroupOrder', monthorder)

xlabel('Month')
ylabel('C:chl a')
%ylim([-20 Inf])
title(append(frac, '-', loc, '-', '50to100m'))
[q, ~, uniqueid] = unique(seasons.month(req));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 4 1 7 6 5]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');
%%
for i = 1:length(labels)
    text(xt(i), -30, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt frac dotidx loc req

%% >100m m depth (plots 8-9)
subplot(3, 3, 9) % 7 will be empty bc MVCO-L2 don't have depth > 50 m
monthorder = {'Feb', 'Aug'}; %, 'Aug', 'Sep', 'Oct','Nov'
loc = 'outershelf'; %'innershelf', 'midshelf', 'outershelf' are options
frac = 'greaterthan20';
dotidx = append("C_chl_ratio_", frac);

req = (phytocchl.depth_m > 100 & ~isnan(phytocchl.(dotidx)) ...
    & ~isoutlier(phytocchl.(dotidx), "median") & (distfromshore == loc));
unique(seasons.month(req))

boxplot(phytocchl.(dotidx)(req), seasons.monthcat(req), 'GroupOrder', monthorder)

xlabel('Month')
ylabel('C:chl a')
ylim([-20 Inf])
title(append(frac, '-', loc, '-', 'over100m'))
[q, ~, uniqueid] = unique(seasons.month(req));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [2 1] ;%[2 1 5 4 3]; % put months in order (aug feb nov oct sep)
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -10, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt frac dotidx loc req

%% 7/28 plotting by shelf distance and depth, c vs. chl by season

cmin = min(cchl100.depth_m);
cmax = max(cchl100.depth_m);
subplot(1, 3, 3)
loc = 'outershelf'; %'innershelf', 'midshelf', 'outershelf' are options
chlfrac = '0_avg'; % options are 0_avg, 20_avg, lessthan5, 5to10, 10to20
chllabel = append('chl_', chlfrac);
cfrac = 'total'; % options are total, 5, 10, 5to10, 10to20, greaterthan20
sumClabel = append('sumC_', cfrac);
scatter(cchl100.(chllabel)(seasons.season == 'Win' & distfromshore == loc), cchl100.(sumClabel)(seasons.season == 'Win' & distfromshore == loc), 40, cchl100.depth_m(seasons.season == 'Win' & distfromshore == loc), 'filled', 'o')
hold on
scatter(cchl100.(chllabel)(seasons.season == 'Spr' & distfromshore == loc), cchl100.(sumClabel)(seasons.season == 'Spr' & distfromshore == loc), 40, cchl100.depth_m(seasons.season == 'Spr' & distfromshore == loc), 'filled', 's')
colormap(flipud(parula))
hold on
scatter(cchl100.(chllabel)(seasons.season == 'Sum' & distfromshore == loc), cchl100.(sumClabel)(seasons.season == 'Sum' & distfromshore == loc), 40, cchl100.depth_m(seasons.season == 'Sum' & distfromshore == loc), 'filled', '^')
hold on
scatter(cchl100.(chllabel)(seasons.season == 'Fal' & distfromshore == loc), cchl100.(sumClabel)(seasons.season == 'Fal' & distfromshore == loc), 40, cchl100.depth_m(seasons.season == 'Fal' & distfromshore == loc), 'filled', 'p')
colormap(flipud(parula))


legend({'Winter', 'Spring','Summer', 'Fall', '', ''}, "Location","northeast")
set(gca, "FontSize", 14)
xlim([0 10])
ylim([0 350])
clim([cmin cmax])
xlh = xlabel('Chl a (µg L^{-1})');
ylabel('Carbon (µg L^{-1})')
title(loc)
xl = xlim;
line(xl, xl*50, 'LineWidth', 1)
line(xl, xl*200, 'LineWidth', 1)
legend({'Winter', 'Spring','Summer', 'Fall', '', ''}, "Location","northeast")

%% do these at the end
c = colorbar;
c.Direction = 'reverse';
c.Label.String = 'depth (m)';
sgtitle('C and chl across shelf')
%set(gca, "PlotBoxAspectRatio", get(subplot(1, 3, 2), 'PlotBoxAspectRatio'))


%% histograms of c:chl dist for diff size classes across shelf
subplot(3, 5, 15)
loc = 'outershelf';
frac = 'greaterthan20';
ratio = append('C_chl_ratio_', frac);
input = cchl100.(ratio)(distfromshore == loc & ~isoutlier(cchl100.(ratio)));
histogram(input, -200:5:200)
xlabel('C:chl a')
ylabel('Count')
title(append(loc, '-', frac))

%% 7/27 scatters of c vs. chl for diff size classes across shelf
tiledlayout(3, 4)
cmin = min(cchl100.depth_m);
cmax = max(cchl100.depth_m);

nexttile;
loc = 'outershelf';
frac = 'greaterthan20';
chlfrac = 'chl_20_avg';
cfrac = 'sumC_greaterthan20';
scatter(cchl100.(chlfrac)(distfromshore == loc), cchl100.(cfrac)(distfromshore == loc), 40, cchl100.depth_m(distfromshore == loc), 'filled')
xlabel('Chl a (µg/L)')
ylabel('Carbon (µg/L)')
xlim([-2 4])
ylim([0 150])
clim([cmin cmax])
xl = xlim;
line(xl, xl*50, 'LineWidth', 1)
line(xl, xl*200, 'LineWidth', 1)
title(append(loc, '-', frac))

c = colorbar;
c.Layout.Tile = 'east';
c.Direction = 'reverse';
c.Label.String = 'depth';


%% 7/27 by season, totals, all samples 
subplot(2, 2, 1)
szn = 'Win';
req = (seasons.season == szn & ~isoutlier(cchl100.chl_0_avg, 'median'));
scatter(cchl100.chl_0_avg(req), cchl100.sumC_total(req), 40, cchl100.depth_m(req), 'filled', 'o');
colormap(flipud(parula))
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
set(gca, 'FontSize', 16)
title('Winter')
ylim([0 300])
xlim([0 5])
xl = xlim;
line(xl, xl*50, 'LineWidth', 1)
line(xl, xl*200, 'LineWidth', 1)
legend('data', 'c:chla = 50', 'c:chla = 200')
hold on
clear req szn

subplot(2, 2, 2)
szn = 'Spr';
req = (seasons.season == szn & ~isoutlier(cchl100.chl_0_avg, 'median'));
scatter(cchl100.chl_0_avg(req), cchl100.sumC_total(req), 40, cchl100.depth_m(req), 'filled', 'o');
colormap(flipud(parula))
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
set(gca, 'FontSize', 16)
title('Spring')
ylim([0 300])
xlim([0 5])
xl = xlim;
line(xl, xl*50, 'LineWidth', 1)
line(xl, xl*200, 'LineWidth', 1)
legend('data', 'c:chla = 50', 'c:chla = 200')
hold on
clear req szn

subplot(2, 2, 3)
szn = 'Sum';
req = (seasons.season == szn & ~isoutlier(cchl100.chl_0_avg, 'median'));
scatter(cchl100.chl_0_avg(req), cchl100.sumC_total(req), 40, cchl100.depth_m(req), 'filled', 'o');
colormap(flipud(parula))
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
set(gca, 'FontSize', 16)
title('Summer')
ylim([0 300])
xlim([0 5])
xl = xlim;
line(xl, xl*50, 'LineWidth', 1)
line(xl, xl*200, 'LineWidth', 1)
legend('data', 'c:chla = 50', 'c:chla = 200')
hold on
clear req szn

subplot(2, 2, 4)
szn = 'Fal';
req = (seasons.season == szn & ~isoutlier(cchl100.chl_0_avg, 'median'));
scatter(cchl100.chl_0_avg(req), cchl100.sumC_total(req), 40, cchl100.depth_m(req), 'filled', 'o');
colormap(flipud(parula))
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
set(gca, 'FontSize', 16)
title('Fall')
ylim([0 300])
xlim([0 5])
xl = xlim;
line(xl, xl*50, 'LineWidth', 1)
line(xl, xl*200, 'LineWidth', 1)
legend('data', 'c:chla = 50', 'c:chla = 200')
hold on
clear req szn

c = colorbar;
c.Direction = 'reverse';
colormap(flipud(parula))
c.Label.String = 'depth (m)';
sgtitle('C vs. chl by season for all total samples, ~isoutlier median')

%% 7/27 BOX and SCATTER WITH ALL MONTHS REPRESENTED (ADD EMPTY MONTHS AS NANS, CONC. ARRAYS)
% overview plot, ≤100 m
monthNames = {'Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'};

concreq = [cchl100.C_chl_ratio_total; nan(12, 1)];
concmon = [seasons.monthcat; monthNames];
concmon = categorical(concmon);
concdpth = [cchl100.depth_m; nan(12, 1)];

ax1 = axes();
scatter(ax1, concmon, concreq, 40, concdpth, "filled", 'o', 'jitter', 'on', 'jitteramount', 0.1, 'Clipping', 'on')
ax1.XAxis.Categories = {'Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'};
set(ax1, 'YLim', [-20 510])
pos = get(ax1, 'Position');
colormap(flipud(parula))
c = colorbar;
c.Label.String = 'Depth (m)';
c.Location = 'eastoutside';
c.Direction = 'reverse';
set(ax1, 'Position', pos)
hold on

ax2 = axes('Position', get(ax1, 'Position'), 'YLim', get(ax1, 'YLim'),  'PlotBoxAspectRatio', get(ax1, 'PlotBoxAspectRatio'), 'Visible','off', 'Color', 'none');
ax2.XLim = [min(double(seasons.monthcat)) - 0.5, max(double(seasons.monthcat)) + 0.5];
ax2.XTick = 1:numel(categories(seasons.monthcat));
hold(ax2, 'on')
boxplot(ax2, concreq, concmon, 'GroupOrder', monthNames, 'OutlierSize', 10,'DataLim', [0 500],'ExtremeMode', 'clip')
ax2.YLim = get(ax1, 'YLim');
dashedLines = findobj(ax2, 'Type', 'Line', 'Tag', '','LineStyle', '--');
delete(dashedLines)
title(ax1, 'C:chl a by month, totals ≤100 m depth', 'FontWeight', 'bold')
set(ax1, 'FontSize', 16)

xt = 1:numel(monthNames);
for i = 1:length(monthNames)
    thisMonth = monthNames{i};
    % logical indexing to count how many match
    n = sum(concmon == thisMonth & ~isnan(concreq));
    % add text below each box
    if n > 0
    text(xt(i), -10, num2str(n),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 13)
    end
end

%% 7/27 subplots with all months represented!

monthPlaceholder = nan(1, 12);
monthNames = {'Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'};

% surface
subplot(2, 3, 1)

loc = 'innershelf'; % 'innershelf', 'midshelf', 'outershelf' are options
frac = 'total';
dotidx = append("C_chl_ratio_", frac);

req = (cchl100.depth_m <= 50 &  ~isnan(cchl100.(dotidx)) &  (distfromshore == loc));

concreq = [cchl100.(dotidx)(req); nan(12, 1)];
concmon = [seasons.month(req); monthNames];
concmon = categorical(concmon);

boxplot(concreq, concmon, 'GroupOrder', monthNames, 'DataLim', [0 350], 'ExtremeMode', 'compress')
xlabel('Month')
ylabel('C:chl a')
ylim([-10 350])
title(append(loc, '-', frac, '-0to50m'))

xt = 1:numel(monthNames);
for i = 1:length(monthNames)
    thisMonth = monthNames{i};
    % logical indexing to count how many match
    n = sum(concmon == thisMonth & ~isnan(concreq));
    % add text below each box
    if n > 0
    text(xt(i), -10, num2str(n),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
    end
end
% 50-100 m
subplot(2, 3, 5)
loc = 'midshelf'; % 'innershelf', 'midshelf', 'outershelf' are options
frac = 'total';
dotidx = append("C_chl_ratio_", frac);

req = (cchl100.depth_m > 50 & cchl100.depth_m <= 100 & ~isnan(cchl100.(dotidx)) & (distfromshore == loc));

concreq = [cchl100.(dotidx)(req); nan(12, 1)];
concmon = [seasons.month(req); monthNames];
concmon = categorical(concmon);

boxplot(concreq, concmon, 'GroupOrder', monthNames)
xlabel('Month')
ylabel('C:chl a')
ylim([-10 350])
title(append(loc, '-', frac, '-50to100m'))

xt = 1:numel(monthNames);
for i = 1:length(monthNames)
    thisMonth = monthNames{i};
    % logical indexing to count how many match
    n = sum(concmon == thisMonth & ~isnan(concreq));
    % add text below each box
    if n > 0
    text(xt(i), -10, num2str(n),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
    end
end

%% or, plot all upper 100m by shelf distance

monthNames = {'Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'};

% surface
subplot(1, 3, 3)

loc = 'outershelf'; % 'innershelf', 'midshelf', 'outershelf' are options
frac = 'total';
dotidx = append("C_chl_ratio_", frac);

req = (~isnan(cchl100.(dotidx)) &  (distfromshore == loc));

concreq = [cchl100.(dotidx)(req); nan(12, 1)];
concmon = [seasons.month(req); monthNames];
concmon = categorical(concmon);

boxplot(concreq, concmon, 'GroupOrder', monthNames, 'DataLim', [0 350], 'ExtremeMode', 'compress')
xlabel('Month')
ylabel('C:chl a')
ylim([-10 350])
title(append(loc, '-', frac))

xt = 1:numel(monthNames);
for i = 1:length(monthNames)
    thisMonth = monthNames{i};
    % logical indexing to count how many match
    n = sum(concmon == thisMonth & ~isnan(concreq));
    % add text below each box
    if n > 0
    text(xt(i), -10, num2str(n),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
    end
end
%% plot each size fraction, median c:chl for each month by shelf distance

monthNames = {'Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'};

loc = 'outershelf'; % 'innershelf', 'midshelf', 'outershelf' are options
frac = 'greaterthan20';
dotidx = append("C_chl_ratio_", frac);

req = (~isnan(cchl100.(dotidx)) & (distfromshore == loc));

concreq = [cchl100.(dotidx)(req); nan(12, 1)];
concmon = [seasons.month(req); monthNames];
concreq = table(concreq);
concreq.mon = concmon;
med = table();
med.mon = monthNames;
med.med = zeros(length(med.mon), 1);
for i = 1:length(med.mon)
    tf = strcmp(concreq.mon, med.mon(i));
    med.med(i) = median(concreq.concreq(tf == 1), 'omitnan');
end
%for i = 1:length(med.med)
%    if med.med(i) < isoutlier(concreq.concreq(tf ==1), 'median') % exclude months with dramatic negatives from summary
%       med.med(i) = NaN;
%    end
%end

med.mon = categorical(med.mon);
med.mon = reordercats(med.mon, monthNames);
plot(med.mon, med.med, '-*', 'LineWidth', 1)
hold on
%%
%legend({'innershelf', 'midshelf', 'outershelf'}, 'FontSize', 22)
set(gca, 'FontSize', 14)
xlabel('Month')
ylabel('C:chl a')
ylim([0 400])
title(frac)

%% plotting to see anomalies in c vs. chl in size fractions

transectStations = {'MVCO', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9', 'L10', 'L11'};

idx = ismember(cchl100.nearest_station, transectStations);
cchltransect = cchl100(idx, :);
clear idx

subplot(2, 2, 1)
scatter(cchltransect.chl_lessthan5, cchltransect.sumC_5)
line(xlim, xlim*50)
ylim([0 150])
line(xlim, xlim*200)
legend('data', 'c:chl 50', 'c:chl 200')
title('lessthan5')

subplot(2, 2, 2)
scatter(cchltransect.chl_5to10, cchltransect.sumC_5to10)
line(xlim, xlim*50)
ylim([0 150])
line(xlim, xlim*200)
legend('data', 'c:chl 50', 'c:chl 200')
title('5to10')

subplot(2, 2, 3)
scatter(cchltransect.chl_10to20, cchltransect.sumC_10to20)
line(xlim, xlim*50)
ylim([0 150])
line(xlim, xlim*200)
legend('data', 'c:chl 50', 'c:chl 200')
title('10to20')

subplot(2, 2, 4)
scatter(cchltransect.chl_20_avg, cchltransect.sumC_greaterthan20)
line(xlim, xlim*50)
ylim([0 150])
line(xlim, xlim*200)
legend('data', 'c:chl 50', 'c:chl 200')
title('greaterthan20')

sgtitle('C vs. chl a for each non-total size fraction, only transect stations')

%% 7/26 looking at more recent chl data to see if there's a relationship between
% rb and chl a concentration (cruises EN695 and onward, accessed from API)

newercruises = {'EN695', 'EN706', 'EN712', 'EN715', 'EN720'};
% can't read in EN727 yet because chl data doesn't exist for it yet

CHL = table();
for i = 1:length(newercruises)
   cruise = string(newercruises(i));
   CHL2 = readtable(append('https://nes-lter-data.whoi.edu/api/chl/', cruise, '.csv'));
   CHL = [CHL; CHL2];
   clear CHL2
end

idx = ismember(CHL(:, {'cruise', 'cast', 'niskin'}), cchltransect(:, {'cruise','cast','niskin'}));
CHL = CHL(idx == 1, :);

tf = strcmp(CHL.filter_size, '>0');
totals = CHL(tf == 1, :);

scatter(totals.chl, totals.rb, 40, totals.depth, 'filled')
xlabel('chl a (µg/L)')
ylabel('rb')
title('totals, EN695 and later')
c = colorbar;
colormap(flipud(parula))
c.Direction = 'reverse';
mdl = fitlm(totals.chl, totals.rb, "linear");
plot(mdl)

%% rb subplots

subplot(2, 3, 1)
cruise = 'EN695';
tf = strcmp(string(totals.cruise), cruise);
scatter(totals.chl(tf == 1), totals.rb(tf == 1), 40, totals.ratio(tf == 1), 'filled')
xlabel('chl a (µg/L)')
ylabel('rb')
ylim([0 500])
title(append('totals, ', cruise))
c = colorbar;
c.Label.String = 'rb:chl a';
%mdl = fitlm(totals.chl(tf == 1), totals.rb(tf == 1), 'linear', Exclude=(totals.flag(tf == 1)));
%plot(mdl)
clear tf cruise

subplot(2, 3, 2)
cruise = 'EN706';
tf = strcmp(string(totals.cruise), cruise);
scatter(totals.chl(tf == 1), totals.rb(tf == 1), 40, totals.ratio(tf == 1), 'filled')
xlabel('chl a (µg/L)')
ylabel('rb')
ylim([0 500])
title(append('totals, ', cruise))
c = colorbar;
c.Label.String = 'rb:chl a';
clear tf cruise

subplot(2, 3, 3)
cruise = 'EN712';
tf = strcmp(string(totals.cruise), cruise);
scatter(totals.chl(tf == 1), totals.rb(tf == 1), 40, totals.ratio(tf == 1), 'filled')
xlabel('chl a (µg/L)')
ylabel('rb')
ylim([0 500])
title(append('totals, ', cruise))
c = colorbar;
c.Label.String = 'rb:chl a';
clear tf cruise

subplot(2, 3, 4)
cruise = 'EN715';
tf = strcmp(string(totals.cruise), cruise);
scatter(totals.chl(tf == 1), totals.rb(tf == 1), 40, totals.ratio(tf == 1), 'filled')
xlabel('chl a (µg/L)')
ylabel('rb')
ylim([0 500])
title(append('totals, ', cruise))
c = colorbar;
c.Label.String = 'rb:chl a';
clear tf cruise

subplot(2, 3, 5)
cruise = 'EN720';
tf = strcmp(string(totals.cruise), cruise);
scatter(totals.chl(tf == 1), totals.rb(tf == 1), 40, totals.ratio(tf == 1), 'filled')
xlabel('chl a (µg/L)')
ylabel('rb')
ylim([0 500])
title(append('totals, ', cruise))
c = colorbar;
c.Label.String = 'rb:chl a';
clear tf cruise

% exclusion criterion: any value for totals.ratio < 100
excluded = totals(totals.ratio < 100, :);
included = totals(totals.ratio > 100, :);

% histogram of excluded values
histogram(excluded.rb)
xlabel('Chl a (µg/L)')
ylabel('count')
title('excluded values based on rb:chl a < 100')
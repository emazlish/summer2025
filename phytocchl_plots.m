% new plotting script for new phytocchl table 7/24/25 - EM

% assigning each row a season designation
m = datetime(phytocchl.date_sampled, 'InputFormat','yyyy-MM-dd HH:mm:ss+00:00');

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

distfromshore = strings(size(phytocchl.nearest_station, 1), 1);
for i = 1:height(distfromshore)
    if ismember(phytocchl.nearest_station(i), {'MVCO','L1', 'L2'}) == 1
        distfromshore(i) = 'innershelf';
    end
    if ismember(phytocchl.nearest_station(i), {'L3', 'L4','L5','L6', 'L7', 'L8', 'L9'}) == 1
        distfromshore(i) = 'midshelf';
    end
    if ismember(phytocchl.nearest_station(i), {'L10', 'L11'}) == 1
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

%% plotting by shelf distance and depth, c vs. chl by season

subplot(1, 3, 3)
loc = 'outershelf'; %'innershelf', 'midshelf', 'outershelf' are options
chlfrac = '0_avg'; % options are 0_avg, greaterthan20, lessthan5, 5to10, 10to20
chllabel = append('chl_', chlfrac);
cfrac = 'total'; % options are total, 5, 10, 5to10, 10to20, greaterthan20
sumClabel = append('sumC_', cfrac);
scatter(phytocchl.(chllabel)(seasons.season == 'Win' & distfromshore == loc & phytocchl.depth_m <= 50), phytocchl.(sumClabel)(seasons.season == 'Win' & distfromshore == loc & phytocchl.depth_m <= 50), 40, phytocchl.depth_m(seasons.season == 'Win' & distfromshore == loc & phytocchl.depth_m <= 50), 'filled', 'o')
hold on
scatter(phytocchl.(chllabel)(seasons.season == 'Sum' & distfromshore == loc & phytocchl.depth_m <= 50), phytocchl.(sumClabel)(seasons.season == 'Sum' & distfromshore == loc & phytocchl.depth_m <= 50), 40, phytocchl.depth_m(seasons.season == 'Sum' & distfromshore == loc & phytocchl.depth_m <= 50), 'filled', '^')
hold on
scatter(phytocchl.(chllabel)(seasons.season == 'Fal' & distfromshore == loc & phytocchl.depth_m <= 50), phytocchl.(sumClabel)(seasons.season == 'Fal' & distfromshore == loc & phytocchl.depth_m <= 50), 40, phytocchl.depth_m(seasons.season == 'Fal' & distfromshore == loc & phytocchl.depth_m <= 50), 'filled', 'p')
colormap(flipud(parula))

legend({'Winter', 'Summer', 'Fall'}, "Location","northeast")
set(gca, "FontSize", 16)
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
title(loc)

% do these at the end
%c = colorbar;
%c.Direction = 'reverse';
%c.Label.String = 'depth (m)';
%sgtitle('C and chl for upper 50 m of ocean across shelf')
%clim([0 50])

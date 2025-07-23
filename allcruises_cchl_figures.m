% plotting all cruises cchl from the allcruises_cchl file
% last modified by EM 20250718 

%% initial plotting attempts

% remove outliers outside 3 st dev. from the mean
[~, d] = rmoutliers(allcruises_cchl.C_chl_ratio_total, "mean");

B = allcruises_cchl(d == 0, :);
outliers = allcruises_cchl(d == 1, :);

% totals
scatter(B.latitude, B.depth_m, 50, B.C_chl_ratio_total,'filled')
xlabel('Latitude (˚N)')
ylabel('Depth (m)')
set(gca, 'Ydir', 'reverse')
set(gca, 'xdir', 'reverse')
c = colorbar('eastoutside');
c.Label.String = 'Carbon:chlorophyll a';
title('C:chl a')
c.Label.FontSize = 20;
set(gca, 'FontSize', 20);

% C vs. chl
scatter(B.chl_0_avg, B.sumC_total, 40, B.depth_m, 'filled')
c = colorbar('eastoutside');
c.Direction = 'reverse';
%clim([0 100])
colormap(flipud(parula)); % makes it so surf is warm and depth is cold
c.Label.String = 'Depth (m)';
title('C:chl a')
set(gca, 'FontSize', 16);
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
xl = xlim;
yl = ylim;
med = median(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total)));
stdev = std(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total)));
%line(xl, xl*50)
hold on
line(xl, (xl*med), 'Color', 'red')
ylim(yl)
hold on
line(xl, xl*(med-stdev), 'Color', 'red', 'LineStyle', '--' )
ylim(yl)
hold on
line(xl, xl*(med+stdev), 'Color', 'red', 'LineStyle', '--' )

% C vs. chl for upper 100 m of ocean
T = B(B.depth_m < 100, :);
scatter(T.chl_0_avg, T.sumC_total, 40, T.depth_m, 'filled')
c = colorbar('eastoutside');
c.Direction = 'reverse';
%clim([0 100])
colormap(flipud(parula)); % makes it so surf is warm and depth is cold
c.Label.String = 'Depth (m)';
title('C:chl a')
set(gca, 'FontSize', 16);
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
xl = xlim;
yl = ylim;
med = median(T.C_chl_ratio_total(~isnan(T.C_chl_ratio_total)));
stdev = std(T.C_chl_ratio_total(~isnan(T.C_chl_ratio_total)));
%line(xl, xl*50)
hold on
line(xl, (xl*med), 'Color', 'red')
ylim(yl)
hold on
line(xl, xl*(med-stdev), 'Color', 'red', 'LineStyle', '--' )
ylim(yl)
hold on
line(xl, xl*(med+stdev), 'Color', 'red', 'LineStyle', '--' )


B.date_sampled = datetime(B.date_sampled, 'InputFormat', 'yyyy-MM-dd HH:mm:ss+00:00');
B.month = month(B.date_sampled, 'shortname');

boxplot(B.C_chl_ratio_total, B.month, 'GroupOrder', {'Jan','Feb', 'May', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
xlabel('Month')
ylabel('C:chl a (µg/L)')
set(gca, 'FontSize', 14)

%% assigning each row a season designation - 7/18/25
m = datetime(allcruises_cchl.date_sampled, 'InputFormat','yyyy-MM-dd HH:mm:ss+00:00');

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
%% add season info to table and pull only values with data - 7/18/25
% 7/22 - need to regenerate B before plotting more because EN715 updated
allcruises_cchl.sumC_total = real(allcruises_cchl.sumC_total); % deal with any complex numbers
allcruises_cchl.C_chl_ratio_total = real(allcruises_cchl.C_chl_ratio_total); % deal with any complex numbers

chlidx = ~isnan(allcruises_cchl.chl_0_avg); % find all chl data that are not nan
cidx = ~isnan(allcruises_cchl.sumC_total); % find all c data that are not nan

[~, d] = rmoutliers(allcruises_cchl.sumC_total(cidx), "median"); % ID which rows have C >3 stdev from mean
[~, e] = rmoutliers(allcruises_cchl.chl_0_avg(chlidx), "median");

bad = false(size(allcruises_cchl.C_chl_ratio_total, 1), 1); % create logical array where 0 = good data
cbad(cidx) = d; % outlier non-nan carbon values are ID'd
cbad = cbad'; % might not be necessary
chlbad(chlidx) = e; % outlier non-nan chl values are ID'd
chlbad = chlbad'; % might not be necessary

% look through non-nan data to find rows with C:chl >3 stdev from mean
bad = cbad | chlbad | isnan(allcruises_cchl.sumC_total) | isnan(allcruises_cchl.chl_0_avg);

B = allcruises_cchl(bad == 0, :); % create array of all cruises that have C:chl ≠ nan
B.season = seasons.season(bad == 0); % grab season IDs for non nan/outliers

outliers = allcruises_cchl(cbad == 1 | chlbad == 1, :); % create array of outliers
nans = isnan(allcruises_cchl.sumC_total) | isnan(allcruises_cchl.chl_0_avg); % index all entries with nans for either c or chl
nantable = allcruises_cchl(nans == 1, :); % create a table containing all nan entries

clear d e cidx cbad chlidx chlbad nans bad
%% plot by season, size fraction, 0-50m, median + MAD
% C vs. chl

% winter
subplot(1, 3, 1)
scatter(B.chl_0_avg(B.season == 'Win' & B.depth_m <= 50), B.sumC_total(B.season == 'Win' & B.depth_m <= 50), 40, B.depth_m(B.season == 'Win' & B.depth_m <= 50), 'filled', 'o');
colormap(flipud(parula))
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
set(gca, 'FontSize', 16)
title('Winter')
ylim([0 200])
hold on

% add median c:chl ratio line and 1 MAD c:chl ratio lines

xl = xlim;
yl = ylim;
med = median(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total(B.season == 'Win' & B.depth_m <= 50))));
mad = mad(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total(B.season == 'Win' & B.depth_m <= 50))), 1);
hold on
line(xl, (xl*med), 'Color', 'red')
ylim(yl)
hold on
line(xl, xl*(med-mad), 'Color', 'red', 'LineStyle', '--' )
ylim(yl)
hold on
line(xl, xl*(med+mad), 'Color', 'red', 'LineStyle', '--' )
clear xl yl mad med
hold on

% summer
subplot(1, 3, 2)
scatter(B.chl_0_avg(B.season == 'Sum' & B.depth_m <= 50), B.sumC_total(B.season == 'Sum' & B.depth_m <= 50), 40, B.depth_m(B.season == 'Sum' & B.depth_m <= 50), 'filled', '^');
colormap(flipud(parula))
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
set(gca, 'FontSize', 16)
title('Summer')

% add median c:chl ratio line and 1 stdev c:chl ratio lines
xl = xlim;
yl = ylim;
med = median(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total(B.season == 'Sum' & B.depth_m <= 50))));
mad = mad(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total(B.season == 'Sum' & B.depth_m <= 50))), 1);
hold on
line(xl, (xl*med), 'Color', 'red')
ylim(yl)
hold on
line(xl, xl*(med-mad), 'Color', 'red', 'LineStyle', '--' )
ylim(yl)
hold on
line(xl, xl*(med+mad), 'Color', 'red', 'LineStyle', '--' )
clear xl yl mad med
hold on

% fall
subplot(1, 3, 3)
scatter(B.chl_0_avg(B.season == 'Fal' & B.depth_m <= 50), B.sumC_total(B.season == 'Fal' & B.depth_m <= 50), 50, B.depth_m(B.season == 'Fal' & B.depth_m <= 50), 'filled', 'p');
c = colorbar('eastoutside');
colormap(flipud(parula))
c.Direction = 'reverse';
c.Label.String = 'Depth (m)';
hold on
set(gca, 'FontSize', 16)
title('Fall')
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')

% add median c:chl ratio line and 1 stdev c:chl ratio lines

xl = xlim;
yl = ylim;
med = median(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total(B.season == 'Fal' & B.depth_m <= 50))));
mad = mad(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total(B.season == 'Fal' & B.depth_m <= 50))), 1);
hold on
line(xl, (xl*med), 'Color', 'red')
ylim(yl)
hold on
line(xl, xl*(med-mad), 'Color', 'red', 'LineStyle', '--' )
ylim(yl)
hold on
line(xl, xl*(med+mad), 'Color', 'red', 'LineStyle', '--' )
clear xl yl mad med
hold on

% entire subplot grid
sgtitle('C:chl a for 0-50m by season, median + MAD', 'FontWeight', 'Bold');

%% plot by season, total size fraction, all depths, median + MAD

% C vs. chl
% winter
subplot(1, 3, 1)
scatter(B.chl_0_avg(B.season == 'Win'), B.sumC_total(B.season == 'Win'), 40, B.depth_m(B.season == 'Win'), 'filled', 'o');
colormap(flipud(parula))
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
set(gca, 'FontSize', 16)
title('Winter')
hold on

% add median c:chl ratio line and 1 stdev c:chl ratio lines

xl = xlim;
yl = ylim;
med = median(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total(B.season == 'Win'))));
mad = mad(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total(B.season == 'Win'))), 1);
hold on
line(xl, (xl*med), 'Color', 'red')
ylim(yl)
hold on
line(xl, xl*(med-mad), 'Color', 'red', 'LineStyle', '--' )
ylim(yl)
hold on
line(xl, xl*(med+mad), 'Color', 'red', 'LineStyle', '--' )
clear xl yl mad med
hold on

% summer
subplot(1, 3, 2)
scatter(B.chl_0_avg(B.season == 'Sum'), B.sumC_total(B.season == 'Sum'), 40, B.depth_m(B.season == 'Sum'), 'filled', '^');
colormap(flipud(parula))
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
set(gca, 'FontSize', 16)
title('Summer')

% add median c:chl ratio line and 1 stdev c:chl ratio lines
xl = xlim;
yl = ylim;
med = median(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total(B.season == 'Sum'))));
mad = mad(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total(B.season == 'Sum'))), 1);
hold on
line(xl, (xl*med), 'Color', 'red')
ylim(yl)
hold on
line(xl, xl*(med-mad), 'Color', 'red', 'LineStyle', '--' )
ylim(yl)
hold on
line(xl, xl*(med+mad), 'Color', 'red', 'LineStyle', '--' )
clear xl yl mad med
hold on

% fall
subplot(1, 3, 3)
scatter(B.chl_0_avg(B.season == 'Fal'), B.sumC_total(B.season == 'Fal'), 40, B.depth_m(B.season == 'Fal'), 'filled', 'p');
c = colorbar('eastoutside');
colormap(flipud(parula))
c.Direction = 'reverse';
c.Label.String = 'Depth (m)';
hold on
set(gca, 'FontSize', 16)
title('Fall')
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')

% add median c:chl ratio line and 1 stdev c:chl ratio lines

xl = xlim;
yl = ylim;
med = median(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total(B.season == 'Fal'))));
mad = mad(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total(B.season == 'Fal'))), 1);
hold on
line(xl, (xl*med), 'Color', 'red')
ylim(yl)
hold on
line(xl, xl*(med-mad), 'Color', 'red', 'LineStyle', '--' )
ylim(yl)
hold on
line(xl, xl*(med+mad), 'Color', 'red', 'LineStyle', '--' )
clear xl yl mad med
hold on

% entire subplot grid
sgtitle('C:chl a for all depths by season for NES-LTER cruises 2018–2024', 'FontWeight', 'Bold');


%% seasons plotted on top of each other as different shapes - 7/18/25


scatter(B.chl_0_avg(B.season == 'Win'), B.sumC_total(B.season == 'Win'), 40, B.depth_m(B.season == 'Win'), 'filled', 'o', 'MarkerEdgeColor', 'k');
colormap(flipud(parula))

hold on
scatter(B.chl_0_avg(B.season == 'Sum'), B.sumC_total(B.season == 'Sum'), 40, B.depth_m(B.season == 'Sum'), 'filled', '^', 'MarkerEdgeColor', 'k');
colormap(flipud(parula))

hold on
scatter(B.chl_0_avg(B.season == 'Fal'), B.sumC_total(B.season == 'Fal'), 40, B.depth_m(B.season == 'Fal'), 'filled', 'p', 'MarkerEdgeColor', 'k');
colormap(flipud(parula))

% general plotting stuff
c = colorbar('eastoutside');
c.Direction = 'reverse';
legend({'Winter', 'Summer', 'Fall'}, "Location","northeast")
set(gca, "FontSize", 16)
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
title("C:chl a by season for all NES-LTER cruises 2018-2024", 'FontWeight', 'bold')

% plotting median and stdev lines
xl = xlim;
yl = ylim;
med = median(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total)));
stdev = std(B.C_chl_ratio_total(~isnan(B.C_chl_ratio_total)));
hold on
line(xl, (xl*med), 'Color', 'red')
ylim(yl)
hold on
line(xl, xl*(med-stdev), 'Color', 'red', 'LineStyle', '--' )
ylim(yl)
hold on
line(xl, xl*(med+stdev), 'Color', 'red', 'LineStyle', '--' )
clear xl yl stdev med

%% plotting only upper 100m of samples to see variation better

T = B(B.depth_m < 100, :);

scatter(T.chl_0_avg(T.season == 'Win'), T.sumC_total(T.season == 'Win'), 40, T.depth_m(T.season == 'Win'), 'filled', 'o')
hold on
scatter(T.chl_0_avg(T.season == 'Sum'), T.sumC_total(T.season == 'Sum'), 40, T.depth_m(T.season == 'Sum'), 'filled', '^')
hold on
scatter(T.chl_0_avg(T.season == 'Fal'), T.sumC_total(T.season == 'Fal'), 40, T.depth_m(T.season == 'Fal'), 'filled', 'p')
c = colorbar;
colormap(flipud(parula))
c.Direction = 'reverse';
legend({'Winter', 'Summer', 'Fall'}, "Location","northeast")
set(gca, "FontSize", 16)
xlabel('Chl a (µg L^{-1})')
ylabel('Carbon (µg L^{-1})')
title("C:chl a by season & upper 100m of ocean, for NES-LTER cruises 2018-2024", 'FontWeight', 'bold')

%% all samples - total fraction, boxplot of month vs. c:chl ratio

[~, d] = rmoutliers(allcruises_cchl.C_chl_ratio_total, "mean");

P = allcruises_cchl(d == 0, :);
cchloutliers = allcruises_cchl(d == 1, :);

P.date_sampled = datetime(P.date_sampled, 'InputFormat', 'yyyy-MM-dd HH:mm:ss+00:00');
P.month = month(P.date_sampled, 'shortname');
[q, ~, uniqueid] = unique(P.month);
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';

%top of box plot
hold on
%P.month = char(P.month);
boxplot(P.C_chl_ratio_total, P.monthcat,'GroupOrder', {'Jan','Feb', 'May','Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
xlabel('Month')
ylabel('C:chl a')
set(gca, 'FontSize', 14, 'LineWidth', 1)
title("C:chl a by month for NES-LTER cruises 2018-2024", 'FontWeight', 'bold')

neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);

hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -10, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 14)
end

clear d labels qsorted q uniqueid counts xt P neworder

%% **doesn't work to plot scatter and boxplot data on top of each other
% if desired, plot each month's data on top of box plot and color by depth

hold on
P.monthcat = categorical(P.month, {'Jan','Feb', 'May','Jul', 'Aug', 'Sep', 'Oct', 'Nov'}, 'Ordinal', 1);
scatter(P.monthcat, P.C_chl_ratio_total, 40, P.depth_m, 'filled',  'o', 'jitter', 'on', 'jitteramount', 0.15) % can scatter actual data points on
colormap(flipud(parula))
c = colorbar;
c.Direction = 'reverse';
xlabel('Month')
ylabel('C:chl a')
set(gca, 'FontSize', 14, 'LineWidth', 1)
title("C:chl a total fraction by month for NES-LTER cruises 2018-2024", 'FontWeight', 'bold')


%% all samples - total fraction, violin plot (not super useful)

x = reordercats(P.month, {'Jan','Feb', 'May','Jul', 'Aug', 'Sep', 'Oct', 'Nov'});
violinplot(x, P.C_chl_ratio_total);

clear A x

%% all samples - any fraction desired, subplots by depth *OUTLIERS INCLUDED

% 0-50 m
subplot(1, 3, 1)
boxplot(P.C_chl_ratio_greaterthan20(P.depth_m <= 50 & ~isnan(P.C_chl_ratio_greaterthan20)), ...
    P.monthcat (P.depth_m <= 50 & ~isnan(P.C_chl_ratio_greaterthan20)), ...
    'GroupOrder', {'Jan','Feb', 'May','Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
xlabel('Month')
ylabel('C:chl a')
title('0-50 m')
%xlim([0 Inf])
[q, ~, uniqueid] = unique(P.month(P.depth_m <= 50 & ~isnan(P.C_chl_ratio_greaterthan20)));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), 2250, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt


%50-100 m

subplot(1, 3, 2)
boxplot(P.C_chl_ratio_greaterthan20(P.depth_m > 50 & P.depth_m <= 100 & ...
    ~isnan(P.C_chl_ratio_greaterthan20)), P.monthcat (P.depth_m > 50 & P.depth_m ...
    <= 100 & ~isnan(P.C_chl_ratio_greaterthan20)),'GroupOrder', ...
    {'Jan','Feb', 'May','Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
xlabel('Month')
ylabel('C:chl a')
title('50-100 m')
[q, ~, uniqueid] = unique(P.month(P.depth_m > 50 & P.depth_m <= 100 & ...
    ~isnan(P.C_chl_ratio_greaterthan20)));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), 1000, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt

% > 100 m

subplot(1, 3, 3)
%over100nan = P.C_chl_ratio_total(P.depth_m > 100)
boxplot(P.C_chl_ratio_greaterthan20(P.depth_m > 100 & ~isnan(P.depth_m) & ...
    ~isnan(P.C_chl_ratio_greaterthan20)), P.monthcat (P.depth_m > 100 & ...
    ~isnan(P.depth_m) & ~isnan(P.C_chl_ratio_greaterthan20)),'GroupOrder', {'Feb', 'Aug', 'Sep'})
xlabel('Month')
ylabel('C:chl a')
title('>100 m')

[q, ~, uniqueid] = unique(P.month(P.depth_m > 100 & ~isnan(P.depth_m) & ...
    ~isnan(P.C_chl_ratio_greaterthan20))); % this should exclude nans
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [2 1 3]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), 50, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt

%% same boxplots as above, but OUTLIERS REMOVED (~isoutlier, 'mean')

% 0-50 m
% & ~isoutlier(P.C_chl_ratio_lessthan5, "mean")
subplot(1, 3, 1)
boxplot(P.C_chl_ratio_greaterthan20(P.depth_m <= 50 & ~isnan(P.C_chl_ratio_greaterthan20) ...
    & ~isoutlier(P.C_chl_ratio_greaterthan20, "median")), P.monthcat (P.depth_m <= 50 ...
    & ~isnan(P.C_chl_ratio_greaterthan20) & ~isoutlier(P.C_chl_ratio_greaterthan20, "median")), ...
    'GroupOrder', {'Jan','Feb', 'May','Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
xlabel('Month')
ylabel('C:chl a')
title('0-50 m')
%xlim([0 Inf])
[q, ~, uniqueid] = unique(P.month(P.depth_m <= 50 & ~isnan(P.C_chl_ratio_greaterthan20) ...
    & ~isoutlier(P.C_chl_ratio_greaterthan20, "median")));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -10, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt


%50-100 m

subplot(1, 3, 2)
boxplot(P.C_chl_ratio_greaterthan20(P.depth_m > 50 & P.depth_m <= 100 & ...
    ~isnan(P.C_chl_ratio_greaterthan20) & ~isoutlier(P.C_chl_ratio_greaterthan20, "median")), ...
    P.monthcat (P.depth_m > 50 & P.depth_m <= 100 & ~isnan(P.C_chl_ratio_greaterthan20) ...
    & ~isoutlier(P.C_chl_ratio_greaterthan20, "median")),'GroupOrder', ...
    {'Jan','Feb', 'May', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
xlabel('Month')
ylabel('C:chl a')
title('50-100 m')
[q, ~, uniqueid] = unique(P.month(P.depth_m > 50 & P.depth_m <= 100 & ...
    ~isnan(P.C_chl_ratio_greaterthan20) & ~isoutlier(P.C_chl_ratio_greaterthan20, "median")));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -10, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt

% > 100 m

subplot(1, 3, 3)
%over100nan = P.C_chl_ratio_total(P.depth_m > 100)
boxplot(P.C_chl_ratio_greaterthan20(P.depth_m > 100 & ~isnan(P.depth_m) & ...
    ~isnan(P.C_chl_ratio_greaterthan20) & ~isoutlier(P.C_chl_ratio_greaterthan20, "median")), ...
    P.monthcat (P.depth_m > 100 & ~isnan(P.depth_m) & ~isnan(P.C_chl_ratio_greaterthan20) ...
    & ~isoutlier(P.C_chl_ratio_greaterthan20, "median")),'GroupOrder', {'Feb', 'Aug'})
xlabel('Month')
ylabel('C:chl a')
title('>100 m')

[q, ~, uniqueid] = unique(P.month(P.depth_m > 100 & ~isnan(P.depth_m) ...
    & ~isnan(P.C_chl_ratio_greaterthan20) & ~isoutlier(P.C_chl_ratio_greaterthan20, "median"))); % this should exclude nans
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [2 1]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -2, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt

%% subplot upper 50m by size fraction, exclude outliers using 'median'
% this section uses dynamic field names, so you don't have to change the
% variable being dot indexed in the table in every single location

% total fraction
subplot(2, 3, 1)

frac = 'total';
dotidx = append("C_chl_ratio_", frac);

boxplot(P.(dotidx)(P.depth_m <= 50 & ~isnan(P.(dotidx)) ...
    & ~isoutlier(P.(dotidx), "median")), P.monthcat (P.depth_m <= 50 ...
    & ~isnan(P.(dotidx)) & ~isoutlier(P.(dotidx), "median")), ...
    'GroupOrder', {'Jan','Feb', 'May','Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
xlabel('Month')
ylabel('C:chl a')
title(frac)
%xlim([0 Inf])
[q, ~, uniqueid] = unique(P.month(P.depth_m <= 50 & ~isnan(P.(dotidx)) ...
    & ~isoutlier(P.(dotidx), "median")));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -10, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt frac dotidx

% <5 fraction
subplot(2, 3, 2)

frac = 'lessthan5';
dotidx = append("C_chl_ratio_", frac);

boxplot(P.(dotidx)(P.depth_m <= 50 & ~isnan(P.(dotidx)) ...
    & ~isoutlier(P.(dotidx), "median")), P.monthcat (P.depth_m <= 50 ...
    & ~isnan(P.(dotidx)) & ~isoutlier(P.(dotidx), "median")), ...
    'GroupOrder', {'Jan','Feb', 'May','Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
xlabel('Month')
ylabel('C:chl a')
title(append(frac, ' µm'))
%xlim([0 Inf])
[q, ~, uniqueid] = unique(P.month(P.depth_m <= 50 & ~isnan(P.(dotidx)) ...
    & ~isoutlier(P.(dotidx), "median")));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -50, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt frac dotidx

% 5-10 fraction
subplot(2, 3, 3)

frac = '5to10';
dotidx = append("C_chl_ratio_", frac);

boxplot(P.(dotidx)(P.depth_m <= 50 & ~isnan(P.(dotidx)) ...
    & ~isoutlier(P.(dotidx), "median")), P.monthcat (P.depth_m <= 50 ...
    & ~isnan(P.(dotidx)) & ~isoutlier(P.(dotidx), "median")), ...
    'GroupOrder', {'Jan','Feb', 'May','Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
xlabel('Month')
ylabel('C:chl a')
title(append(frac, ' µm'))
%xlim([0 Inf])
[q, ~, uniqueid] = unique(P.month(P.depth_m <= 50 & ~isnan(P.(dotidx)) ...
    & ~isoutlier(P.(dotidx), "median")));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -80, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt frac dotidx

% 10-20 fraction
subplot(2, 3, 4)

frac = '10to20';
dotidx = append("C_chl_ratio_", frac);

boxplot(P.(dotidx)(P.depth_m <= 50 & ~isnan(P.(dotidx)) ...
    & ~isoutlier(P.(dotidx), "median")), P.monthcat (P.depth_m <= 50 ...
    & ~isnan(P.(dotidx)) & ~isoutlier(P.(dotidx), "median")), ...
    'GroupOrder', {'Jan','Feb', 'May','Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
xlabel('Month')
ylabel('C:chl a')
title(append(frac, ' µm'))
%xlim([0 Inf])
[q, ~, uniqueid] = unique(P.month(P.depth_m <= 50 & ~isnan(P.(dotidx)) ...
    & ~isoutlier(P.(dotidx), "median")));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -60, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt frac dotidx

% >20 fraction
subplot(2, 3, 5)
frac = 'greaterthan20';
dotidx = append("C_chl_ratio_", frac);

boxplot(P.(dotidx)(P.depth_m <= 50 & ~isnan(P.(dotidx)) ...
    & ~isoutlier(P.(dotidx), "median")), P.monthcat (P.depth_m <= 50 ...
    & ~isnan(P.(dotidx)) & ~isoutlier(P.(dotidx), "median")), ...
    'GroupOrder', {'Jan','Feb', 'May','Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
xlabel('Month')
ylabel('C:chl a')
title(append(frac, ' µm'))
%xlim([0 Inf])
[q, ~, uniqueid] = unique(P.month(P.depth_m <= 50 & ~isnan(P.(dotidx)) ...
    & ~isoutlier(P.(dotidx), "median")));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -10, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt frac dotidx

%% 50 - 100m, all sizes

subplot(2, 3, 5)
frac = 'greaterthan20';
dotidx = append("C_chl_ratio_", frac);

boxplot(P.(dotidx)(P.depth_m > 50 & P.depth_m <= 100 & ...
    ~isnan(P.(dotidx)) & ~isoutlier(P.(dotidx), "median")), ...
    P.monthcat (P.depth_m > 50 & P.depth_m <= 100 & ~isnan(P.(dotidx)) ...
    & ~isoutlier(P.(dotidx), "median")),'GroupOrder', ...
    {'Jan','Feb', 'May', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
xlabel('Month')
ylabel('C:chl a')
title(append(frac, ' µm'))
[q, ~, uniqueid] = unique(P.month(P.depth_m > 50 & P.depth_m <= 100 & ...
    ~isnan(P.(dotidx)) & ~isoutlier(P.(dotidx), "median")));
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [3 2 5 4 1 8 7 6]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -50, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt frac dotidx

%% >100 m depth, all sizes

subplot(2, 3, 2)
frac = 'lessthan5';
dotidx = append("C_chl_ratio_", frac);
boxplot(P.(dotidx)(P.depth_m > 100 & ~isnan(P.depth_m) & ...
    ~isnan(P.(dotidx)) & ~isoutlier(P.(dotidx), "median")), ...
    P.monthcat (P.depth_m > 100 & ~isnan(P.depth_m) & ~isnan(P.(dotidx)) ...
    & ~isoutlier(P.(dotidx), "median")),'GroupOrder', {'Feb', 'Aug'})
xlabel('Month')
ylabel('C:chl a')
title(append(frac, ' µm'))

[q, ~, uniqueid] = unique(P.month(P.depth_m > 100 & ~isnan(P.depth_m) ...
    & ~isnan(P.(dotidx)) & ~isoutlier(P.(dotidx), "median"))); % this should exclude nans
counts = histcounts(uniqueid);
counts = num2cell(counts);
q(:, 2) = counts';
neworder = [2 1];%[2 3 1 6 5 4]; % put months in order
qsorted = q(neworder, :);
hold on
labels = string(qsorted(:, 2));

xt = get(gca, 'Xtick');

for i = 1:length(labels)
    text(xt(i), -2, labels{i},'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10)
end
clear q counts uniqueid qsorted labels xt
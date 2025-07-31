% nitracline and euphotic depth things

%cruise = "EN715"; % this is the only thing you have to change to run script

% load nutrient data
tbl = readtable('/Users/emazlish/Downloads/nes_lter_surface_features_winter_summer_v20250116_bofu_added_emma.csv');

envivars = tbl(:, {'cruise', 'cast', 'station', 'date', 'latitude', 'longitude', 'nitracline_bofu', 'Euphotic_depth'});
cchl100.cruise = string(cchl100.cruise);

cruises = unique(cchl100.cruise);

% make cchl100 cruise names lowercase to match envivars case
for i = 1:length(cchl100.cruise)
    q = strcmp(cruises(:, 1), cchl100.cruise(i));
    if any(q)
        [~, idx] = ismember(cchl100.cruise(i), cruises(:, 1));
         c = cruises(idx, 2);
         cchl100.newcruise(i) = c;
     end
clear c idx q
end

% get rid of spaces in envivars cruise names
for i = 1:length(envivars.cruise)
    q = strcmp(ugh(:, 1), envivars.cruise(i));
    if any(q)
        [~, idx] = ismember(envivars.cruise(i), ugh(:, 1));
         c = ugh(idx, 2);
         envivars.newcruise(i) = c;
     end
clear c idx q
end

envivars.newcruise = string(envivars.newcruise);

% get row numbers from envivars that match each row of cchl100
[~, idx] = ismember(cchl100(:, {'newcruise', 'cast'}), envivars(:, {'newcruise', 'cast'}));

% grab nitracline info for each row
cchl100.nitracline = envivars.nitracline_bofu(idx);

% grab euphotic depth info for each row
cchl100.euphoticDepth = envivars.Euphotic_depth(idx);

% do z/z nitracline normalization
cchl100.normNitracline = cchl100.depth_m ./ cchl100.nitracline;

% do z/z euphotic normalization
cchl100.normEuphotic = cchl100.depth_m ./ cchl100.euphoticDepth;

%% 4 boxes, one per season. x = c:chl, y  = depth
subplot(2, 2, 1)
frac = 'C_chl_ratio_total';
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)));
scatter(cchl100.(frac)(req), cchl100.depth_m(req))
set(gca, 'YDir', 'reverse', 'FontSize', 14)
xlabel('C:chl a')
ylabel('Depth (m)')
title('Winter, totals')

subplot(2, 2, 2)
frac = 'C_chl_ratio_total';
req = (seasons.season == 'Spr' & ~isoutlier(cchl100.(frac)));
scatter(cchl100.(frac)(req), cchl100.depth_m(req))
set(gca, 'YDir', 'reverse', 'FontSize', 14)
xlabel('C:chl a')
ylabel('Depth (m)')
title('Spring, totals')

subplot(2, 2, 3)
frac = 'C_chl_ratio_total';
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)));
scatter(cchl100.(frac)(req), cchl100.depth_m(req))
set(gca, 'YDir', 'reverse', 'FontSize', 14)
xlabel('C:chl a')
ylabel('Depth (m)')
title('Summer, totals')

subplot(2, 2, 4)
frac = 'C_chl_ratio_total';
req = (seasons.season == 'Fal' & ~isoutlier(cchl100.(frac)));
scatter(cchl100.(frac)(req), cchl100.depth_m(req))
set(gca, 'YDir', 'reverse', 'FontSize', 14)
xlabel('C:chl a')
ylabel('Depth (m)')
title('Fall, totals')

%% 4 boxes, one per season. x = lat, y = depth, c = c:chl

frac = 'C_chl_ratio_total';
cmin = min(cchl100.(frac));
cmax = max(cchl100.(frac));

tiledlayout(2, 2)
nexttile;
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)));
scatter(cchl100.latitude(req), cchl100.depth_m(req), 40, cchl100.(frac)(req), 'filled')
clim = ([cmin cmax]);
set(gca, 'XDir', 'reverse', 'YDir', 'reverse', 'FontSize', 14)
xlabel('Latitude ˚N')
ylabel('Depth (m)')
title('Winter, totals')

nexttile;
req = (seasons.season == 'Spr' & ~isoutlier(cchl100.(frac)));
scatter(cchl100.latitude(req), cchl100.depth_m(req), 40, cchl100.(frac)(req), 'filled')
clim = ([cmin cmax]);
set(gca, 'XDir', 'reverse', 'YDir', 'reverse', 'FontSize', 14)
xlabel('Latitude ˚N')
ylabel('Depth (m)')
title('Spring, totals')

nexttile;
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)));
scatter(cchl100.latitude(req), cchl100.depth_m(req), 40, cchl100.(frac)(req), 'filled')
clim = ([cmin cmax]);
set(gca, 'XDir', 'reverse', 'YDir', 'reverse', 'FontSize', 14)
xlabel('Latitude ˚N')
ylabel('Depth (m)')
title('Summer, totals')

nexttile;
req = (seasons.season == 'Fal' & ~isoutlier(cchl100.(frac)));
scatter(cchl100.latitude(req), cchl100.depth_m(req), 40, cchl100.(frac)(req), 'filled')
clim = ([cmin cmax]);
set(gca, 'XDir', 'reverse', 'YDir', 'reverse', 'FontSize', 14)
xlabel('Latitude ˚N')
ylabel('Depth (m)')
title('Fall, totals')

c = colorbar;
c.Layout.Tile = 'east';
c.Direction = 'normal';
c.Label.String = 'C:chl a';

%% 2 boxes, summer/winter. x = c:chl, y = depth/euphotic depth
frac = 'C_chl_ratio_total';

tiledlayout(2, 1)
nexttile;
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)));
scatter(cchl100.(frac)(req), cchl100.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
title('Winter, totals')


nexttile;
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)));
scatter(cchl100.(frac)(req), cchl100.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
title('Summer, totals')

%% 2 boxes, summer/winter. x = c:chl, y = depth/depth of nitracline

frac = 'C_chl_ratio_total';
cmin = min(cchl100.(frac));
cmax = max(cchl100.(frac));

tiledlayout(2, 1)
nexttile;
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)));
scatter(cchl100.(frac)(req), cchl100.normNitracline(req), 'filled')
set(gca,  'FontSize', 14)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 4])
xlim([0 200])
title('Winter, totals')


nexttile;
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)));
scatter(cchl100.(frac)(req), cchl100.normNitracline(req), 'filled')
clim = ([cmin cmax]);
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 4])
xlim([0 200])

title('Summer, totals')

%% 6 boxes. season & distance from shore, x = cchl, y = z/znitracline

frac = 'C_chl_ratio_total';

tiledlayout(2, 3)
nexttile;
loc = 'innershelf';
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.normNitracline(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 3])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.normNitracline(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 3])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.normNitracline(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 3])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'innershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.normNitracline(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 3])
xlim([0 200])
title(append('Summer, totals-', loc))

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.normNitracline(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 3])
xlim([0 200])
title(append('Summer, totals-', loc))

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.normNitracline(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('Z/Znitracline')
mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
plot(mdl)
ylim([0 3])
xlim([0 200])
title(append('Summer, totals-', loc))

%% 6 boxes. season & distance from shore, x = cchl, y = z/zeuphotic

frac = 'C_chl_ratio_total';

tiledlayout(2, 3)
nexttile;
loc = 'innershelf';
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'innershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
xlim([0 200])
title(append('Summer, totals-', loc))

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
xlim([0 200])
title(append('Summer, totals-', loc))

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
xlim([0 200])
title(append('Summer, totals-', loc))

%% plotting by temperature

frac = 'C_chl_ratio_total';
tiledlayout(2, 3)
nexttile;
loc = 'innershelf';
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.potemp090c(req), 'filled')
set(gca, 'FontSize', 14)
%mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
%plot(mdl)
xlabel('C:chl a')
ylabel('temperature')
ylim([0 30])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.potemp090c(req), 'filled')
set(gca, 'FontSize', 14)
%mdl = fitlm(cchl100.(frac)(req), cchl100.normNitracline(req), 'linear');
%plot(mdl)
xlabel('C:chl a')
ylabel('temperature')
ylim([0 30])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Win' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.potemp090c(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('temperature')
ylim([0 30])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'innershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.potemp090c(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('temperature')
ylim([0 30])
xlim([0 200])
title(append('Summer, totals-', loc))

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.potemp090c(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('temperature')
ylim([0 30])
xlim([0 200])
title(append('Summer, totals-', loc))

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchl100.(frac)) & distfromshore == loc);
scatter(cchl100.(frac)(req), cchl100.potemp090c(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('temperature')
ylim([0 30])
xlim([0 200])
title(append('Summer, totals-', loc))
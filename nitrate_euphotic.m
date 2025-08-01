% nitracline and euphotic depth things

% load nutrient data
tbl = readtable('/Users/emazlish/Downloads/nes_lter_surface_features_winter_summer_v20250116_bofu_added_emma.csv');

envivars = tbl(:, {'cruise', 'cast', 'station', 'date', 'latitude', 'longitude', 'nitracline_bofu', 'Euphotic_depth'});
cchlQC.cruise = string(cchlQC.cruise);

cruises = unique(cchlQC.cruise);

% make cchl100 cruise names lowercase to match envivars case
for i = 1:length(cchlQC.cruise)
    q = strcmp(cruises(:, 1), cchlQC.cruise(i));
    if any(q)
        [~, idx] = ismember(cchlQC.cruise(i), cruises(:, 1));
         c = cruises(idx, 2);
         cchlQC.newcruise(i) = c;
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

% get row numbers from envivars that match each row of cchlQC
[~, idx] = ismember(cchlQC(:, {'newcruise', 'cast'}), envivars(:, {'newcruise', 'cast'}));

% grab nitracline info for each row
cchlQC.nitracline = envivars.nitracline_bofu(idx);

% grab euphotic depth info for each row
cchlQC.euphoticDepth = envivars.Euphotic_depth(idx);

% do z/z nitracline normalization
cchlQC.normNitracline = cchlQC.depth_m ./ cchlQC.nitracline;

% do z/z euphotic normalization
cchlQC.normEuphotic = cchlQC.depth_m ./ cchlQC.euphoticDepth;

% add mixed layer depth
cchlQC.mld = envivars.mld(idx);

%% 4 boxes, one per season. x = c:chl, y  = depth
subplot(2, 2, 1)
frac = 'C_chl_ratio_total';
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)));
scatter(cchlQC.(frac)(req), cchlQC.depth_m(req))
set(gca, 'YDir', 'reverse', 'FontSize', 14)
xlabel('C:chl a')
ylabel('Depth (m)')
title('Winter, totals')

subplot(2, 2, 2)
frac = 'C_chl_ratio_total';
req = (seasons.season == 'Spr' & ~isoutlier(cchlQC.(frac)));
scatter(cchlQC.(frac)(req), cchlQC.depth_m(req))
set(gca, 'YDir', 'reverse', 'FontSize', 14)
xlabel('C:chl a')
ylabel('Depth (m)')
title('Spring, totals')

subplot(2, 2, 3)
frac = 'C_chl_ratio_total';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)));
scatter(cchlQC.(frac)(req), cchlQC.depth_m(req))
set(gca, 'YDir', 'reverse', 'FontSize', 14)
xlabel('C:chl a')
ylabel('Depth (m)')
title('Summer, totals')

subplot(2, 2, 4)
frac = 'C_chl_ratio_total';
req = (seasons.season == 'Fal' & ~isoutlier(cchlQC.(frac)));
scatter(cchlQC.(frac)(req), cchlQC.depth_m(req))
set(gca, 'YDir', 'reverse', 'FontSize', 14)
xlabel('C:chl a')
ylabel('Depth (m)')
title('Fall, totals')

%% 4 boxes, one per season. x = lat, y = depth, c = c:chl

frac = 'C_chl_ratio_total';
cmin = min(cchlQC.(frac));
cmax = max(cchlQC.(frac));

tiledlayout(2, 2)
nexttile;
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)));
scatter(cchlQC.latitude(req), cchlQC.depth_m(req), 40, cchlQC.(frac)(req), 'filled')
clim = ([cmin cmax]);
set(gca, 'XDir', 'reverse', 'YDir', 'reverse', 'FontSize', 14)
xlabel('Latitude ˚N')
ylabel('Depth (m)')
title('Winter, totals')

nexttile;
req = (seasons.season == 'Spr' & ~isoutlier(cchlQC.(frac)));
scatter(cchlQC.latitude(req), cchlQC.depth_m(req), 40, cchlQC.(frac)(req), 'filled')
clim = ([cmin cmax]);
set(gca, 'XDir', 'reverse', 'YDir', 'reverse', 'FontSize', 14)
xlabel('Latitude ˚N')
ylabel('Depth (m)')
title('Spring, totals')

nexttile;
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)));
scatter(cchlQC.latitude(req), cchlQC.depth_m(req), 40, cchlQC.(frac)(req), 'filled')
clim = ([cmin cmax]);
set(gca, 'XDir', 'reverse', 'YDir', 'reverse', 'FontSize', 14)
xlabel('Latitude ˚N')
ylabel('Depth (m)')
title('Summer, totals')

nexttile;
req = (seasons.season == 'Fal' & ~isoutlier(cchlQC.(frac)));
scatter(cchlQC.latitude(req), cchlQC.depth_m(req), 40, cchlQC.(frac)(req), 'filled')
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
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)));
scatter(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
title('Winter, totals')


nexttile;
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)));
scatter(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
title('Summer, totals')

%% 2 boxes, summer/winter. x = c:chl, y = depth/depth of nitracline

frac = 'C_chl_ratio_total';
cmin = min(cchlQC.(frac));
cmax = max(cchlQC.(frac));

tiledlayout(2, 1)
nexttile;
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)));
scatter(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'filled')
set(gca,  'FontSize', 14)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 4])
xlim([0 200])
title('Winter, totals')


nexttile;
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)));
scatter(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'filled')
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
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 3])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 3])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 3])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'innershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 3])
xlim([0 200])
title(append('Summer, totals-', loc))

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Znitracline')
ylim([0 3])
xlim([0 200])
title(append('Summer, totals-', loc))

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('Z/Znitracline')
mdl = fitlm(cchlQC.(frac)(req), cchlQC.normNitracline(req), 'linear');
plot(mdl)
ylim([0 3])
xlim([0 200])
title(append('Summer, totals-', loc))

%% 6 boxes. season & distance from shore, x = cchl, y = z/zeuphotic

frac = 'C_chl_ratio_total';

tiledlayout(2, 3)
nexttile;
loc = 'innershelf';
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'innershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
xlim([0 200])
title(append('Summer, totals-', loc))

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('Z/Zeuphotic')
ylim([0 3])
xlim([0 200])
title(append('Summer, totals-', loc))

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.normEuphotic(req), 'linear');
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
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.potemp090c(req), 'filled')
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
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.potemp090c(req), 'filled')
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
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.potemp090c(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('temperature')
ylim([0 30])
xlim([0 200])
title(append('Winter, totals-', loc))

nexttile;
loc = 'innershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.potemp090c(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('temperature')
ylim([0 30])
xlim([0 200])
title(append('Summer, totals-', loc))

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.potemp090c(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('temperature')
ylim([0 30])
xlim([0 200])
title(append('Summer, totals-', loc))

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.potemp090c(req), 'filled')
set(gca, 'FontSize', 14)
xlabel('C:chl a')
ylabel('temperature')
ylim([0 30])
xlim([0 200])
title(append('Summer, totals-', loc))

%% plotting by mixed layer depth

% winter
frac = 'C_chl_ratio_total';
tiledlayout(1, 3)
nexttile;
loc = 'innershelf';
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.mld(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.mld(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('mld (m)')
ylim([0 200])
%xlim([0 200])
title(append('Winter, totals-', loc))
coeffs = mdl.Coefficients.Estimate;
eqnStr = sprintf('y = %.2fx + %.2f', coeffs(2), coeffs(1));
legend({'Data', eqnStr, '95% conf. bounds'}, 'Location', 'best')

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.mld(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.mld(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('mld (m)')
ylim([0 200])
%xlim([0 200])
title(append('Winter, totals-', loc))
coeffs = mdl.Coefficients.Estimate;
eqnStr = sprintf('y = %.2fx + %.2f', coeffs(2), coeffs(1));
legend({'Data', eqnStr, '95% conf. bounds'}, 'Location', 'best')

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Win' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.mld(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.mld(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('mld (m)')
ylim([0 200])
%xlim([0 200])
title(append('Winter, totals-', loc))
coeffs = mdl.Coefficients.Estimate;
eqnStr = sprintf('y = %.2fx + %.2f', coeffs(2), coeffs(1));
legend({'Data', eqnStr, '95% conf. bounds'}, 'Location', 'best')

%% spring
frac = 'C_chl_ratio_total';
tiledlayout(1, 3)
nexttile;
loc = 'innershelf';
req = (seasons.season == 'Spr' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.mld(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.mld(req), 'linear');
plot(mdl)
coeffs = mdl.Coefficients.Estimate;
eqnStr = sprintf('y = %.2fx + %.2f', coeffs(2), coeffs(1));
xlabel('C:chl a')
ylabel('mld (m)')
ylim([0 100])
xlim([0 200])
title(append('Spring, totals-', loc))
legend({'Data', eqnStr, '95% conf. bounds'}, 'Location', 'best')

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Spr' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.mld(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.mld(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('mld (m)')
ylim([0 100])
xlim([0 200])
title(append('Spring, totals-', loc))
coeffs = mdl.Coefficients.Estimate;
eqnStr = sprintf('y = %.2fx + %.2f', coeffs(2), coeffs(1));
legend({'Data', eqnStr, '95% conf. bounds'}, 'Location', 'best')

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Spr' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.mld(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.mld(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('mld (m)')
ylim([0 100])
xlim([0 200])
title(append('Spring, totals-', loc))
coeffs = mdl.Coefficients.Estimate;
eqnStr = sprintf('y = %.2fx + %.2f', coeffs(2), coeffs(1));
legend({'Data', eqnStr, '95% conf. bounds'}, 'Location', 'best')

%% summer
frac = 'C_chl_ratio_total';
tiledlayout(1, 3)
nexttile;
loc = 'innershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.mld(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.mld(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('mld (m)')
ylim([0 50])
xlim([0 200])
title(append('Summer, totals-', loc))
coeffs = mdl.Coefficients.Estimate;
eqnStr = sprintf('y = %.2fx + %.2f', coeffs(2), coeffs(1));
legend({'Data', eqnStr, '95% conf. bounds'}, 'Location', 'best')

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.mld(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.mld(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('mld (m)')
ylim([0 50])
xlim([0 200])
title(append('Summer, totals-', loc))
coeffs = mdl.Coefficients.Estimate;
eqnStr = sprintf('y = %.2fx + %.2f', coeffs(2), coeffs(1));
legend({'Data', eqnStr, '95% conf. bounds'}, 'Location', 'best')

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Sum' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.mld(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.mld(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('mld (m)')
ylim([0 50])
xlim([0 200])
title(append('Summer, totals-', loc))
coeffs = mdl.Coefficients.Estimate;
eqnStr = sprintf('y = %.2fx + %.2f', coeffs(2), coeffs(1));
legend({'Data', eqnStr, '95% conf. bounds'}, 'Location', 'best')

%% fall
frac = 'C_chl_ratio_total';
tiledlayout(1, 3)
nexttile;
loc = 'innershelf';
req = (seasons.season == 'Fal' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.mld(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.mld(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('mld (m)')
ylim([0 100])
xlim([0 200])
title(append('Fall, totals-', loc))
coeffs = mdl.Coefficients.Estimate;
eqnStr = sprintf('y = %.2fx + %.2f', coeffs(2), coeffs(1));
legend({'Data', eqnStr, '95% conf. bounds'}, 'Location', 'best')

nexttile;
loc = 'midshelf';
req = (seasons.season == 'Fal' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.mld(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.mld(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('mld (m)')
ylim([0 100])
xlim([0 200])
title(append('Fall, totals-', loc))
coeffs = mdl.Coefficients.Estimate;
eqnStr = sprintf('y = %.2fx + %.2f', coeffs(2), coeffs(1));
legend({'Data', eqnStr, '95% conf. bounds'}, 'Location', 'best')

nexttile;
loc = 'outershelf';
req = (seasons.season == 'Fal' & ~isoutlier(cchlQC.(frac)) & distfromshore == loc);
scatter(cchlQC.(frac)(req), cchlQC.mld(req), 'filled')
set(gca, 'FontSize', 14)
mdl = fitlm(cchlQC.(frac)(req), cchlQC.mld(req), 'linear');
plot(mdl)
xlabel('C:chl a')
ylabel('mld (m)')
ylim([0 100])
xlim([0 200])
title(append('Fall, totals-', loc))
coeffs = mdl.Coefficients.Estimate;
eqnStr = sprintf('y = %.2fx + %.2f', coeffs(2), coeffs(1));
legend({'Data', eqnStr, '95% conf. bounds'}, 'Location', 'best')

%% or plot a single shelf region by season


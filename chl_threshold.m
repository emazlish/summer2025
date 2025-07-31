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

allCHL = table();
for i = 1:length(allcruises)
   cruise = string(allcruises(i));
   CHL2 = readtable(append('https://nes-lter-data.whoi.edu/api/chl/', cruise, '.csv'));
   allCHL = [allCHL; CHL2];
   clear CHL2
end

idx = ismember(CHL(:, {'cruise', 'cast', 'niskin'}), cchltransect(:, {'cruise','cast','niskin'}));
CHL = CHL(idx == 1, :);
%%
tf = strcmp(allCHL.filter_size, '>20');
sizefrac = allCHL(tf == 1, :);
%sizefrac.rbchl = sizefrac.rb ./ sizefrac.chl;

sizefrac = sortrows(sizefrac, {'date', 'cruise', 'cast', 'niskin'});

% Step 2: Assign sample numbers based on cruise/cast/niskin grouping
[groupID, ~] = findgroups(sizefrac(:, {'date', 'cruise', 'cast', 'niskin'}));
sizefrac.sample_num = groupID;

%repA = sizefrac(strcmp(sizefrac.replicate, 'a') == 1, :);
%repB = sizefrac(strcmp(sizefrac.replicate, 'b') == 1, :);
sizefrac.replicate = categorical(sizefrac.replicate);
scatter(sizefrac.sample_num, sizefrac.chl, 40, sizefrac.replicate, 'filled')
%% qc chl data: if chl is within 2x blank, it gets removed, also bad flags are removed

exclude = allCHL(allCHL.chl <= (1.5 * allCHL.blank), :);
CHLqc = allCHL(allCHL.chl > (1.5 * allCHL.blank), :);

%%
scatter(sizefrac.chl, sizefrac.phaeo, 40, sizefrac.rbchl, 'filled')
xlabel('chl a (µg/L)')
ylabel('phaeo')
title('<10, EN695 and later')
c = colorbar;
colormap(flipud(parula))
c.Direction = 'reverse';
c.Label.String = 'rb:chl';
set(gca, 'FontSize', 14)
%mdl = fitlm(sizefrac.chl, sizefrac.rb, "linear");
%plot(mdl)

%sizefrac.rb_rbblank = sizefrac.rb ./ sizefrac.rb_blank;

%% rb subplots

subplot(2, 3, 1)
cruise = 'EN695';
tf = strcmp(string(sizefrac.cruise), cruise);
scatter(sizefrac.chl(tf == 1), sizefrac.rb(tf == 1), 40, sizefrac.ratio(tf == 1), 'filled')
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
tf = strcmp(string(sizefrac.cruise), cruise);
scatter(sizefrac.chl(tf == 1), sizefrac.rb(tf == 1), 40, sizefrac.ratio(tf == 1), 'filled')
xlabel('chl a (µg/L)')
ylabel('rb')
ylim([0 500])
title(append('totals, ', cruise))
c = colorbar;
c.Label.String = 'rb:chl a';
clear tf cruise

subplot(2, 3, 3)
cruise = 'EN712';
tf = strcmp(string(sizefrac.cruise), cruise);
scatter(sizefrac.chl(tf == 1), sizefrac.rb(tf == 1), 40, sizefrac.ratio(tf == 1), 'filled')
xlabel('chl a (µg/L)')
ylabel('rb')
ylim([0 500])
title(append('totals, ', cruise))
c = colorbar;
c.Label.String = 'rb:chl a';
clear tf cruise

subplot(2, 3, 4)
cruise = 'EN715';
tf = strcmp(string(sizefrac.cruise), cruise);
scatter(sizefrac.chl(tf == 1), sizefrac.rb(tf == 1), 40, sizefrac.ratio(tf == 1), 'filled')
xlabel('chl a (µg/L)')
ylabel('rb')
ylim([0 500])
title(append('totals, ', cruise))
c = colorbar;
c.Label.String = 'rb:chl a';
clear tf cruise

subplot(2, 3, 5)
cruise = 'EN720';
tf = strcmp(string(sizefrac.cruise), cruise);
scatter(sizefrac.chl(tf == 1), sizefrac.rb(tf == 1), 40, sizefrac.ratio(tf == 1), 'filled')
xlabel('chl a (µg/L)')
ylabel('rb')
ylim([0 500])
title(append('totals, ', cruise))
c = colorbar;
c.Label.String = 'rb:chl a';
clear tf cruise

% exclusion criterion: any value for totals.ratio < 100
excluded = sizefrac(sizefrac.ratio < 100, :);
included = sizefrac(sizefrac.ratio > 100, :);

% histogram of excluded values
histogram(excluded.rb)
xlabel('Chl a (µg/L)')
ylabel('count')
title('excluded values based on rb:chl a < 100')

%% hist with qc'd chl values grayed out

% totals
subplot(3, 5, 11)
loc = 'outershelf';
input = totals.C_chl_ratio_total(totals.dist == loc & totals.chl_0_avg > 0.1);
histogram(input, -200:5:200, 'FaceColor', 'blue')
hold on
input2 = totals.C_chl_ratio_total(totals.dist == loc & totals.chl_0_avg <= 0.1);
histogram(input2, -200:5:200, 'FaceColor', 'red')
xlabel('C:chl a')
ylabel('Count')
title(append('totals-', loc))

% lessthan5
subplot(3, 5, 12)
loc = 'outershelf';
input = lessthan5.C_chl_ratio_lessthan5(lessthan5.dist == loc & lessthan5.chl_lessthan5 > 0.1);
histogram(input, -200:5:200, 'FaceColor', 'blue')
hold on
input2 = lessthan5.C_chl_ratio_lessthan5(lessthan5.dist == loc & lessthan5.chl_lessthan5 <= 0.1);
histogram(input2, -200:5:200, 'FaceColor', 'red')
xlabel('C:chl a')
ylabel('Count')
title(append('lessthan5-', loc))

% 5to10
subplot(3, 5, 13)
loc = 'outershelf';
input = fivetoten.C_chl_ratio_5to10(fivetoten.dist == loc & fivetoten.chl_5to10 > 0.1);
histogram(input, -200:5:200, 'FaceColor', 'blue')
hold on
input2 = fivetoten.C_chl_ratio_5to10(fivetoten.dist == loc & fivetoten.chl_5to10 <= 0.1);
histogram(input2, -200:5:200, 'FaceColor', 'red')
xlabel('C:chl a')
ylabel('Count')
title(append('5to10-', loc))

% 10 to 20
subplot(3, 5, 14)
loc = 'outershelf';
input = tentotwenty.C_chl_ratio_10to20(tentotwenty.dist == loc & tentotwenty.chl_10to20 > 0.1);
histogram(input, -200:5:200, 'FaceColor', 'blue')
hold on
input2 = tentotwenty.C_chl_ratio_10to20(tentotwenty.dist == loc & tentotwenty.chl_10to20 <= 0.1);
histogram(input2, -200:5:200, 'FaceColor', 'red')
xlabel('C:chl a')
ylabel('Count')
title(append('10to20-', loc))

% greater than 20
subplot(3, 5, 15)
loc = 'outershelf';
input = greaterthan20.C_chl_ratio_greaterthan20(greaterthan20.dist == loc & greaterthan20.chl_20_avg > 0.1);
histogram(input, -200:5:200, 'FaceColor', 'blue')
hold on
input2 = greaterthan20.C_chl_ratio_greaterthan20(greaterthan20.dist == loc & greaterthan20.chl_20_avg <= 0.1);
histogram(input2, -200:5:200, 'FaceColor', 'red')
xlabel('C:chl a')
ylabel('Count')
title(append('10to20-', loc))

%% scatter of c vs. chl by shelf and size

% totals
subplot(3, 5, 11)
loc = 'outershelf';
req = totals.dist == loc & totals.chl_0_avg > 0.1;
scatter(totals.sumC_total(req), totals.chl_0_avg(req), 40, 'filled')
hold on
req2 = totals.dist == loc & totals.chl_0_avg <= 0.1;
scatter(totals.sumC_total(req2), totals.chl_0_avg(req2), 40, 'filled')

% lessthan5
subplot(3, 5, 12)
loc = 'outershelf';
req = lessthan5.dist == loc & lessthan5.chl_lessthan5 > 0.1;
scatter(lessthan5.sumC_5(req), lessthan5.chl_lessthan5(req), 40, 'filled')
hold on
req2 = lessthan5.dist == loc & lessthan5.chl_lessthan5 <= 0.1;
scatter(lessthan5.sumC_5(req2), lessthan5.chl_lessthan5(req2), 40, 'filled')
title(append('lessthan5-', loc))
xlabel('chl a')
ylabel('carbon')

% 5 to 10
subplot(3, 5, 13)
loc = 'outershelf';
req = fivetoten.dist == loc & fivetoten.chl_5to10 > 0.1;
scatter(fivetoten.sumC_5to10(req), fivetoten.chl_5to10(req), 40, 'filled')
hold on
req2 = fivetoten.dist == loc & fivetoten.chl_5to10 <= 0.1;
scatter(fivetoten.sumC_5to10(req2), fivetoten.chl_5to10(req2), 40, 'filled')
title(append('lessthan5-', loc))
xlabel('chl a')
ylabel('carbon')

% 10 to 20
subplot(3, 5, 14)
loc = 'outershelf';
req = tentotwenty.dist == loc & tentotwenty.chl_10to20 > 0.1;
scatter(tentotwenty.sumC_10to20(req), tentotwenty.chl_10to20(req), 40, 'filled')
hold on
req2 = tentotwenty.dist == loc & tentotwenty.chl_10to20 <= 0.1;
scatter(tentotwenty.sumC_10to20(req2), tentotwenty.chl_10to20(req2), 40, 'filled')
title(append('10to20-', loc))
xlabel('chl a')
ylabel('carbon')

% greaterthan20
subplot(3, 5, 15)
loc = 'outershelf';
req = greaterthan20.dist == loc & greaterthan20.chl_20_avg > 0.1;
scatter(greaterthan20.sumC_greaterthan20(req), greaterthan20.chl_20_avg(req), 40, 'filled')
hold on
req2 = greaterthan20.dist == loc & greaterthan20.chl_20_avg <= 0.1;
scatter(greaterthan20.sumC_greaterthan20(req2), greaterthan20.chl_20_avg(req2), 40, 'filled')
title(append('greaterthan20-', loc))
xlabel('chl a')
ylabel('carbon')
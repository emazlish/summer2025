% generate a table for C:chl that excludes chl values <= 1.5*blank

cd /Users/emazlish/CodeGithub/summer2025/

% load carbon data from Attune/IFCB merge
load('/Volumes/Lab_data/Attune/cruise_data/IFCB_Attune_merge/summary_files_discrete/Attune_IFCB_merge_class_summary_discrete_carbon2_phyto_allcruises.mat')

% load label options
load /Users/emazlish/Library/CloudStorage/OneDrive-BowdoinCollege/WHOI_2025/Code/label_options.mat

% load all chl files for cruises from the API (this .mat includes cruises
% EN608-EN720 in chronological order)
load /Users/emazlish/Library/CloudStorage/OneDrive-BowdoinCollege/WHOI_2025/Datasets/new_Cchltable_outputs/allchl_API.mat

% need to add 'label options' to this file
%% arrange CHL data into one big file

cruises = ismember(CHL.cruise, phytoC_all.cruise);
CHL = CHL(cruises, :);
clear cruises

% create string with cruises after EN687 with C and CHL data, or whichever
% hasn't been loaded yet
newercruises = {'EN695', 'EN706', 'EN712', 'EN715', 'EN720', 'EN727'};
% can't read in EN727 yet because chl data doesn't exist for it yet

for i = 1:length(newercruises)
    cruise = string(newercruises(i));
    CHL2 = readtable(append('https://nes-lter-data.whoi.edu/api/chl/', cruise, '.csv'));
    CHL2.iode_quality_flag = CHL2.quality_flag;
    CHL2.date_time_utc = CHL2.date;
    CHL2(:, ["vol_filtered", "tau_calibration", "fd_calibration", "rb", "ra", "blank", "rb_blank", "ra_blank", "quality_flag", "date"]) = [];
    varsToAdd = CHL.Properties.VariableNames(~ismember(CHL.Properties.VariableNames, CHL2.Properties.VariableNames));
    t = table(nan(size(CHL2, 1), 1));
    t = addvars(t, nan(size(CHL2, 1), 1));
    t = addvars(t, nan(size(CHL2, 1), 1));
    t = addvars(t, nan(size(CHL2, 1), 1));
    t.Properties.VariableNames = varsToAdd;
    CHL2 = [CHL2, t];
    CHL2 = CHL2(:, string(CHL.Properties.VariableNames));
    CHL2.date_time_utc = datetime(CHL2.date_time_utc, 'InputFormat', 'yyyy-MM-dd HH:mm:ss+00:00');
    CHL2.method_contributor = num2cell(CHL2.method_contributor);
    CHL2.project_id = num2cell(CHL2.project_id);
    CHL2.nearest_station = num2cell(CHL2.nearest_station);
    CHL2.distance = num2cell(CHL2.distance);
    CHL = [CHL; CHL2];
    clear CHL2 t
end

clear newercruises varsToAdd
%% create label for each permutation of chl filter size and replicate

label_options = strings(15, 2);

label_options(1:15, 1) = ["chl>0&<10_a"; "chl>0&<10_b"; "chl>0&<10_c"; "chl>0&<200_a"; "chl>0&<200_b"; "chl>0&<200_c"; "chl>10&<200_a"; "chl>10&<200_b"; "chl>10&<200_c"; "chl>20&<200_a"; "chl>20&<200_b"; "chl>20&<200_c"; "chl>5&<200_a"; "chl>5&<200_b"; "chl>5&<200_c"];
label_options(1:15, 2) =  ["chl<10_a"; "chl<10_b"; "chl<10_c"; "chl>0_a"; "chl>0_b"; "chl>0_c"; "chl_greaterthan10_a"; "chl_greaterthan10_b"; "chl_greaterthan10_c"; "chl>20_a"; "chl>20_b"; "chl>20_c"; "chl>5_a"; "chl>5_b";"chl>5_c"];
label_options = cellstr(label_options);

allCHL.label = append('chl', allCHL.filter_size, '_', allCHL.replicate);
allCHL.label2 = (allCHL.label);

for i = 1:size(allCHL, 1)
    allCHL.label2(i) = allCHL.label(i);
    q = strcmp(label_options(:, 1), allCHL.label2(i));
    if any(q)
        [~, idx] = ismember(allCHL.label(i), label_options(:, 1));
         c = label_options(idx, 2);
        allCHL.label2(i) = c;
     end
    clear idx q c
end

%% QC chl data and unstack CHL by labels

% remove chl samples that have bad quality flags (flag = 3 or 4)
flagidx = (allCHL.quality_flag == 3 | allCHL.quality_flag == 4);
allCHL = allCHL(flagidx == 0, :);
clear flagidx

% remove samples that have rb <= 3 * blank

idx = allCHL.rb <= 3 * allCHL.blank;
CHLqc = allCHL(idx == 0, :);
clear idx

% separate samples based on whether chl value is above/below 0.1 µg/L
chlToKeep = CHLqc(CHLqc.chl > 0.1, :); % above threshold, get used in chlToKeep
belowThresh = CHLqc(CHLqc.chl <= 0.1, :); % below threshold, don't get used

%% unstack by label
q = chlToKeep(:, ["cruise", "cast", "niskin", "chl", "label2"]); % pull cruise, cast, niskin, chl and label from CHL
e = unstack(q, "chl", "label2"); % unstack chl data by sorting by label

clear q
%%
phytoC_all.cruise = string(phytoC_all.cruise);
e.cruise = string(e.cruise);

% figure out which index in chl table (e) each phytoC row belongs to
[~, ia] = ismember(e(:, {'cruise', 'cast', 'niskin'}), ...
                   phytoC_all(:, {'cruise', 'cast', 'niskin'}));

valid = ia > 0; % removes 0 values which throw off the index

% sort phytoC by ordering of chl table (i.e. chronologically)
sorted_phytoC = phytoC_all(ia(valid), :);

clear ia valid

%% join carbon and chl data tables into C_chl master table

cchlQC = innerjoin(sorted_phytoC, e, 'Keys', {'cruise', 'cast', 'niskin'});

cchlQC = sortrows(cchlQC, {'date_sampled', 'cast', 'niskin'});

cchlQC.dt = [];
cchlQC.mdate = [];

%% only keep samples that have both Attune and IFCB information
noIFCB = isnan(cchlQC.C_100_Inf);
cchlQC = cchlQC(noIFCB == 0, :);

%% Chl average calculations between two replicates

% average chl between two replicates for >0 µm chl
% including the "2" in the second position finds mean by row
cchlQC.chl_0_avg = mean([cchlQC.chl_0_a cchlQC.chl_0_b], 2, 'omitnan');

% average chl between two replicates for >5 µm chl
% including the "2" in the second position finds mean by row

cchlQC.chl_5_avg = mean([cchlQC.chl_5_a cchlQC.chl_5_b], 2, 'omitnan');

% average chl between two replicates for <10 µm chl
% including the "2" in the second position finds mean by row
cchlQC.chl_10_avg = mean([cchlQC.chl_10_a cchlQC.chl_10_b], 2, 'omitnan');

% average chl between two replicates for >20 µm chl
% including the "2" in the second position finds mean by row
cchlQC.chl_20_avg = mean([cchlQC.chl_20_a cchlQC.chl_20_b], 2, 'omitnan');

%% Chl calculations for size fractions that are not directly given by chl

% <5 µm chl fraction
cchlQC.chl_lessthan5 = cchlQC.chl_0_avg - cchlQC.chl_5_avg;

% 5-10 µm chl fraction
cchlQC.chl_5to10 = cchlQC.chl_10_avg - cchlQC.chl_lessthan5;

% <20 µm chl fraction
cchlQC.chl_lessthan20 = cchlQC.chl_0_avg - cchlQC.chl_20_avg;

% 10-20 µm chl fraction
cchlQC.chl_10to20 = cchlQC.chl_lessthan20 - cchlQC.chl_10_avg;

% 5 to 20 µm chl fraction
cchlQC.chl_5to20 = cchlQC.chl_5_avg - cchlQC.chl_20_avg;
%% sum C for each size bin and add as new column

% total C
cchlQC.sumC_total = real(table2array(sum(cchlQC(:, 11:22), 2, 'omitnan')));

% < 5 µm fraction
cchlQC.sumC_5 = real(table2array(sum(cchlQC(:,11:13), 2, 'omitnan')));

% < 10 µm fraction
cchlQC.sumC_10 = real(table2array(sum(cchlQC(:,11:15), 2, 'omitnan')));

% < 20 µm fraction
cchlQC.sumC_lessthan20 = real(table2array(sum(cchlQC(:,11:17), 2, 'omitnan')));

% > 20 µm fraction
cchlQC.sumC_greaterthan20 = real(table2array(sum(cchlQC(:,18:22), 2, 'omitnan')));

% 5-10 µm fraction
cchlQC.sumC_5to10 = real(table2array(sum(cchlQC(:,14:15), 2, 'omitnan')));

% 10-20 µm fraction
cchlQC.sumC_10to20 = real(table2array(sum(cchlQC(:,16:17), 2, 'omitnan')));

% 5-20 µm fraction
cchlQC.sumC_5to20 = real(table2array(sum(cchlQC(:, 14:17), 2, 'omitnan')));

% greater than 5 fraction
cchlQC.sumC_greaterthan5 = real(table2array(sum(cchlQC(:, 14:22), 2, 'omitnan')));
%% C:chl ratios

% C:chl for >0 µm fraction (total)
cchlQC.C_chl_ratio_total = cchlQC.sumC_total ./ cchlQC.chl_0_avg;

% C:chl for <10 µm fraction
cchlQC.C_chl_ratio_10 = cchlQC.sumC_10 ./ cchlQC.chl_10_avg;

% C:chl for <5 µm fraction
cchlQC.C_chl_ratio_lessthan5 = cchlQC.sumC_5 ./ cchlQC.chl_lessthan5;

% C:chl for <20 µm fraction
cchlQC.C_chl_ratio_lessthan20 = cchlQC.sumC_lessthan20 ./ cchlQC.chl_lessthan20;

% C:chl for >20 µm fraction
cchlQC.C_chl_ratio_greaterthan20 = cchlQC.sumC_greaterthan20 ./ cchlQC.chl_20_avg;

% C:chl for 10-20 fraction
cchlQC.C_chl_ratio_10to20 = cchlQC.sumC_10to20 ./ cchlQC.chl_10to20;

% C:chl for 5-10 µm fraction
cchlQC.C_chl_ratio_5to10 = cchlQC.sumC_5to10 ./ cchlQC.chl_5to10;

% C:chl for 5–20µm fraction
cchlQC.C_chl_ratio_5to20 = cchlQC.sumC_5to20 ./ cchlQC.chl_5to20;

% C:chl for >5 µm fraction
cchlQC.C_chl_ratio_greaterthan5 = cchlQC.sumC_greaterthan5 ./ cchlQC.chl_5_avg;

%% remove depths deeper than 100m

cchlQC = cchlQC(cchlQC.depth_m <= 100, :);

%% save files

save QC_data_0731.mat cchlQC allCHL chlToKeep belowThresh phytoC_all

%% assign months and distfromshore

m = datetime(cchlQC.date_sampled, 'InputFormat','yyyy-MM-dd HH:mm:ss+00:00');

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
seasons.monthcat = categorical(seasons.month);

distfromshore = strings(size(cchlQC.nearest_station, 1), 1);
for i = 1:height(distfromshore)
    if ismember(cchlQC.nearest_station(i), {'MVCO','L1', 'L2'}) == 1
        distfromshore(i) = 'innershelf';
    end
    if ismember(cchlQC.nearest_station(i), {'L3', 'L4','L5','L6', 'L7', 'L8', 'L9'}) == 1
        distfromshore(i) = 'midshelf';
    end
    if ismember(cchlQC.nearest_station(i), {'L10', 'L11'}) == 1
        distfromshore(i) = 'outershelf';
    end
end

%% 7/30 BOX and SCATTER WITH ALL MONTHS REPRESENTED (ADD EMPTY MONTHS AS NANS, CONC. ARRAYS)
% overview plot, ≤100 m
monthNames = {'Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'};

frac = 'C_chl_ratio_total';
concreq = [cchlQC.(frac); nan(12, 1)];
concmon = [seasons.monthcat; monthNames];
concmon = categorical(concmon);
concdpth = [cchlQC.depth_m; nan(12, 1)];

ax1 = axes();
scatter(ax1, concmon, concreq, 40, concdpth, "filled", 'o', 'jitter', 'on', 'jitteramount', 0.1, 'Clipping', 'on')
ax1.XAxis.Categories = {'Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'};
set(ax1, 'YLim', [-20 350])
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
boxplot(ax2, concreq, concmon, 'GroupOrder', monthNames, 'OutlierSize', 10,'DataLim', [0 350],'ExtremeMode', 'clip')
ax2.YLim = get(ax1, 'YLim');
dashedLines = findobj(ax2, 'Type', 'Line', 'Tag', '','LineStyle', '--');
delete(dashedLines)
%title(ax1, append('C:chl a by month, ', frac,' ≤100 m depth'), 'FontWeight', 'bold')
set(ax1, 'FontSize', 16)
xlabel(ax1, 'Month')
ylabel(ax1, 'C:chl a')

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
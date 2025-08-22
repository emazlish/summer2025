% This script allows you to pull chlorophyll data from the NESLTER
% REST API and carbon merge data from the IFCB and Attune to generate
% a table with carbon-to-chlorophyll ratios of different size 
% fractions.
% Last edited by Emma Mazlish on 20250822

% load the data - change filepaths as necessary

% load carbon data from Attune/IFCB merge
load('/Volumes/Lab_data/Attune/cruise_data/IFCB_Attune_merge/summary_files_discrete/Attune_IFCB_merge_class_summary_discrete_carbon2_phyto_allcruises.mat')

%% load the chlorophyll data for each cruise of interest from API and

% can pull cruises that match phytoC data
cruises = unique(phytoC_all.cruise);

% or, pull data for specific cruises
cruises = {'EN695', 'EN706', 'EN712', 'EN715', 'EN720', 'EN727'}; % for example

CHL = table(); % initialize CHL master table

% access each cruise file on the API and concatenate data to existing table
for i = 1:length(cruises)
    cruise = string(cruises(i));
    CHL2 = readtable(append('https://nes-lter-data.whoi.edu/api/chl/', cruise, '.csv'));
    
    % choose which variables to keep
    CHL2 = CHL2(:, {'cruise','cast', 'niskin', 'replicate', 'vol_filtered', ...
        'filter_size', 'rb', 'ra', 'blank', 'chl', 'phaeo','quality_flag', ...
        'date','longitude', 'latitude', 'depth'});

    CHL = [CHL; CHL2]; % concatenate with existing table
    clear CHL2 t
end

clear cruises varsToAdd

%% create label dictionary for each chl replicate/filter size option

label_options = strings(15, 2);

label_options(1:15, 1) = ["chl>0&<10_a"; "chl>0&<10_b"; "chl>0&<10_c"; "chl>0&<200_a"; "chl>0&<200_b"; "chl>0&<200_c"; "chl>10&<200_a"; "chl>10&<200_b"; "chl>10&<200_c"; "chl>20&<200_a"; "chl>20&<200_b"; "chl>20&<200_c"; "chl>5&<200_a"; "chl>5&<200_b"; "chl>5&<200_c"];
label_options(1:15, 2) =  ["chl<10_a"; "chl<10_b"; "chl<10_c"; "chl>0_a"; "chl>0_b"; "chl>0_c"; "chl_greaterthan10_a"; "chl_greaterthan10_b"; "chl_greaterthan10_c"; "chl>20_a"; "chl>20_b"; "chl>20_c"; "chl>5_a"; "chl>5_b";"chl>5_c"];
label_options = cellstr(label_options);

CHL.label = append('chl', CHL.filter_size, '_', CHL.replicate); % create initial label for each row
CHL.label2 = (CHL.label); % duplicate that label

% create a second label (CHL.label2) that has a MATLAB-readable version of
% CHL.label for use in the unstack function later
for i = 1:size(CHL, 1)
    CHL.label2(i) = CHL.label(i);
    q = strcmp(label_options(:, 1), CHL.label2(i));
    if any(q)
        [~, idx] = ismember(CHL.label(i), label_options(:, 1));
         c = label_options(idx, 2);
        CHL.label2(i) = c;
     end
    clear idx q c
end

%% QC chl data as desired

% remove chl samples that have bad quality flags (flag = 3 or 4)
flagidx = (CHL.quality_flag == 3 | CHL.quality_flag == 4);
CHL = CHL(flagidx == 0, :);
clear flagidx

% remove samples that have rb <= 3 * blank
idx = CHL.rb <= 3 * CHL.blank;
CHLqc = CHL(idx == 0, :);
clear idx

% separate samples based on whether chl value is above/below 0.1 µg/L
chlToKeep = CHLqc(CHLqc.chl > 0.1, :); % above threshold, get used in chlToKeep
belowThresh = CHLqc(CHLqc.chl <= 0.1, :); % below threshold, don't get used

%% unstack chl by label
q = chlToKeep(:, ["cruise", "cast", "niskin", "chl", "label2"]); % pull cruise, cast, niskin, chl and MATLAB-readable label from CHL
e = unstack(q, "chl", "label2"); % unstack chl data by sorting by label
clear q

%% sort and arrange phytoC table to align with CHL table
phytoC_all.cruise = string(phytoC_all.cruise); % standardize data format
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

cchlQC = sortrows(cchlQC, {'date_sampled', 'cast', 'niskin'}); % order chronologically

clear e

%% only keep samples that have both Attune and IFCB information
noIFCB = isnan(cchlQC.C_100_Inf);
cchlQC = cchlQC(noIFCB == 0, :);

%% Chl average calculations between two replicates

% average chl between two replicates for >0 µm chl
% including the "2" in the second position finds mean by row
cchlQC.chl_0_avg = mean([cchlQC.chl_0_a cchlQC.chl_0_b], 2, 'omitnan');

% average chl between two replicates for >5 µm chl
cchlQC.chl_5_avg = mean([cchlQC.chl_5_a cchlQC.chl_5_b], 2, 'omitnan');

% average chl between two replicates for <10 µm chl
cchlQC.chl_10_avg = mean([cchlQC.chl_10_a cchlQC.chl_10_b], 2, 'omitnan');

% average chl between two replicates for >20 µm chl
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

% 5-20 µm chl fraction
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

% > 5 fraction
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

%% remove depths deeper than 100m if desired

cchlQC = cchlQC(cchlQC.depth_m <= 100, :);

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

% add to cchlQC table
cchlQC.month = seasons.month;
cchlQC.monthcat = seasons.monthcat; % can be useful to have months as categorical

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

cchlQC.distfromshore = distfromshore; % add to cchlQC table

%% save files

save [filename].mat cchlQC CHL chlToKeep belowThresh phytoC_all
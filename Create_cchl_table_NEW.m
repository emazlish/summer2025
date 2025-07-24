% new (20250723) C_chl table creation using updated IFCB/Attune merge dataset 
% this script looks at carbon across the IFCB and Attune for discrete
% samples as well as chlorophyll for an entire cruise to generate C:chl 
% ratios for the NES

% file is saved in github access folder
cd /Users/emazlish/CodeGithub/summer2025/

%cruise = "EN715"; % this is the only thing you have to change to run script

% load carbon data from Attune/IFCB merge
load('/Volumes/Lab_data/Attune/cruise_data/IFCB_Attune_merge/summary_files_discrete/Attune_IFCB_merge_class_summary_discrete_carbon2_phyto_allcruises.mat')

% load label options
load /Users/emazlish/Library/CloudStorage/OneDrive-BowdoinCollege/WHOI_2025/Code/label_options.mat

% load CHL data from API for cruises more recent than EN687
%CHL2 = readtable(append('https://nes-lter-data.whoi.edu/api/chl/', cruise, '.csv'));

% load CHL data from EDI for cruises through EN687
CHL = readtable('/Users/emazlish/Library/CloudStorage/OneDrive-BowdoinCollege/WHOI_2025/Datasets/20250604_practice_POC_CHL/nes-lter-chl-transect-rosette.csv');

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

CHL.label = append('chl', CHL.filter_size, '_', CHL.replicate);
CHL.label2 = (CHL.label);

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

%% QC chl data and unstack CHL by labels

% remove chl samples that have bad quality flags (flag = 3 or 4)
flagidx = (CHL.iode_quality_flag == 3 | CHL.iode_quality_flag == 4);
CHL = CHL(flagidx == 0, :);
clear flagidx

% unstack by label
q = CHL(:, ["cruise", "cast", "niskin", "chl", "label2"]); % pull cruise, cast, niskin, chl and label from CHL
e = unstack(q, "chl", "label2"); % unstack chl data by sorting by label

clear q
%% sort carbon table in chronological order by date to match CHL table

% turn both cruise columns into strings (easier for ismember)
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

C_chl = innerjoin(sorted_phytoC, e, 'Keys', {'cruise', 'cast', 'niskin'});

C_chl = sortrows(C_chl, {'date_sampled', 'cast', 'niskin'});

% clear extraneous columns
C_chl.dt = []; 
C_chl.mdate = [];
C_chl.chl_10_c = [];
C_chl.chl_0_c = [];
C_chl.chl_5_c = [];
C_chl.chl_20_c = [];
C_chl.chl_greaterthan10_c = [];
clear e sorted_phytoC

%% sum C for each size bin and add as new column

% total C
C_chl.sumC_total = real(table2array(sum(C_chl(:, 11:22), 2, 'omitnan')));

% < 5 µm fraction
C_chl.sumC_5 = real(table2array(sum(C_chl(:,11:13), 2, 'omitnan')));

% < 10 µm fraction
C_chl.sumC_10 = real(table2array(sum(C_chl(:,11:15), 2, 'omitnan')));

% < 20 µm fraction
C_chl.sumC_lessthan20 = real(table2array(sum(C_chl(:,11:17), 2, 'omitnan')));

% > 20 µm fraction
C_chl.sumC_greaterthan20 = real(table2array(sum(C_chl(:,18:22), 2, 'omitnan')));

% 5-10 µm fraction
C_chl.sumC_5to10 = real(table2array(sum(C_chl(:,14:15), 2, 'omitnan')));

% 10-20 µm fraction
C_chl.sumC_10to20 = real(table2array(sum(C_chl(:,16:17), 2, 'omitnan')));

%% Chl average calculations between two replicates

% average chl between two replicates for >0 µm chl
% including the "2" in the second position finds mean by row
C_chl.chl_0_avg = mean([C_chl.chl_0_a C_chl.chl_0_b], 2, 'omitnan');

% average chl between two replicates for >5 µm chl
% including the "2" in the second position finds mean by row

C_chl.chl_5_avg = mean([C_chl.chl_5_a C_chl.chl_5_b], 2, 'omitnan');

% average chl between two replicates for <10 µm chl
% including the "2" in the second position finds mean by row
C_chl.chl_10_avg = mean([C_chl.chl_10_a C_chl.chl_10_b], 2, 'omitnan');

% average chl between two replicates for >20 µm chl
% including the "2" in the second position finds mean by row
C_chl.chl_20_avg = mean([C_chl.chl_20_a C_chl.chl_20_b], 2, 'omitnan');

%% Chl calculations for size fractions that are not directly given by chl

% <5 µm chl fraction
C_chl.chl_lessthan5 = C_chl.chl_0_avg - C_chl.chl_5_avg;

% 5-10 µm chl fraction
C_chl.chl_5to10 = C_chl.chl_10_avg - C_chl.chl_lessthan5;

% <20 µm chl fraction
C_chl.chl_lessthan20 = C_chl.chl_0_avg - C_chl.chl_20_avg;

% 10-20 µm chl fraction
C_chl.chl_10to20 = C_chl.chl_lessthan20 - C_chl.chl_10_avg;

%% C:chl ratios

% C:chl for >0 µm fraction (total)
C_chl.C_chl_ratio_total = C_chl.sumC_total ./ C_chl.chl_0_avg;

% C:chl for <10 µm fraction
C_chl.C_chl_ratio_10 = C_chl.sumC_10 ./ C_chl.chl_10_avg;

% C:chl for <5 µm fraction
C_chl.C_chl_ratio_lessthan5 = C_chl.sumC_5 ./ C_chl.chl_lessthan5;

% C:chl for <20 µm fraction
C_chl.C_chl_ratio_lessthan20 = C_chl.sumC_lessthan20 ./ C_chl.chl_lessthan20;

% C:chl for >20 µm fraction
C_chl.C_chl_ratio_greaterthan20 = C_chl.sumC_greaterthan20 ./ C_chl.chl_20_avg;

% C:chl for 10-20 fraction
C_chl.C_chl_ratio_10to20 = C_chl.sumC_10to20 ./ C_chl.chl_10to20;

% C:chl for 5-10 µm fraction
C_chl.C_chl_ratio_5to10 = C_chl.sumC_5to10 ./ C_chl.chl_5to10;

%% save file to appropriate directory

phytocchl = C_chl; % call it whatever you want (can change if doing separate files for diatoms, dinos etc)

cd /Users/emazlish/Library/CloudStorage/OneDrive-BowdoinCollege/WHOI_2025/Datasets/new_Cchltable_outputs/

save allcruises_phytocchl.mat phytocchl % remember to change name of file
write a method section for a paper based on this script:  % load_excel_tables.m
% Script to load area and particle tables, annotate with metadata, and merge
clear all
clc
close all
%%

% Initialize master table
masterTable = table();

% Set default directory to user's home
if ispc
    currentDir = getenv('USERPROFILE');
else
    currentDir = getenv('HOME');
end

for i = 1:100
    % Select area file
    [file1, path1] = uigetfile('*.*', sprintf('Select Area File #%d', i), currentDir);
    if isequal(file1, 0)
        fprintf('No file selected. Ending process.\n');
        break;
    end
    areaFile = fullfile(path1, file1);
    currentDir = path1;  % Update path to last selected location

    % Select particles file
    [file2, path2] = uigetfile('*.*', sprintf('Select Particles File #%d', i), currentDir);
    if isequal(file2, 0)
        fprintf('No file selected. Ending process.\n');
        break;
    end
    particleFile = fullfile(path2, file2);
    currentDir = path2;  % Update path again

    % Try to read both files into tables
    try
        areaTbl = readtable(areaFile);
        particleTbl = readtable(particleFile);
    catch 
        warning('Error reading files: %s', ME.message);
        continue;
    end

    % Prompt user for metadata
    fprintf('\nProcessing files:\n  Area: %s\n  Particles: %s\n', areaFile, particleFile);
    mouse = input('Enter Mouse ID (numeric): ');
    axon  = input('Enter Axon ID (numeric): ');
    day   = input('Enter Day (numeric): ');

    % Add metadata + Type column
    areaTbl.Mouse = repmat(mouse, height(areaTbl), 1);
    areaTbl.Axon  = repmat(axon, height(areaTbl), 1);
    areaTbl.Day   = repmat(day, height(areaTbl), 1);
    areaTbl.Type  = repmat(2, height(areaTbl), 1); % 2 = Area

    particleTbl.Mouse = repmat(mouse, height(particleTbl), 1);
    particleTbl.Axon  = repmat(axon, height(particleTbl), 1);
    particleTbl.Day   = repmat(day, height(particleTbl), 1);
    particleTbl.Type  = repmat(1, height(particleTbl), 1); % 1 = Particle

    % Merge into master table
    masterTable = [masterTable; areaTbl; particleTbl];
end

% Show summary
fprintf('\nFinal table has %d rows and %d columns.\n', height(masterTable), width(masterTable));
disp(masterTable);

% Optionally save result
% writetable(masterTable, 'merged_table.xlsx');
%%
masterTable_cell{1,1} = masterTable_adolescent ;
masterTable_cell{1,2} = masterTable_adult ;

%%

clearvars -except masterTable_cell

colors_ado_adu = {[.5 .7 .2],[ .2 .4 .2]} ; 

%%
% File: analyze_synapse_data.m


group_names = {'adolescent', 'adult'};
mice_counts = [3, 4];
num_days = 4;

synapse_per_day = struct();
axon_avg_per_group = struct();
synapse_per_axon = struct();
mouse_avg_per_group = struct();

for g = 1:2
    tbl = masterTable_cell{1, g};

    if istable(tbl)
        tbl_struct = table2struct(tbl);
    else
        tbl_struct = tbl;
    end

    mouse_ids = unique([tbl_struct.Mouse]);
    axon_id_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    synapse_id_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    synapse_counter = 1;
    global_axon_id = 1;  % Unique axon ID counter across all mice

    % Assign new axon IDs across mice
    for m = 1:length(mouse_ids)
        mouse = mouse_ids(m);
        mouse_idx = find([tbl_struct.Mouse] == mouse);
        axons = unique([tbl_struct(mouse_idx).Axon]);

        for a = 1:length(axons)
            original_id = axons(a);
            key = sprintf('M%d_A%d', mouse, original_id);
            axon_id_map(key) = global_axon_id;
            global_axon_id = global_axon_id + 1;
        end
    end

    % Assign new axon and synapse IDs to each entry
    new_axon_ids = zeros(length(tbl_struct), 1);
    new_synapse_ids = zeros(length(tbl_struct), 1);
    for i = 1:length(tbl_struct)
        syn = tbl_struct(i);
        key_ax = sprintf('M%d_A%d', syn.Mouse, syn.Axon);
        key_syn = sprintf('M%d_A%d_S%d', syn.Mouse, syn.Axon, syn.Var1);

        new_axon_ids(i) = axon_id_map(key_ax);

        if ~isKey(synapse_id_map, key_syn)
            synapse_id_map(key_syn) = synapse_counter;
            synapse_counter = synapse_counter + 1;
        end
        new_synapse_ids(i) = synapse_id_map(key_syn);
    end

    % Build axon fluorescence matrix
    total_axons = global_axon_id - 1;
    axon_fluorescence = NaN(total_axons, num_days);

    % Build synapse matrix
    total_synapses = synapse_counter - 1;
    synapse_matrix = NaN(total_synapses, num_days);
    synapse_area = NaN(total_synapses, num_days);
    synapse_max = NaN(total_synapses, num_days);
    synapse_IntDen = NaN(total_synapses, num_days);
    Synapse_ID = NaN(total_synapses, num_days);  % Initialize

    for i = 1:length(tbl_struct)
        syn = tbl_struct(i);
        day = syn.Day;
        syn_idx = new_synapse_ids(i);
        axon_idx = new_axon_ids(i);

        if day >= 1 && day <= num_days
            synapse_matrix(syn_idx, day) = syn.Mean;
            Synapse_ID(syn_idx, day) = axon_idx;  % use updated axon ID
            synapse_area(syn_idx, day) = syn.Area;
            synapse_max(syn_idx, day) = syn.Max;
            synapse_IntDen(syn_idx, day) = syn.IntDen;
        end
    end

    synapse_per_day.(group_names{g}) = synapse_matrix;
    synapse_matrix_cell{1,g} = synapse_matrix;
    synapse_axon_cell{1,g} = Synapse_ID;
    synapse_area_cell{1,g} = synapse_area;
    synapse_max_cell{1,g} = synapse_max;
    synapse_IntDen_cell{1,g} = synapse_IntDen;

    % Compute axon fluorescence average
    for i = 1:length(tbl_struct)
        axon_idx = new_axon_ids(i);
        day = tbl_struct(i).Day;
        if day >= 1 && day <= num_days
            if isnan(axon_fluorescence(axon_idx, day))
                axon_fluorescence(axon_idx, day) = tbl_struct(i).Mean;
            else
                axon_fluorescence(axon_idx, day) = ...
                    mean([axon_fluorescence(axon_idx, day), tbl_struct(i).Mean]);
            end
        end
    end

    % Compute per-mouse average
    mouse_fluorescence = NaN(length(mouse_ids), num_days);
    for m = 1:length(mouse_ids)
        m_id = mouse_ids(m);
        m_idx = find([tbl_struct.Mouse] == m_id);
        for d = 1:num_days
            values = [tbl_struct(m_idx).Day] == d;
            day_vals = [tbl_struct(m_idx(values)).Mean];
            if ~isempty(day_vals)
                mouse_fluorescence(m, d) = mean(day_vals, 'omitnan');
            end
        end
    end

    % Store results
    mouse_avg_per_group.(group_names{g}) = mouse_fluorescence;
    axon_avg_per_group.(group_names{g}) = axon_fluorescence;
    synapse_per_axon.(group_names{g}) = synapse_matrix;
    synapse_per_axon_ID.(group_names{g}) = Synapse_ID;
end

%%
group_names = {'adolescent', 'adult'};
num_days = 4;
TOR_values_group = cell(1, 2);  % per-axon TORs for each group

% === Compute TOR per axon for each group ===
for g = 1:2
    group = group_names{g};
    syn_matrix = synapse_per_axon.(group);       % synapse x day
    syn_axon_ids = synapse_per_axon_ID.(group);  % synapse x day

    % Get unique axon IDs
    axon_ids = unique(syn_axon_ids(:));
    axon_ids(isnan(axon_ids)) = [];
    num_axons = length(axon_ids);
    TOR_per_axon = NaN(num_axons, num_days - 1);  % 3 transitions

    % Compute TOR per axon based on synapse presence
    for a = 1:num_axons
        ax_id = axon_ids(a);
        syn_idx = any(syn_axon_ids == ax_id, 2);  % synapses belonging to this axon

        for d = 1:3
            present_d1 = ~isnan(syn_matrix(syn_idx, d));
            present_d2 = ~isnan(syn_matrix(syn_idx, d+1));

            lost = sum(present_d1 & ~present_d2);
            gained = sum(~present_d1 & present_d2);
            total = sum(present_d1) + sum(present_d2);

            if total > 0
                TOR_per_axon(a, d) = (lost + gained) / total;
            end
        end
    end

    TOR_values_group{g} = TOR_per_axon;
end

% === Statistical comparison per TOR (per-axon basis) ===
p = NaN(1, 3); h = NaN(1, 3);
mean_tor = NaN(2, 3);  % rows: group 1/2, cols: TOR1-3
sem_tor = NaN(2, 3);

for d = 1:3
    tor1 = TOR_values_group{1}(:, d);  % adolescent
    tor2 = TOR_values_group{2}(:, d);  % adult

    tor1(tor1 == 0) = nan ; 
    tor2(tor2 == 0) = nan ; 

    % Store mean and SEM
    mean_tor(1, d) = nanmean(tor1);
    mean_tor(2, d) = nanmean(tor2);
    sem_tor(1, d) = nanstd(tor1) / sqrt(length(tor1));
    sem_tor(2, d) = nanstd(tor2) / sqrt(length(tor2));

    % Perform two-sample t-test
    [p(d), ~, ~] = ranksum(tor1, tor2);
end
figure
errorbar([1 3 5], mean_tor(1, :),sem_tor(1, :),'linestyle','none','Color',colors_ado_adu{1});
hold on
bar([1 3 5], mean_tor(1, :),'linestyle','none','FaceColor',colors_ado_adu{1},'EdgeColor','none','BarWidth', 0.4);
hold on
errorbar([2 4 6], mean_tor(2, :),sem_tor(2, :),'linestyle','none','Color',colors_ado_adu{2});
hold on
bar([2 4 6], mean_tor(2, :),'linestyle','none','FaceColor',colors_ado_adu{2},'EdgeColor','none','BarWidth', 0.4);
box off
ylabel('Turnover Rate (%)');
xlim([0 7])
ylim([0 1])


%% lost stable and gained bouttons 


% === Compute overall stable/lost/gained ratios per synapse, averaged across 3 transitions ===
stable_ratios = cell(1,2);
lost_ratios = cell(1,2);
gained_ratios = cell(1,2);

for g = 1:2
    syn_matrix = synapse_per_axon.(group_names{g});
    presence = ~isnan(syn_matrix);  % [syn x day]
    n_syn = size(syn_matrix, 1);

    N_stable = zeros(n_syn, 1);
    N_lost   = zeros(n_syn, 1);
    N_gained = zeros(n_syn, 1);

    for d = 1:3
        pre = presence(:, d);
        post = presence(:, d+1);

        N_stable = N_stable + (pre & post);
        N_lost   = N_lost   + (pre & ~post);
        N_gained = N_gained + (~pre & post);
    end

    % Compute per-synapse ratios across 3 days
    stable_ratios{g} = N_stable ./ (N_stable + N_lost + eps);      % +eps avoids /0
    lost_ratios{g}   = N_lost   ./ (N_stable + N_lost + eps);
    gained_ratios{g} = N_gained ./ (N_stable + N_gained + eps);
end

% Drop NaNs
s1 = stable_ratios{1}(~isnan(stable_ratios{1}));
s2 = stable_ratios{2}(~isnan(stable_ratios{2}));

l1 = lost_ratios{1}(~isnan(lost_ratios{1}));
l2 = lost_ratios{2}(~isnan(lost_ratios{2}));

g1 = gained_ratios{1}(~isnan(gained_ratios{1}));
g2 = gained_ratios{2}(~isnan(gained_ratios{2}));

% Mean ± SEM
mean_stable = [mean(s1), mean(s2)];
sem_stable  = [std(s1)/sqrt(length(s1)), std(s2)/sqrt(length(s2))];

mean_lost = [mean(l1), mean(l2)];
sem_lost  = [std(l1)/sqrt(length(l1)), std(l2)/sqrt(length(l2))];

mean_gained = [mean(g1), mean(g2)];
sem_gained  = [std(g1)/sqrt(length(g1)), std(g2)/sqrt(length(g2))];

% t-tests
[~, p_stable] = ttest2(s1, s2);
[~, p_lost]   = ttest2(l1, l2);
[~, p_gained] = ttest2(g1, g2);


figure
errorbar([1], mean_stable(1) ,sem_stable(1),'linestyle','none','Color',colors_ado_adu{1});
hold on
bar([1], mean_stable(1) ,'linestyle','none','FaceColor',colors_ado_adu{1},'EdgeColor','none','BarWidth', 0.8);
hold on
errorbar([2], mean_stable(2) ,sem_stable(2),'linestyle','none','Color',colors_ado_adu{2});
hold on
bar([2], mean_stable(2) ,'linestyle','none','FaceColor',colors_ado_adu{2},'EdgeColor','none','BarWidth', 0.8);
hold on

errorbar([3], mean_lost(1) ,sem_lost(1),'linestyle','none','Color',colors_ado_adu{1});
hold on
bar([3], mean_lost(1) ,'linestyle','none','FaceColor',colors_ado_adu{1},'EdgeColor','none','BarWidth', 0.8);
hold on
errorbar([4], mean_lost(2) ,sem_lost(2),'linestyle','none','Color',colors_ado_adu{2});
hold on
bar([4], mean_lost(2) ,'linestyle','none','FaceColor',colors_ado_adu{2},'EdgeColor','none','BarWidth', 0.8);
hold on

errorbar([5], mean_gained(1) ,sem_gained(1),'linestyle','none','Color',colors_ado_adu{1});
hold on
bar([5], mean_gained(1) ,'linestyle','none','FaceColor',colors_ado_adu{1},'EdgeColor','none','BarWidth', 0.8);
hold on
errorbar([6], mean_gained(2) ,sem_gained(2),'linestyle','none','Color',colors_ado_adu{2});
hold on
bar([6], mean_gained(2) ,'linestyle','none','FaceColor',colors_ado_adu{2},'EdgeColor','none','BarWidth', 0.8);
hold on
box off
ylabel('% of spines');
xlim([0 7])
ylim([0 1])



%% CV of mean density
for g = 1:2
        synapse_IntDen_cell{1, g}(isnan(synapse_IntDen_cell{1, g})) = 0 ;

    var_fluo{1,g} = nan (2500,1) ;
    var_fluo_vec = nanstd(synapse_IntDen_cell{1, g}')' ./ nanmean(synapse_IntDen_cell{1, g}')' ;
    var_fluo{1,g}(1:size( var_fluo_vec,1),1) = var_fluo_vec ;
     var_fluo{1,g} = double( var_fluo{1,g}) ;
     var_fluo{1,g}(var_fluo{1,g} == 2) = nan ; 
     
end

  figure
violinplot([var_fluo{1,1}, var_fluo{1,2}],{'adolescent ','adult'},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;
box off
ylabel('CV of density')
box off;
hold on;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

[p,h,stats] = ranksum(var_fluo{1,1},var_fluo{1,2})

%% CV of area
for g = 1:2
        %synapse_area_cell{1, g}(isnan(synapse_area_cell{1, g})) = 0 ;
    var_fluo{1,g} = nan (2500,1) ;
    var_fluo_vec = nanstd(synapse_area_cell{1, g}')' ./ nanmean(synapse_area_cell{1, g}')' ;
    var_fluo{1,g}(1:size( var_fluo_vec,1),1) = var_fluo_vec ;
     var_fluo{1,g} = double( var_fluo{1,g}) ;
    % var_fluo{1,g}(var_fluo{1,g} >= 2) = nan ; 
     
end

  figure
violinplot([var_fluo{1,1}, var_fluo{1,2}],{'adolescent ','adult'},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;
box off
ylabel('CV of area')

box off;
hold on;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

[p,h,stats] = ranksum(var_fluo{1,1},var_fluo{1,2})

%% CV of mean fluorescence
for g = 1:2
        synapse_matrix_cell{1, g}(isnan(synapse_matrix_cell{1, g})) = 0 ;
    var_fluo{1,g} = nan (2500,1) ;
    var_fluo_vec = nanstd(synapse_matrix_cell{1, g}')' ./ nanmean(synapse_matrix_cell{1, g}')' ;
    var_fluo{1,g}(1:size( var_fluo_vec,1),1) = var_fluo_vec ;
    var_fluo{1,g} = double( var_fluo{1,g}) ;
    var_fluo{1,g}(var_fluo{1,g} >= 2) = nan ;

end

figure
violinplot([var_fluo{1,1}, var_fluo{1,2}],{'adolescent ','adult'},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;
box off
ylabel('CV of mean fluorescence')
ylim([0 2])
box off;
hold on;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');


%%
% Extract tables
T1 = masterTable_cell{1,1};
T2 = masterTable_cell{1,2};

% Add group labels
T1.Group = repmat("Group1", height(T1), 1);
T2.Group = repmat("Group2", height(T2), 1);


% Combine tables
T = [T1; T2];

% Convert categorical variables if needed
T.Mouse = categorical(T.Mouse);
T.Axon = categorical(T.Axon);
T.Day = categorical(T.Day);
T.Type = categorical(T.Type);
T.Group = categorical(T.Group);

% Fit Linear Mixed Effects model
% Replace 'YourVariable' with the name of the variable you want to compare

lme = fitlme(T, 'Area ~ Group + (1|Mouse) + (1|Axon) + (1|Day) + (1|Type)') ;

lme = fitlme(T, 'IntDen ~ Group + (1|Mouse) + (1|Axon) + (1|Day) + (1|Type)') ;

lme = fitlme(T, 'Mean ~ Group + (1|Mouse) + (1|Axon) + (1|Day) + (1|Type)') ;
%%


clearvars -except masterTable_cell

colors_ado_adu = {[.5 .7 .2],[ .2 .4 .2]} ; 

%%
group_names = {'adolescent', 'adult'};
mice_counts = [3, 4];
num_days = 4;

synapse_per_day = struct();
axon_avg_per_group = struct();
synapse_per_axon = struct();
mouse_avg_per_group = struct();

for g = 1:2
    tbl = masterTable_cell{1, g};

    if istable(tbl)
        tbl_struct = table2struct(tbl);
    else
        tbl_struct = tbl;
    end

    mouse_ids = unique([tbl_struct.Mouse]);
    axon_id_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    synapse_id_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    synapse_counter = 1;
    global_axon_id = 1;  % Unique axon ID counter across all mice

    % Assign new axon IDs across mice
    for m = 1:length(mouse_ids)
        mouse = mouse_ids(m);
        mouse_idx = find([tbl_struct.Mouse] == mouse);
        axons = unique([tbl_struct(mouse_idx).Axon]);

        for a = 1:length(axons)
            original_id = axons(a);
            key = sprintf('M%d_A%d', mouse, original_id);
            axon_id_map(key) = global_axon_id;
            global_axon_id = global_axon_id + 1;
        end
    end

    % Assign new axon and synapse IDs to each entry
    new_axon_ids = zeros(length(tbl_struct), 1);
    new_synapse_ids = zeros(length(tbl_struct), 1);
    for i = 1:length(tbl_struct)
        syn = tbl_struct(i);
        key_ax = sprintf('M%d_A%d', syn.Mouse, syn.Axon);
        key_syn = sprintf('M%d_A%d_S%d', syn.Mouse, syn.Axon, syn.Var1);

        new_axon_ids(i) = axon_id_map(key_ax);

        if ~isKey(synapse_id_map, key_syn)
            synapse_id_map(key_syn) = synapse_counter;
            synapse_counter = synapse_counter + 1;
        end
        new_synapse_ids(i) = synapse_id_map(key_syn);
    end

    % Build axon fluorescence matrix
    total_axons = global_axon_id - 1;
    axon_fluorescence = NaN(total_axons, num_days);

    % Build synapse matrix
    total_synapses = synapse_counter - 1;
    synapse_matrix = NaN(total_synapses, num_days);
    synapse_area = NaN(total_synapses, num_days);
    synapse_max = NaN(total_synapses, num_days);
    synapse_IntDen = NaN(total_synapses, num_days);
    Synapse_ID = NaN(total_synapses, num_days);  % Initialize

    for i = 1:length(tbl_struct)
        syn = tbl_struct(i);
        day = syn.Day;
        syn_idx = new_synapse_ids(i);
        axon_idx = new_axon_ids(i);

        if day >= 1 && day <= num_days
            synapse_matrix(syn_idx, day) = syn.Mean;
            Synapse_ID(syn_idx, day) = axon_idx;  % use updated axon ID
            synapse_area(syn_idx, day) = syn.Area;
            synapse_max(syn_idx, day) = syn.Max;
            synapse_IntDen(syn_idx, day) = syn.IntDen;
        end
    end

    synapse_per_day.(group_names{g}) = synapse_matrix;
    synapse_matrix_cell{1,g} = synapse_matrix;
    synapse_axon_cell{1,g} = Synapse_ID;
    synapse_area_cell{1,g} = synapse_area;
    synapse_max_cell{1,g} = synapse_max;
    synapse_IntDen_cell{1,g} = synapse_IntDen;

    % Compute axon fluorescence average
    for i = 1:length(tbl_struct)
        axon_idx = new_axon_ids(i);
        day = tbl_struct(i).Day;
        if day >= 1 && day <= num_days
            if isnan(axon_fluorescence(axon_idx, day))
                axon_fluorescence(axon_idx, day) = tbl_struct(i).Mean;
            else
                axon_fluorescence(axon_idx, day) = ...
                    mean([axon_fluorescence(axon_idx, day), tbl_struct(i).Mean]);
            end
        end
    end

    % Compute per-mouse average
    mouse_fluorescence = NaN(length(mouse_ids), num_days);
    for m = 1:length(mouse_ids)
        m_id = mouse_ids(m);
        m_idx = find([tbl_struct.Mouse] == m_id);
        for d = 1:num_days
            values = [tbl_struct(m_idx).Day] == d;
            day_vals = [tbl_struct(m_idx(values)).Mean];
            if ~isempty(day_vals)
                mouse_fluorescence(m, d) = mean(day_vals, 'omitnan');
            end
        end
    end

    % Store results
    mouse_avg_per_group.(group_names{g}) = mouse_fluorescence;
    axon_avg_per_group.(group_names{g}) = axon_fluorescence;
    synapse_per_axon.(group_names{g}) = synapse_matrix;
    synapse_per_axon_ID.(group_names{g}) = Synapse_ID;
end
%% stable lost and gained synapses

% === Compute overall stable/lost/gained ratios per synapse, averaged across 3 transitions ===
stable_ratios = cell(1,2);
lost_ratios = cell(1,2);
gained_ratios = cell(1,2);

for g = 1:2
    syn_matrix = synapse_per_axon.(group_names{g});
    presence = ~isnan(syn_matrix);  % [syn x day]
    n_syn = size(syn_matrix, 1);

    N_stable = zeros(n_syn, 1);
    N_lost   = zeros(n_syn, 1);
    N_gained = zeros(n_syn, 1);

    for d = 1:3
        pre = presence(:, d);
        post = presence(:, d+1);

        N_stable = N_stable + (pre & post);
        N_lost   = N_lost   + (pre & ~post);
        N_gained = N_gained + (~pre & post);
    end

    % Compute per-synapse ratios across 3 days
    stable_ratios{g} = N_stable ./ (N_stable + N_lost + eps);      % +eps avoids /0
    lost_ratios{g}   = N_lost   ./ (N_stable + N_lost + eps);
    gained_ratios{g} = N_gained ./ (N_stable + N_gained + eps);
end

%%

% Drop NaNs
s1 = stable_ratios{1}(~isnan(stable_ratios{1}));
s2 = stable_ratios{2}(~isnan(stable_ratios{2}));

l1 = lost_ratios{1}(~isnan(lost_ratios{1}));
l2 = lost_ratios{2}(~isnan(lost_ratios{2}));

g1 = gained_ratios{1}(~isnan(gained_ratios{1}));
g2 = gained_ratios{2}(~isnan(gained_ratios{2}));

% Mean ± SEM
mean_stable = [mean(s1), mean(s2)];
sem_stable  = [std(s1)/sqrt(length(s1)), std(s2)/sqrt(length(s2))];

mean_lost = [mean(l1), mean(l2)];
sem_lost  = [std(l1)/sqrt(length(l1)), std(l2)/sqrt(length(l2))];

mean_gained = [mean(g1), mean(g2)];
sem_gained  = [std(g1)/sqrt(length(g1)), std(g2)/sqrt(length(g2))];

% t-tests
[~, p_stable] = ttest2(s1, s2);
[~, p_lost]   = ttest2(l1, l2);
[~, p_gained] = ttest2(g1, g2);

%%
% === Define data ===
means = [
    mean_stable(1), mean_stable(2);
    mean_lost(1),   mean_lost(2);
    mean_gained(1), mean_gained(2)
];  % rows = type (stable/lost/gained), columns = group (ado/adu)

sems = [
    sem_stable(1), sem_stable(2);
    sem_lost(1),   sem_lost(2);
    sem_gained(1), sem_gained(2)
];

% === Colors ===
colors_ado_adu = {[0.5 0.7 0.2], [0.2 0.4 0.2]};

% === Create bar plot ===
figure;
hold on;
b = bar(means, 'grouped', 'BarWidth', 0.6);

% Apply custom colors
for k = 1:2
    b(k).FaceColor = colors_ado_adu{k};
end

% === Add error bars ===
ngroups = size(means, 1);  % 3: stable, lost, gained
nbars = size(means, 2);    % 2: adolescent, adult
groupwidth = min(0.8, nbars/(nbars + 1.5));

for i = 1:ngroups
    for j = 1:nbars
        % Center each bar
        x = i - groupwidth/2 + (2*j-1) * groupwidth / (2*nbars);
        errorbar(x, means(i,j), sems(i,j), 'k', 'linestyle', 'none', 'LineWidth', 1.2);
    end
end

% === Axis formatting ===
set(gca, 'XTick', 1:3, 'XTickLabel', {'Stable', 'Lost', 'Gained'});
ylabel('Fraction of Synapses');
legend({'Adolescent', 'Adult'}, 'Location', 'northeast');
ylim([0, 1]);  % or adjust as needed
box off;
hold off;
%% TOR 
% === Parameters ===
group_names = {'adolescent', 'adult'};
num_days = 4;
TOR_values_group = cell(1, 2);  % TOR per axon per group

% === Compute TOR per axon ===
for g = 1:2
    group = group_names{g};

    syn_matrix = synapse_per_axon.(group);       % synapse x day
    syn_axon_ids = synapse_per_axon_ID.(group);  % synapse x day

    % Get unique axons from synapse ID matrix
    all_axons = unique(syn_axon_ids(:));
    all_axons(isnan(all_axons)) = [];
    num_axons = length(all_axons);

    TOR_per_axon = NaN(num_axons, num_days - 1);  % 3 transitions

    for a = 1:num_axons
        axon_id = all_axons(a);
        syn_idx = any(syn_axon_ids == axon_id, 2);  % synapses belonging to this axon
        syn_present = ~isnan(syn_matrix(syn_idx, :));  % presence matrix: syn x day

        for d = 1:3
            present_d1 = syn_present(:, d);
            present_d2 = syn_present(:, d+1);

            lost   = sum(present_d1 & ~present_d2);
            gained = sum(~present_d1 & present_d2);
            total  = sum(present_d1) + sum(present_d2);

            if total > 0
                TOR_per_axon(a, d) = (lost + gained) / total;
            end
        end
    end

    TOR_values_group{g} = TOR_per_axon;  % store for stats
end
%%
% Initialize outputs
p = NaN(1, 3); h = NaN(1, 3);
mean_tor = NaN(2, 3);  % rows = group, cols = D1–D2, D2–D3, D3–D4
sem_tor = NaN(2, 3);

% Compare TOR per axon across groups
for d = 1:3
    tor_ado = TOR_values_group{1}(:, d);  % adolescent
    tor_adu = TOR_values_group{2}(:, d);  % adult

    tor_ado(tor_ado == 0) = nan;
    tor_adu(tor_adu == 0) = nan;

    mean_tor(1, d) = nanmean(tor_ado);
    mean_tor(2, d) = nanmean(tor_adu);

    sem_tor(1, d) = nanstd(tor_ado) / sqrt(length(tor_ado));
    sem_tor(2, d) = nanstd(tor_adu) / sqrt(length(tor_adu));

    [h(d), p(d)] = ttest2(tor_ado, tor_adu);  % significance test
end
%%
% === Define TOR data (already computed) ===
% mean_tor: 2x3 [group, day_transition]
% sem_tor:  2x3

% Transitions: D1–D2, D2–D3, D3–D4
tor_labels = {'Day 1–2', 'Day 2–3', 'Day 3–4'};

% Bar chart data: rows = transitions, cols = groups
means = mean_tor';  % 3x2: [transitions x group]
sems  = sem_tor';

% === Create bar chart ===
figure;
hold on;
b = bar(means, 'grouped', 'BarWidth', 0.6);

% Apply colors
colors_ado_adu = {[0.5 0.7 0.2], [0.2 0.4 0.2]};
for k = 1:2
    b(k).FaceColor = colors_ado_adu{k};
end

% Add error bars
ngroups = size(means, 1);  % 3 transitions
nbars = size(means, 2);    % 2 groups
groupwidth = min(0.8, nbars / (nbars + 1.5));

for i = 1:ngroups
    for j = 1:nbars
        x = i - groupwidth/2 + (2*j - 1) * groupwidth / (2*nbars);
        errorbar(x, means(i, j), sems(i, j), 'k', 'linestyle', 'none', 'LineWidth', 1.2);
    end
end
ylim([0 1])
box off
ylabel('TOR')

%%
synapse_density_per_axon = cell(1,2);  % {ado, adu}

for g = 1:2
    tbl = masterTable_cell{1,g};
    if istable(tbl)
        tbl_struct = table2struct(tbl);
    else
        tbl_struct = tbl;
    end

    % Identify axons
    mouse_ids = unique([tbl_struct.Mouse]);
    global_axon_id = 1;
    axon_id_map = containers.Map('KeyType', 'char', 'ValueType', 'double');

    % === Assign unique axon ID per mouse+axon
    for m = 1:length(mouse_ids)
        m_id = mouse_ids(m);
        idx = find([tbl_struct.Mouse] == m_id);
        axons = unique([tbl_struct(idx).Axon]);
        for a = 1:length(axons)
            key = sprintf('M%d_A%d', m_id, axons(a));
            axon_id_map(key) = global_axon_id;
            global_axon_id = global_axon_id + 1;
        end
    end

    num_axons = axon_id_map.Count;
    density_matrix = NaN(num_axons, num_days);

    for i = 1:length(tbl_struct)
        row = tbl_struct(i);
        if row.Type ~= 1 && row.Type ~= 2
            continue
        end
        key = sprintf('M%d_A%d', row.Mouse, row.Axon);
        ax_id = axon_id_map(key);
        day = row.Day;
        if day < 1 || day > num_days
            continue
        end

        % Use sum of Mean for synapses; Mean for axon
        if row.Type == 1
            % Sum of fluorescence over all synapses
            if isnan(density_matrix(ax_id, day))
                density_matrix(ax_id, day) = row.Mean;
            else
                density_matrix(ax_id, day) = density_matrix(ax_id, day) + row.Mean;
            end
        elseif row.Type == 2
            % Divide synapse sum by axon mean
            % Store denominator in separate matrix
            if ~exist('axon_denominator', 'var')
                axon_denominator = NaN(num_axons, num_days);
            end
            axon_denominator(ax_id, day) = row.Mean;
        end
    end

    % Normalize synapse sum by axon segment intensity
    density_matrix = density_matrix ./ (axon_denominator + eps);
    synapse_density_per_axon{g} = density_matrix;
end
%%
syn_density_matrix = cell(1,2);  % per group

for g = 1:2
    tbl = masterTable_cell{1,g};
    if istable(tbl)
        tbl_struct = table2struct(tbl);
    else
        tbl_struct = tbl;
    end

    % Build axon ID map
    mouse_ids = unique([tbl_struct.Mouse]);
    global_axon_id = 1;
    axon_id_map = containers.Map('KeyType', 'char', 'ValueType', 'double');

    for m = 1:length(mouse_ids)
        m_id = mouse_ids(m);
        idx = find([tbl_struct.Mouse] == m_id);
        axons = unique([tbl_struct(idx).Axon]);
        for a = 1:length(axons)
            key = sprintf('M%d_A%d', m_id, axons(a));
            axon_id_map(key) = global_axon_id;
            global_axon_id = global_axon_id + 1;
        end
    end

    num_axons = axon_id_map.Count;
    density_matrix = NaN(num_axons, num_days);
    axon_volume = NaN(num_axons, num_days);
    bouton_count = zeros(num_axons, num_days);

    % Accumulate data
    for i = 1:length(tbl_struct)
        row = tbl_struct(i);
        if row.Type ~= 1 && row.Type ~= 2
            continue
        end
        key = sprintf('M%d_A%d', row.Mouse, row.Axon);
        ax_id = axon_id_map(key);
        day = row.Day;
        if day < 1 || day > num_days
            continue
        end

        if row.Type == 1
            bouton_count(ax_id, day) = bouton_count(ax_id, day) + 1;
        elseif row.Type == 2
            axon_volume(ax_id, day) = row.Mean;  % proxy for axon length
        end
    end

    % Compute density: boutons per axon unit
    density_matrix = bouton_count ./ (axon_volume + eps);
    syn_density_matrix{g} = density_matrix;
end
%%
density_TOR_group = cell(1,2);

for g = 1:2
    density = syn_density_matrix{g};  % [axon x day]
    n_axons = size(density, 1);
    density_TOR = NaN(n_axons, 3);  % 3 transitions

    for d = 1:3
        rho_t = density(:, d);
        rho_tp1 = density(:, d+1);

        valid = ~isnan(rho_t) & ~isnan(rho_tp1);
        numerator = abs(rho_tp1(valid) - rho_t(valid));
        denominator = rho_tp1(valid) + rho_t(valid);

        density_TOR(valid, d) = numerator ./ (denominator + eps);
    end

    density_TOR_group{g} = density_TOR;
end

%%
mean_dtor = NaN(2,3);
sem_dtor = NaN(2,3);
p_dtor = NaN(1,3);

for d = 1:3
    ado = density_TOR_group{1}(:, d);
    adu = density_TOR_group{2}(:, d);

    ado = ado(~isnan(ado));
    adu = adu(~isnan(adu));

    mean_dtor(1,d) = mean(ado);
    mean_dtor(2,d) = mean(adu);

    sem_dtor(1,d) = std(ado)/sqrt(length(ado));
    sem_dtor(2,d) = std(adu)/sqrt(length(adu));

    [~, p_dtor(d)] = ttest2(ado, adu);
end
%%
% Data
means = mean_dtor';
sems = sem_dtor';
labels = {'D1–D2', 'D2–D3', 'D3–D4'};
colors_ado_adu = {[0.5 0.7 0.2], [0.2 0.4 0.2]};

% Plot
figure;
hold on;
b = bar(means, 'grouped', 'BarWidth', 0.6);

% Colors
for k = 1:2
    b(k).FaceColor = colors_ado_adu{k};
end

% Error bars
ngroups = size(means, 1);
nbars = size(means, 2);
groupwidth = min(0.8, nbars / (nbars + 1.5));

for i = 1:ngroups
    for j = 1:nbars
        x = i - groupwidth/2 + (2*j - 1) * groupwidth / (2*nbars);
        errorbar(x, means(i,j), sems(i,j), 'k', 'linestyle', 'none', 'LineWidth', 1.2);
    end
end

% Final styling
set(gca, 'XTick', 1:3, 'XTickLabel', labels);
ylim([0 1])
ylabel('Synapse Density TOR (Boutons / μm proxy)');
legend({'Adolescent', 'Adult'}, 'Location', 'northeast');
title('Density-Based Turnover Ratio of Synapses');
box off;
hold off;

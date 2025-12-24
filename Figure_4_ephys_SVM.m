%%% - this script includes shuffleled data - graphs might slightley differ
%%% from figures 
% Clear workspace and keep relevant variable
clearvars -except GNG_rec_all_cell
close all
clc
%% Define parameters
startRange = -0.2;
stopRange = 0.6;
binSize = 0.001;
smoothSize = 20;
min_spikes = 100;
rasterScale = 5;
colors_ado_adu = {[.5 .7 .2],[ .2 .4 .2]} ;
window = [-0.2; 0.6]; % time window according to stimulus onset
directory = 'Z:\Shared\Benne\Praegel_et_al_2025\figures';

%% extract event_times and spike_times of all units
clc
tic
for g   = 1:numel(GNG_rec_all_cell) % run per group
    GNG_rec_all = cell2table(GNG_rec_all_cell(g)); %convert to table
    GNG_rec_all = GNG_rec_all.Var1{1,1};

    Recs = unique(string(1:height(GNG_rec_all)),'stable'); % number of recordings

    for i = 1:length(Recs) % run per recording

        et{g,i} = GNG_rec_all(i).eventTimes
        spikes = GNG_rec_all(i).Units; %extract spikes
        spikeTimes = spikes.st; %extract spike times
        clusters = spikes.clu; %extract cluster ids per spike
        clusterIDs = unique(GNG_rec_all(i).Units.Unit_table.clu); % extract the unique clusterIDs
        clusterIDs = clusterIDs(clusterIDs >0);


        st_pop{g,i} = spikeTimes;

        % run this for loop to determine significantly modulated units
        % across all events.
        for c = 1:length(clusterIDs)  %run per cluster = per neuron
            st{g,i}{c,1} = spikeTimes(clusters == clusterIDs(c));
        end
    end
end
toc
%% create FR_array_matrix
clc
tic

% test if neurons are modulated by the tone
for g   = 1:size(et,1) % run per group

    for i = 1:size(et,2) % run per recording

        spikeTimes = st{g,i}  ; %extract spike times


        et_single = et{g, i};
        % Extract indices for different trial types
        hit_easy_idx = et_single.tr_resp_idx(1, ~isnan(et_single.tr_resp_idx(1, :)))';
        fa_easy_idx = et_single.tr_resp_idx(2, ~isnan(et_single.tr_resp_idx(2, :)))';
        miss_easy_idx = et_single.tr_resp_idx(3, ~isnan(et_single.tr_resp_idx(3, :)))';
        cr_easy_idx = et_single.tr_resp_idx(4, ~isnan(et_single.tr_resp_idx(4, :)))';

        hit_hard_idx = et_single.tr_resp_idx(7, ~isnan(et_single.tr_resp_idx(7, :)))';
        fa_hard_idx = et_single.tr_resp_idx(8, ~isnan(et_single.tr_resp_idx(8, :)))';
        miss_hard_idx = et_single.tr_resp_idx(9, ~isnan(et_single.tr_resp_idx(9, :)))';
        cr_hard_idx = et_single.tr_resp_idx(10, ~isnan(et_single.tr_resp_idx(10, :)))';


        %------------- Lick bias -----------
        %after licks
        licks = sort([fa_easy_idx' fa_hard_idx' ]) ;
        after_licks = licks + 1 ;
        after_licks = after_licks(1:end-2) ;
        % for specific stim
        after_licks = intersect(after_licks,sort([hit_easy_idx]))' ;
        events_after_licks = et{g, i}.all(after_licks) ;
        %after no licks
        nlicks = sort([cr_easy_idx'  cr_hard_idx' ]) ;
        after_nlicks = nlicks + 1 ;
        after_nlicks = after_nlicks(1:end-2) ;
        % for specific stim
        after_nlicks = intersect(after_nlicks,sort([hit_easy_idx]))' ;
        events_after_nlicks = et{g, i}.all  (after_nlicks) ;


        for c = 1:length(st{g,i})  %run per cluster = per neuron

            %sort out spike times outside of the event time range
            spikeTime = spikeTimes{c,1}(spikeTimes{c,1}>min(events_after_licks+startRange) & spikeTimes{c,1}<max(events_after_licks+stopRange));
            if length(spikeTime) >=  min_spikes

                %extract the basics
                [~, psth,GNG_analysis.binArray{g,i}{c,1},binCenters] = GNG_binning_PSTH(startRange,stopRange, binSize, events_after_licks, spikeTime);
                % [rasterX{g,i}{c,1}, rasterY{g,i}{c,1}] = GNG_raster (GNG_analysis.binArray{g,i}{c,1}, binSize, binCenters, rasterScale);
                % extract the smoothed basics
                [~,psth_Smoothed, GNG_analysis.frArray_Smoothed{g,i}{c,1},spikeCounts] ...
                    = GNG_smoothed_PSTH (startRange, stopRange, binCenters,GNG_analysis.binArray{g,i}{c,1}, binSize, smoothSize);

                psth_lick = psth_Smoothed (1,200:400) ;
                psth_lick_all{g,i}(c,:) = psth_Smoothed ;
            end


            %sort out spike times outside of the event time range
            spikeTime = spikeTimes{c,1}(spikeTimes{c,1}>min(events_after_nlicks+startRange) & spikeTimes{c,1}<max(events_after_nlicks+stopRange));
            if length(spikeTime) >=  min_spikes

                %extract the basics
                [~, psth,GNG_analysis.binArray{g,i}{c,1},binCenters] = GNG_binning_PSTH(startRange,stopRange, binSize, events_after_nlicks, spikeTime);
                % [rasterX{g,i}{c,1}, rasterY{g,i}{c,1}] = GNG_raster (GNG_analysis.binArray{g,i}{c,1}, binSize, binCenters, rasterScale);
                % extract the smoothed basics
                [~,psth_Smoothed, GNG_analysis.frArray_Smoothed{g,i}{c,1},spikeCounts] ...
                    = GNG_smoothed_PSTH (startRange, stopRange, binCenters,GNG_analysis.binArray{g,i}{c,1}, binSize, smoothSize);

                psth_nlick = psth_Smoothed (1,200:400) ;
                                psth_nlick_all{g,i}(c,:) = psth_Smoothed ;

            end


            %look at both excited and inhibited neurons or population decoding
            [p,~,~] = ranksum(psth_lick , psth_nlick,0.05,'tail','both');
            if  p < 0.05 % sort out non-excited units
                clusterIDs_modh{g,i}(c,:) = c ;
            end

        end


        cn =  clusterIDs_modh{g,i}(clusterIDs_modh{g,i}>1)
        for cn2 = 1:length(cn)
            st_mod{g,i}{cn2,1} =  st{g,i}{cn(cn2),1} ;
        end

        st_pop{g, i} = cell2mat(st_mod{g,i}) ;
    end
end
%% psth examples
figure
plot(psth_nlick_all{g,i}(c,:),'Color',[ 0.5 0 .5])
hold on
plot(psth_lick_all{g,i}(c,:),'Color',[ 1  1 0])
xlabel('time')
xticks([0 200 400 600 800])
xticks([-200 0 200 400 600])
ylabel('Fr(Hz)')
box off 

%% percentage of modulated units


for g = 1:2
    for i = 1:7
        n_total = size(st{g,i}, 1);
        n_mod   = size(st_mod{g,i}, 1);
        mod_all(g,i,1) = n_total;
        mod_all(g,i,2) = n_mod;
    end
end

for g = 1:2
    total_all = sum(mod_all(g,:,1), 2);
    total_mod = sum(mod_all(g,:,2), 2);
    mod_group(g,1) = total_all - total_mod;
    mod_group(g,2) = total_mod;
end

% ---- Pie charts ----
figure('Name','Modified vs Not-Modified by Group');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile;
pie(mod_group(1,:), {'Not Modified','Modified'});
title('Group 1');

nexttile;
pie(mod_group(2,:), {'Not Modified','Modified'});
title('Group 2');

% ---- 2x2 Chi-square test of independence ----
O = mod_group;                         % observed
N = sum(O,'all');
assert(N > 0, 'No observations to test.');
row_sums = sum(O,2);
col_sums = sum(O,1);
E = (row_sums * col_sums) / N;         % expected
% Why: Chi-square requires nonzero expected counts.
if any(E(:) == 0)
    warning('Some expected counts are zero; chi-square may be invalid.');
end
chi2_stat = sum((O - E).^2 ./ E,'all');
df = (size(O,1)-1) * (size(O,2)-1);    % =1
p_value = 1 - chi2cdf(chi2_stat, df);

% Effect size (phi for 2x2)
phi = sqrt(chi2_stat / N);

%% SVM decoder - bias 
clc
tic

for g   = 1:size(et,1) % run per group

    for i = 1:size(et,2) % run per recording
        spikeTimes = st_pop{g, i};

        et_single = et{g, i};
        % Extract indices for different trial types
        hit_easy_idx = et_single.tr_resp_idx(1, ~isnan(et_single.tr_resp_idx(1, :)))';
        fa_easy_idx = et_single.tr_resp_idx(2, ~isnan(et_single.tr_resp_idx(2, :)))';
        miss_easy_idx = et_single.tr_resp_idx(3, ~isnan(et_single.tr_resp_idx(3, :)))';
        cr_easy_idx = et_single.tr_resp_idx(4, ~isnan(et_single.tr_resp_idx(4, :)))';

        hit_hard_idx = et_single.tr_resp_idx(7, ~isnan(et_single.tr_resp_idx(7, :)))';
        fa_hard_idx = et_single.tr_resp_idx(8, ~isnan(et_single.tr_resp_idx(8, :)))';
        miss_hard_idx = et_single.tr_resp_idx(9, ~isnan(et_single.tr_resp_idx(9, :)))';
        cr_hard_idx = et_single.tr_resp_idx(10, ~isnan(et_single.tr_resp_idx(10, :)))';


        %------------- Lick bias -----------
        %after licks
        licks = sort([fa_easy_idx' fa_hard_idx' ]) ;
        after_licks = licks + 1 ;
        after_licks = after_licks(1:end-2) ;
        % for specific stim
        after_licks = intersect(after_licks,sort([hit_easy_idx ; miss_easy_idx]))' ;
        events_after_licks = et{g, i}.all(after_licks) ;
        %after no licks
        nlicks = sort([cr_easy_idx'  cr_hard_idx' ]) ;
        after_nlicks = nlicks + 1 ;
        after_nlicks = after_nlicks(1:end-2) ;
        % for specific stim
        after_nlicks = intersect(after_nlicks,sort([hit_easy_idx; miss_easy_idx]))' ;
        events_after_nlicks = et{g, i}.all  (after_nlicks) ;


        if ~isempty(events_after_licks)
            if ~isempty(events_after_nlicks)


                %-------------------------------------
                % Filter spike times within the event time range
                spikeTimes_after_licks = spikeTimes(spikeTimes > min(events_after_licks + startRange) & spikeTimes < max(events_after_licks + stopRange));
                spikeTimes_after_nlicks = spikeTimes(spikeTimes > min(events_after_nlicks + startRange) & spikeTimes < max(events_after_nlicks + stopRange));

                % Compute binned and smoothed firing rates
                [~, ~, binArray_alicks{g, i}, binCenters] = GNG_binning_PSTH(startRange, stopRange, binSize, events_after_licks, spikeTimes_after_licks);
                [~, ~, binArray_anlicks{g, i}, ~] = GNG_binning_PSTH(startRange, stopRange, binSize, events_after_nlicks, spikeTimes_after_nlicks);

                [~, ~, frArray_Smoothed_alicks{g, i}, ~] = GNG_smoothed_PSTH(startRange, stopRange, binCenters, binArray_alicks{g, i}, binSize, smoothSize);
                [~, ~, frArray_Smoothed_anlicks{g, i}, ~] = GNG_smoothed_PSTH(startRange, stopRange, binCenters, binArray_anlicks{g, i}, binSize, smoothSize);
            end
        end

    end
end
 %% population decoding accuracy - bias 
 clc;

tic
    for g = 1:size(et, 1)
        for i = 1:size(et, 2)
            et_single = et{g, i};

                    events_A = frArray_Smoothed_alicks{g, i} ;
                    events_B = frArray_Smoothed_anlicks{g, i} ;
                    
            

            if size(events_A,1) > 10 && size(events_B,1) > 10


                min_trials = min(size(events_A, 1), size(events_B, 1));
             
                events_A = events_A(1:min_trials, :);
                events_B = events_B(1:min_trials, :);


               events_A = (events_A - mean(events_A)) ./ std(events_A);
               events_B = (events_B - mean(events_B)) ./ std(events_B);

                 baseline = events_A(:,1:200);
                events_A = (events_A - mean(mean(baseline))) ./ std(mean(baseline));

                 baseline = events_B(:,1:200);
                events_B = (events_B - mean(mean(baseline))) ./ std(mean(baseline));

                data = [events_A; events_B];
                labels = [ones(min_trials, 1); 2 * ones(min_trials, 1)];

                binsize = 50 ;
                bint_it = 1:25:800 ;

                accuracy_over_time = zeros(1, size(bint_it,2));
                accuracy_shuffled = zeros(1, size(bint_it,2));
                cv = cvpartition(labels, 'KFold', 10);

                for t = 1:length(bint_it(1:end-1))
                    data_time_bin = mean(data(:, bint_it(t):bint_it(t)+binsize-1),2);
                    accuracies = zeros(cv.NumTestSets, 1);
                    accuracies_shuffled = zeros(cv.NumTestSets, 1);

                    for fold = 1:cv.NumTestSets

                        train_idx = cv.training(fold);
                        test_idx = cv.test(fold);
                        train_data = data_time_bin(train_idx);
                        if sum(train_data) > 0
                            train_labels = labels(train_idx);
                            test_data = data_time_bin(test_idx);
                            test_labels = labels(test_idx);
                            if sum((train_labels==1)) >0
                                if  sum((train_labels==2 )) >0
                                    % SVM
                                    model = fitcsvm(train_data, train_labels, 'KernelFunction', 'rbf', 'Standardize', true) ;
                                  
                                elseif sum(train_data(train_labels==2)) == 0
                                    accuracies(fold)= nan(1,1) ;
                                end
                            elseif sum(train_data(train_labels==1)) == 0
                                accuracies(fold)= nan(1,1) ;
                            end

                            predictions = predict(model, test_data);
                            accuracies(fold) = nanmean(predictions == test_labels);


                            for fold_2 = 1:cv.NumTestSets

                                % Shuffled model
                                shuffled_labels = train_labels(randperm(length(train_labels)));
                                for perm = 1:cv.NumTestSets
                                    shuffled_labels = train_labels(randperm(length(shuffled_labels)));
                                end
                                if sum((shuffled_labels==1)) >0
                                    if  sum((shuffled_labels==2 )) >0
                                        %  %SVM
                                         model_shuffled = fitcsvm(train_data, shuffled_labels, 'KernelFunction', 'rbf', 'Standardize', true) ;
                                     
                                        predictions_shuffled = predict(model_shuffled, test_data);

                                    elseif sum(train_data(shuffled_labels==2)) == 0
                                        accuracies_shuffled_2(fold_2)= nan(1,1) ;
                                    end
                                elseif sum(train_data(shuffled_labels==1)) == 0
                                    accuracies_shuffled_2(fold_2)= nan(1,1) ;
                                end

                            accuracies_shuffled_2(fold_2) = nanmean(predictions_shuffled == test_labels);
                        end

                            accuracies_shuffled(fold) = nanmean(accuracies_shuffled_2) ; 

                        elseif sum(train_data) ==  0 
                            accuracies(fold) = nan(1,1) ;
                            accuracies_shuffled(fold) = nan(1,1) ;

                        end 
                    end
                    accuracy_over_time(t) = nanmean(accuracies);
                    accuracy_shuffled(t) = nanmean(accuracies_shuffled);
                end

                all_accuracies_over_time(i, :) = accuracy_over_time;
                all_accuracies_shuffled(i, :) = accuracy_shuffled;

                all_accuracies_over_time_cell{1, g}(i, :) = accuracy_over_time;
                all_accuracies_shuffled_cell{1, g}(i, :) = accuracy_shuffled;

            elseif size(events_A,1) <= 10  && size(events_B,1)  <= 10
                   all_accuracies_over_time_cell{1, g}(i, :) = nan(1,size(bint_it,2));
                all_accuracies_shuffled_cell{1, g}(i, :) = nan(1,size(bint_it,2));

            end
        end
    end

    for g = 1:size(et, 1)

         % Aggregate accuracy
        mean_all = smoothdata(nanmean( all_accuracies_over_time_cell{1, g}));
        stderr_all = smoothdata(nanstd( all_accuracies_over_time_cell{1, g}) / sqrt(size(et, 2)));
        mean_shuffled = smoothdata(nanmean(all_accuracies_shuffled_cell{1, g}));
        stderr_shuffled = smoothdata(nanstd(all_accuracies_shuffled_cell{1, g}) / sqrt(size(et, 2)));

        mean_all_cell{1, g} = mean_all;
        stderr_all_cell{1, g} = stderr_all;
        mean_shuffled_cell{1, g} = mean_shuffled;
        stderr_shuffled_cell{1, g} = stderr_shuffled;
    end

toc

%% absoulte decoding accuracy across time - bias

        time_axis = linspace(startRange, stopRange, size(bint_it,2));
time_axis = time_axis(1:end-4)
close all
figure
     for g = 1:2
      abs_mean_all_cell{1, g} = smoothdata(abs( mean_all_cell{1, g} -0.5)+0.5,5) ;
            abs_mean_shuffled_cell{1, g} =smoothdata(abs( mean_shuffled_cell{1, g} -0.5) +0.5,5) ;

        % Plot real vs. shuffled accuracy
        plot(time_axis(1:end-1), abs_mean_all_cell{1, g}(1:end-1), 'Color', colors_ado_adu{g}, 'LineWidth', 2);
        hold on;
        plot(time_axis(1:end-1), abs_mean_shuffled_cell{1, g}(1:end-1), '--', 'Color',  colors_ado_adu{g}, 'LineWidth', 2);

        % Shaded error for real data
        fill([time_axis(1:end-1), fliplr(time_axis(1:end-1))], ...
            [abs_mean_all_cell{1, g}(1:end-1) + stderr_all_cell{1, g}(1:end-1), fliplr(abs_mean_all_cell{1, g}(1:end-1) - stderr_all_cell{1, g}(1:end-1))], ...
            colors_ado_adu{g}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

        % Shaded error for shuffled data
        fill([time_axis(1:end-1), fliplr(time_axis(1:end-1))], ...
            [abs_mean_shuffled_cell{1, g}(1:end-1) + stderr_shuffled_cell{1, g}(1:end-1), fliplr(abs_mean_shuffled_cell{1, g}(1:end-1) - stderr_shuffled_cell{1, g}(1:end-1))], ...
            [0.5, 0.5, 0.5], 'FaceAlpha', 0.1, 'EdgeColor', 'none');

        box off;
        xlabel('Time (s)');
        ylabel('Decoding Accuracy');
        legend('Actual', 'Shuffled');
    end

    for bin = 1:size(time_axis,2)-3
        [p(bin,:),~,~] = ranksum(all_accuracies_over_time_cell{1, 1}(:,bin),...
            all_accuracies_over_time_cell{1, 2}(:,bin),'alpha',0.05,'tail','both') ;
    end

    idx_sig  = find (p < 0.05) ;
    p = p(~isnan(p)) ;
    p = - p ; % change to minus for visualization of p-value
    figure
    imagesc(p')
    xticks([])
    xline(200, 'Color','w','LineWidth',2,'LineStyle','--')
    xline(300, 'Color','w','LineWidth',2,'LineStyle','--')
    yticks([])
    h = colorbar
    h.FontSize = 20
    h.Axes.CLim = [-1 0];
    ax = gca;
    ax.XAxis.FontSize =20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;



%% SVM decoder - choice 
clc
tic

for g   = 1:size(et,1) % run per group

    for i = 1:size(et,2) % run per recording
        spikeTimes = st_pop{g, i};

        et_single = et{g, i};
        % Extract indices for different trial types
        hit_easy_idx = et_single.tr_resp_idx(1, ~isnan(et_single.tr_resp_idx(1, :)))';
        fa_easy_idx = et_single.tr_resp_idx(2, ~isnan(et_single.tr_resp_idx(2, :)))';
        miss_easy_idx = et_single.tr_resp_idx(3, ~isnan(et_single.tr_resp_idx(3, :)))';
        cr_easy_idx = et_single.tr_resp_idx(4, ~isnan(et_single.tr_resp_idx(4, :)))';

        hit_hard_idx = et_single.tr_resp_idx(7, ~isnan(et_single.tr_resp_idx(7, :)))';
        fa_hard_idx = et_single.tr_resp_idx(8, ~isnan(et_single.tr_resp_idx(8, :)))';
        miss_hard_idx = et_single.tr_resp_idx(9, ~isnan(et_single.tr_resp_idx(9, :)))';
        cr_hard_idx = et_single.tr_resp_idx(10, ~isnan(et_single.tr_resp_idx(10, :)))';


        %------------- Lick bias -----------
        %after licks
        licks = sort([fa_easy_idx' fa_hard_idx' ]) ;
        after_licks = licks  ;
        after_licks = after_licks(1:end) ;
        % for specific stim
        events_after_licks = et{g, i}.all(after_licks)
        %after no licks
        nlicks = sort([cr_easy_idx'  cr_hard_idx' ]) ;
        after_nlicks = nlicks  ;
        after_nlicks = after_nlicks(1:end) ;
        % for specific stim
        events_after_nlicks = et{g, i}.all  (after_nlicks) ;


        if ~isempty(events_after_licks)
            if ~isempty(events_after_nlicks)


                %-------------------------------------
                % Filter spike times within the event time range
                spikeTimes_after_licks = spikeTimes(spikeTimes > min(events_after_licks + startRange) & spikeTimes < max(events_after_licks + stopRange));
                spikeTimes_after_nlicks = spikeTimes(spikeTimes > min(events_after_nlicks + startRange) & spikeTimes < max(events_after_nlicks + stopRange));

                % Compute binned and smoothed firing rates
                [~, ~, binArray_alicks{g, i}, binCenters] = GNG_binning_PSTH(startRange, stopRange, binSize, events_after_licks, spikeTimes_after_licks);
                [~, ~, binArray_anlicks{g, i}, ~] = GNG_binning_PSTH(startRange, stopRange, binSize, events_after_nlicks, spikeTimes_after_nlicks);

                [~, ~, frArray_Smoothed_alicks{g, i}, ~] = GNG_smoothed_PSTH(startRange, stopRange, binCenters, binArray_alicks{g, i}, binSize, smoothSize);
                [~, ~, frArray_Smoothed_anlicks{g, i}, ~] = GNG_smoothed_PSTH(startRange, stopRange, binCenters, binArray_anlicks{g, i}, binSize, smoothSize);
            end
        end

    end
end
 %% population decoding accuracy - choice 
 clc;

tic
    for g = 1:size(et, 1)
        for i = 1:size(et, 2)
            et_single = et{g, i};

                    events_A = frArray_Smoothed_alicks{g, i} ;
                    events_B = frArray_Smoothed_anlicks{g, i} ;
                    
            

            if size(events_A,1) > 10 && size(events_B,1) > 10


                min_trials = min(size(events_A, 1), size(events_B, 1));
             
                events_A = events_A(1:min_trials, :);
                events_B = events_B(1:min_trials, :);


               events_A = (events_A - mean(events_A)) ./ std(events_A);
               events_B = (events_B - mean(events_B)) ./ std(events_B);

                 baseline = events_A(:,1:200);
                events_A = (events_A - mean(mean(baseline))) ./ std(mean(baseline));

                 baseline = events_B(:,1:200);
                events_B = (events_B - mean(mean(baseline))) ./ std(mean(baseline));

                data = [events_A; events_B];
                labels = [ones(min_trials, 1); 2 * ones(min_trials, 1)];

                binsize = 50 ;
                bint_it = 1:25:800 ;

                accuracy_over_time = zeros(1, size(bint_it,2));
                accuracy_shuffled = zeros(1, size(bint_it,2));
                cv = cvpartition(labels, 'KFold', 10);

                for t = 1:length(bint_it(1:end-1))
                    data_time_bin = mean(data(:, bint_it(t):bint_it(t)+binsize-1),2);
                    accuracies = zeros(cv.NumTestSets, 1);
                    accuracies_shuffled = zeros(cv.NumTestSets, 1);

                    for fold = 1:cv.NumTestSets

                        train_idx = cv.training(fold);
                        test_idx = cv.test(fold);
                        train_data = data_time_bin(train_idx);
                        if sum(train_data) > 0
                            train_labels = labels(train_idx);
                            test_data = data_time_bin(test_idx);
                            test_labels = labels(test_idx);
                            if sum((train_labels==1)) >0
                                if  sum((train_labels==2 )) >0
                                    % SVM
                                    model = fitcsvm(train_data, train_labels, 'KernelFunction', 'rbf', 'Standardize', true) ;
                                  
                                elseif sum(train_data(train_labels==2)) == 0
                                    accuracies(fold)= nan(1,1) ;
                                end
                            elseif sum(train_data(train_labels==1)) == 0
                                accuracies(fold)= nan(1,1) ;
                            end

                            predictions = predict(model, test_data);
                            accuracies(fold) = nanmean(predictions == test_labels);


                            for fold_2 = 1:cv.NumTestSets

                                % Shuffled model
                                shuffled_labels = train_labels(randperm(length(train_labels)));
                                for perm = 1:cv.NumTestSets
                                    shuffled_labels = train_labels(randperm(length(shuffled_labels)));
                                end
                                if sum((shuffled_labels==1)) >0
                                    if  sum((shuffled_labels==2 )) >0
                                        %  %SVM
                                         model_shuffled = fitcsvm(train_data, shuffled_labels, 'KernelFunction', 'rbf', 'Standardize', true) ;
                                     
                                        predictions_shuffled = predict(model_shuffled, test_data);

                                    elseif sum(train_data(shuffled_labels==2)) == 0
                                        accuracies_shuffled_2(fold_2)= nan(1,1) ;
                                    end
                                elseif sum(train_data(shuffled_labels==1)) == 0
                                    accuracies_shuffled_2(fold_2)= nan(1,1) ;
                                end

                            accuracies_shuffled_2(fold_2) = nanmean(predictions_shuffled == test_labels);
                        end

                            accuracies_shuffled(fold) = nanmean(accuracies_shuffled_2) ; 

                        elseif sum(train_data) ==  0 
                            accuracies(fold) = nan(1,1) ;
                            accuracies_shuffled(fold) = nan(1,1) ;

                        end 
                    end
                    accuracy_over_time(t) = nanmean(accuracies);
                    accuracy_shuffled(t) = nanmean(accuracies_shuffled);
                end

                all_accuracies_over_time(i, :) = accuracy_over_time;
                all_accuracies_shuffled(i, :) = accuracy_shuffled;

                all_accuracies_over_time_cell{1, g}(i, :) = accuracy_over_time;
                all_accuracies_shuffled_cell{1, g}(i, :) = accuracy_shuffled;

            elseif size(events_A,1) <= 10  && size(events_B,1)  <= 10
                   all_accuracies_over_time_cell{1, g}(i, :) = nan(1,size(bint_it,2));
                all_accuracies_shuffled_cell{1, g}(i, :) = nan(1,size(bint_it,2));

            end
        end
    end

    for g = 1:size(et, 1)

         % Aggregate accuracy
        mean_all = smoothdata(nanmean( all_accuracies_over_time_cell{1, g}));
        stderr_all = smoothdata(nanstd( all_accuracies_over_time_cell{1, g}) / sqrt(size(et, 2)));
        mean_shuffled = smoothdata(nanmean(all_accuracies_shuffled_cell{1, g}));
        stderr_shuffled = smoothdata(nanstd(all_accuracies_shuffled_cell{1, g}) / sqrt(size(et, 2)));

        mean_all_cell{1, g} = mean_all;
        stderr_all_cell{1, g} = stderr_all;
        mean_shuffled_cell{1, g} = mean_shuffled;
        stderr_shuffled_cell{1, g} = stderr_shuffled;
    end

toc

%% absoulte decoding accuracy across time - choice


        time_axis = linspace(startRange, stopRange, size(bint_it,2));
time_axis = time_axis(1:end-4)
close all
figure
     for g = 1:2
      abs_mean_all_cell{1, g} = smoothdata(abs( mean_all_cell{1, g} -0.5)+0.5,5) ;
            abs_mean_shuffled_cell{1, g} =smoothdata(abs( mean_shuffled_cell{1, g} -0.5) +0.5,5) ;

        % Plot real vs. shuffled accuracy
        plot(time_axis(1:end-1), abs_mean_all_cell{1, g}(1:end-1), 'Color', colors_ado_adu{g}, 'LineWidth', 2);
        hold on;
        plot(time_axis(1:end-1), abs_mean_shuffled_cell{1, g}(1:end-1), '--', 'Color',  colors_ado_adu{g}, 'LineWidth', 2);

        % Shaded error for real data
        fill([time_axis(1:end-1), fliplr(time_axis(1:end-1))], ...
            [abs_mean_all_cell{1, g}(1:end-1) + stderr_all_cell{1, g}(1:end-1), fliplr(abs_mean_all_cell{1, g}(1:end-1) - stderr_all_cell{1, g}(1:end-1))], ...
            colors_ado_adu{g}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

        % Shaded error for shuffled data
        fill([time_axis(1:end-1), fliplr(time_axis(1:end-1))], ...
            [abs_mean_shuffled_cell{1, g}(1:end-1) + stderr_shuffled_cell{1, g}(1:end-1), fliplr(abs_mean_shuffled_cell{1, g}(1:end-1) - stderr_shuffled_cell{1, g}(1:end-1))], ...
            [0.5, 0.5, 0.5], 'FaceAlpha', 0.1, 'EdgeColor', 'none');

        box off;
        xlabel('Time (s)');
        ylabel('Decoding Accuracy');
        legend('Actual', 'Shuffled');
    end

    for bin = 1:size(time_axis,2)-3
        [p(bin,:),~,~] = ranksum(all_accuracies_over_time_cell{1, 1}(:,bin),...
            all_accuracies_over_time_cell{1, 2}(:,bin),'alpha',0.05,'tail','both') ;
    end

    idx_sig  = find (p < 0.05) ;
    p = p(~isnan(p)) ;
    p = - p ; % change to minus for visualization of p-value
    figure
    imagesc(p')
    xticks([])
    xline(200, 'Color','w','LineWidth',2,'LineStyle','--')
    xline(300, 'Color','w','LineWidth',2,'LineStyle','--')
    yticks([])
    h = colorbar
    h.FontSize = 20
    h.Axes.CLim = [-1 0];
    ax = gca;
    ax.XAxis.FontSize =20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;
%% SVM decoder - ITI 
clc
tic

for g   = 1:size(et,1) % run per group

    for i = 1:size(et,2) % run per recording
        spikeTimes = st_pop{g, i};

        et_single = et{g, i};
        % Extract indices for different trial types
        hit_easy_idx = et_single.tr_resp_idx(1, ~isnan(et_single.tr_resp_idx(1, :)))';
        fa_easy_idx = et_single.tr_resp_idx(2, ~isnan(et_single.tr_resp_idx(2, :)))';
        miss_easy_idx = et_single.tr_resp_idx(3, ~isnan(et_single.tr_resp_idx(3, :)))';
        cr_easy_idx = et_single.tr_resp_idx(4, ~isnan(et_single.tr_resp_idx(4, :)))';

        hit_hard_idx = et_single.tr_resp_idx(7, ~isnan(et_single.tr_resp_idx(7, :)))';
        fa_hard_idx = et_single.tr_resp_idx(8, ~isnan(et_single.tr_resp_idx(8, :)))';
        miss_hard_idx = et_single.tr_resp_idx(9, ~isnan(et_single.tr_resp_idx(9, :)))';
        cr_hard_idx = et_single.tr_resp_idx(10, ~isnan(et_single.tr_resp_idx(10, :)))';


        %------------- Lick bias -----------
        %after licks
        licks = sort([fa_easy_idx' fa_hard_idx' ]) ;
        after_licks = licks +1 ;
        after_licks = after_licks(1:end-2) ;
        % for specific stim
        after_licks = intersect(after_licks,sort([hit_easy_idx ; miss_easy_idx]))' ;
        events_after_licks = et{g, i}.all(after_licks)-2000 ;
        %after no licks
        nlicks = sort([cr_easy_idx'  cr_hard_idx' ]) ;
        after_nlicks = nlicks +1 ;
        after_nlicks = after_nlicks(1:end-2) ;
        % for specific stim
        after_nlicks = intersect(after_nlicks,sort([hit_easy_idx; miss_easy_idx]))' ;
        events_after_nlicks = et{g, i}.all  (after_nlicks)-2000 ;


        if ~isempty(events_after_licks)
            if ~isempty(events_after_nlicks)


                %-------------------------------------
                % Filter spike times within the event time range
                spikeTimes_after_licks = spikeTimes(spikeTimes > min(events_after_licks + startRange) & spikeTimes < max(events_after_licks + stopRange));
                spikeTimes_after_nlicks = spikeTimes(spikeTimes > min(events_after_nlicks + startRange) & spikeTimes < max(events_after_nlicks + stopRange));

                % Compute binned and smoothed firing rates
                [~, ~, binArray_alicks{g, i}, binCenters] = GNG_binning_PSTH(startRange, stopRange, binSize, events_after_licks, spikeTimes_after_licks);
                [~, ~, binArray_anlicks{g, i}, ~] = GNG_binning_PSTH(startRange, stopRange, binSize, events_after_nlicks, spikeTimes_after_nlicks);

                [~, ~, frArray_Smoothed_alicks{g, i}, ~] = GNG_smoothed_PSTH(startRange, stopRange, binCenters, binArray_alicks{g, i}, binSize, smoothSize);
                [~, ~, frArray_Smoothed_anlicks{g, i}, ~] = GNG_smoothed_PSTH(startRange, stopRange, binCenters, binArray_anlicks{g, i}, binSize, smoothSize);
            end
        end

    end
end
 %% population decoding accuracy - ITI 
 clc;

tic
    for g = 1:size(et, 1)
        for i = 1:size(et, 2)
            et_single = et{g, i};

                    events_A = frArray_Smoothed_alicks{g, i} ;
                    events_B = frArray_Smoothed_anlicks{g, i} ;
                    
            

            if size(events_A,1) > 10 && size(events_B,1) > 10


                min_trials = min(size(events_A, 1), size(events_B, 1));
             
                events_A = events_A(1:min_trials, :);
                events_B = events_B(1:min_trials, :);


               events_A = (events_A - mean(events_A)) ./ std(events_A);
               events_B = (events_B - mean(events_B)) ./ std(events_B);

                 baseline = events_A(:,1:200);
                events_A = (events_A - mean(mean(baseline))) ./ std(mean(baseline));

                 baseline = events_B(:,1:200);
                events_B = (events_B - mean(mean(baseline))) ./ std(mean(baseline));

                data = [events_A; events_B];
                labels = [ones(min_trials, 1); 2 * ones(min_trials, 1)];

                binsize = 50 ;
                bint_it = 1:25:800 ;

                accuracy_over_time = zeros(1, size(bint_it,2));
                accuracy_shuffled = zeros(1, size(bint_it,2));
                cv = cvpartition(labels, 'KFold', 10);

                for t = 1:length(bint_it(1:end-1))
                    data_time_bin = mean(data(:, bint_it(t):bint_it(t)+binsize-1),2);
                    accuracies = zeros(cv.NumTestSets, 1);
                    accuracies_shuffled = zeros(cv.NumTestSets, 1);

                    for fold = 1:cv.NumTestSets

                        train_idx = cv.training(fold);
                        test_idx = cv.test(fold);
                        train_data = data_time_bin(train_idx);
                        if sum(train_data) > 0
                            train_labels = labels(train_idx);
                            test_data = data_time_bin(test_idx);
                            test_labels = labels(test_idx);
                            if sum((train_labels==1)) >0
                                if  sum((train_labels==2 )) >0
                                    % SVM
                                    model = fitcsvm(train_data, train_labels, 'KernelFunction', 'rbf', 'Standardize', true) ;
                                  
                                elseif sum(train_data(train_labels==2)) == 0
                                    accuracies(fold)= nan(1,1) ;
                                end
                            elseif sum(train_data(train_labels==1)) == 0
                                accuracies(fold)= nan(1,1) ;
                            end

                            predictions = predict(model, test_data);
                            accuracies(fold) = nanmean(predictions == test_labels);


                            for fold_2 = 1:cv.NumTestSets

                                % Shuffled model
                                shuffled_labels = train_labels(randperm(length(train_labels)));
                                for perm = 1:cv.NumTestSets
                                    shuffled_labels = train_labels(randperm(length(shuffled_labels)));
                                end
                                if sum((shuffled_labels==1)) >0
                                    if  sum((shuffled_labels==2 )) >0
                                        %  %SVM
                                         model_shuffled = fitcsvm(train_data, shuffled_labels, 'KernelFunction', 'rbf', 'Standardize', true) ;
                                     
                                        predictions_shuffled = predict(model_shuffled, test_data);

                                    elseif sum(train_data(shuffled_labels==2)) == 0
                                        accuracies_shuffled_2(fold_2)= nan(1,1) ;
                                    end
                                elseif sum(train_data(shuffled_labels==1)) == 0
                                    accuracies_shuffled_2(fold_2)= nan(1,1) ;
                                end

                            accuracies_shuffled_2(fold_2) = nanmean(predictions_shuffled == test_labels);
                        end

                            accuracies_shuffled(fold) = nanmean(accuracies_shuffled_2) ; 

                        elseif sum(train_data) ==  0 
                            accuracies(fold) = nan(1,1) ;
                            accuracies_shuffled(fold) = nan(1,1) ;

                        end 
                    end
                    accuracy_over_time(t) = nanmean(accuracies);
                    accuracy_shuffled(t) = nanmean(accuracies_shuffled);
                end

                all_accuracies_over_time(i, :) = accuracy_over_time;
                all_accuracies_shuffled(i, :) = accuracy_shuffled;

                all_accuracies_over_time_cell{1, g}(i, :) = accuracy_over_time;
                all_accuracies_shuffled_cell{1, g}(i, :) = accuracy_shuffled;

            elseif size(events_A,1) <= 10  && size(events_B,1)  <= 10
                   all_accuracies_over_time_cell{1, g}(i, :) = nan(1,size(bint_it,2));
                all_accuracies_shuffled_cell{1, g}(i, :) = nan(1,size(bint_it,2));

            end
        end
    end

    for g = 1:size(et, 1)

         % Aggregate accuracy
        mean_all = smoothdata(nanmean( all_accuracies_over_time_cell{1, g}));
        stderr_all = smoothdata(nanstd( all_accuracies_over_time_cell{1, g}) / sqrt(size(et, 2)));
        mean_shuffled = smoothdata(nanmean(all_accuracies_shuffled_cell{1, g}));
        stderr_shuffled = smoothdata(nanstd(all_accuracies_shuffled_cell{1, g}) / sqrt(size(et, 2)));

        mean_all_cell{1, g} = mean_all;
        stderr_all_cell{1, g} = stderr_all;
        mean_shuffled_cell{1, g} = mean_shuffled;
        stderr_shuffled_cell{1, g} = stderr_shuffled;
    end

toc

%% absoulte decoding accuracy across time - ITI

        time_axis = linspace(startRange, stopRange, size(bint_it,2));
time_axis = time_axis(1:end-4)
close all
figure
     for g = 1:2
      abs_mean_all_cell{1, g} = smoothdata(abs( mean_all_cell{1, g} -0.5)+0.5,5) ;
            abs_mean_shuffled_cell{1, g} =smoothdata(abs( mean_shuffled_cell{1, g} -0.5) +0.5,5) ;

        % Plot real vs. shuffled accuracy
        plot(time_axis(1:end-1), abs_mean_all_cell{1, g}(1:end-1), 'Color', colors_ado_adu{g}, 'LineWidth', 2);
        hold on;
        plot(time_axis(1:end-1), abs_mean_shuffled_cell{1, g}(1:end-1), '--', 'Color',  colors_ado_adu{g}, 'LineWidth', 2);

        % Shaded error for real data
        fill([time_axis(1:end-1), fliplr(time_axis(1:end-1))], ...
            [abs_mean_all_cell{1, g}(1:end-1) + stderr_all_cell{1, g}(1:end-1), fliplr(abs_mean_all_cell{1, g}(1:end-1) - stderr_all_cell{1, g}(1:end-1))], ...
            colors_ado_adu{g}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

        % Shaded error for shuffled data
        fill([time_axis(1:end-1), fliplr(time_axis(1:end-1))], ...
            [abs_mean_shuffled_cell{1, g}(1:end-1) + stderr_shuffled_cell{1, g}(1:end-1), fliplr(abs_mean_shuffled_cell{1, g}(1:end-1) - stderr_shuffled_cell{1, g}(1:end-1))], ...
            [0.5, 0.5, 0.5], 'FaceAlpha', 0.1, 'EdgeColor', 'none');

        box off;
        xlabel('Time (s)');
        ylabel('Decoding Accuracy');
        legend('Actual', 'Shuffled');
    end

    for bin = 1:size(time_axis,2)-3
        [p(bin,:),~,~] = ranksum(all_accuracies_over_time_cell{1, 1}(:,bin),...
            all_accuracies_over_time_cell{1, 2}(:,bin),'alpha',0.05,'tail','both') ;
    end

    idx_sig  = find (p < 0.05) ;
    p = p(~isnan(p)) ;
    p = - p ; % change to minus for visualization of p-value
    figure
    imagesc(p')
    xticks([])
    xline(200, 'Color','w','LineWidth',2,'LineStyle','--')
    xline(300, 'Color','w','LineWidth',2,'LineStyle','--')
    yticks([])
    h = colorbar
    h.FontSize = 20
    h.Axes.CLim = [-1 0];
    ax = gca;
    ax.XAxis.FontSize =20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;
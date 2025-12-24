clearvars -except GNG_rec_all_cell
%%
startRange = -0.2
stopRange = 0.6
binSize = 0.001
smoothSize = 20
min_spikes = 100 ;
rasterScale = 5

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
frArray_Smoothed_event_type_exc=  cell(size(et,1),size(et,2)) ;
frArray_Smoothed_event_type_inh=  cell(size(et,1),size(et,2)) ;

% test if neurons are modulated by the tone
for g   = 1:size(et,1) % run per group

    for i = 1:size(et,2) % run per recording

        spikeTimes = st{g,i}  ; %extract spike times

        eventTimes_all = et{g, i}.all  ;

        for c = 1:length(st{g,i})  %run per cluster = per neuron

            %sort out spike times outside of the event time range
            spikeTime = spikeTimes{c,1}(spikeTimes{c,1}>min(eventTimes_all+startRange) & spikeTimes{c,1}<max(eventTimes_all+stopRange));
            if length(spikeTime) >=  min_spikes

                %extract the basics
                [~, psth,binArray{g,i}{c,1},binCenters] = GNG_binning_PSTH(startRange,stopRange, binSize, eventTimes_all, spikeTime);

                [rasterX{g,i}{c,1}, rasterY{g,i}{c,1}] = GNG_raster (binArray{g,i}{c,1}, binSize, binCenters, rasterScale);

                % extract the smoothed basics
                [~,psth_Smoothed{g,i}(c,:), frArray_Smoothed{g,i}{c,1},spikeCounts] ...
                    = GNG_smoothed_PSTH (startRange, stopRange, binCenters,binArray{g,i}{c,1}, binSize, smoothSize);


            end
        end
    end
end
% take out the cells that didn't pass the thresholds of minimal spikes
for g   = 1:size(et,1) % run per group
    for i = 1:size(et,2) % run per recording
        frArray_Smoothed_new{g,i} = frArray_Smoothed{g,i}(~cellfun('isempty',frArray_Smoothed{g,i}))
    end
end
for g   = 1:size(et,1) % run per group
    for i = 1:size(et,2) % run per recording
        frArray_Smoothed{g,i} = frArray_Smoothed_new{g,i}
    end
end
 %% GLM - kernel regression model
clc
close all
tic
for g   = 1:size(et,1) % run per group

    for i =1:size(et,2) % run per recording

        for c = 1:height(frArray_Smoothed{g,i})
           
            % Extract Firing Rate Data for Single Neuron
            fr_rec = frArray_Smoothed{g, i}; % Firing rate matrix
            fr_single = fr_rec{c, 1}; % Firing rate for the selected neuron
            bin_rec = binArray{g, i};
            bin_single = bin_rec{c, 1};
            % Check size
            n_trials = size(fr_single, 1);
            n_bins = size(fr_single, 2);

            % Extract Event Times and Trial Types
            et_single = et{g, i};

            % Extract indices for different t (2, ~isnan(et_single.tr_resp_idx(2, :)))';
            hit_easy_idx = et_single.tr_resp_idx(1, ~isnan(et_single.tr_resp_idx(1, :)))';
            fa_easy_idx = et_single.tr_resp_idx(2, ~isnan(et_single.tr_resp_idx(2, :)))';
            miss_easy_idx = et_single.tr_resp_idx(3, ~isnan(et_single.tr_resp_idx(3, :)))';
            cr_easy_idx = et_single.tr_resp_idx(4, ~isnan(et_single.tr_resp_idx(4, :)))';

            hit_hard_idx = et_single.tr_resp_idx(7, ~isnan(et_single.tr_resp_idx(7, :)))';
            fa_hard_idx = et_single.tr_resp_idx(8, ~isnan(et_single.tr_resp_idx(8, :)))';
            miss_hard_idx = et_single.tr_resp_idx(9, ~isnan(et_single.tr_resp_idx(9, :)))';
            cr_hard_idx = et_single.tr_resp_idx(10, ~isnan(et_single.tr_resp_idx(10, :)))';


            % Extract Licks and Align to Stimulus Onset
            licks = et_single.times_vec_licks; % Lick times
            licks_aligned = cell(n_trials, 1);
            for t = 1:n_trials
                stim_onset = et_single.all(t); % Stimulus onset time for trial t (column vector)
                licks_aligned{t} = licks(licks >= stim_onset & licks <= stim_onset + n_bins * binSize) - stim_onset;
            end

            % Create Kernels Aligned to Stimulus, Action, and Bias Predictors
            easy_go_stim_kernel = [] ;
            easy_ngo_stim_kernel = [] ;
            hard_go_stim_kernel = [] ;
            hard_ngo_stim_kernel = [] ;
            lick_kernel = [] ;
            reward_kernel = [] ;
            bias_Go_kernel= [] ;
            bias_NoGo_kernel= [] ;

            for t = 1:n_trials 
                %baseline
                baseline_onset = 1;
                baseline_offset = 100 ;
                trial_duration = length(fr_single); % Trial duration in ms
                bin_size = 1;        % Bin size in ms
                num_bins = trial_duration / bin_size; % Total bins
                time_vector = linspace(0, trial_duration, num_bins);

                %stim_kernel
                stimulus_onset = 200;  % Stimulus starts at 200 ms
                stimulus_offset = 300; % Stimulus ends at 400 ms
                reinforcement = 800; 
                % Initialize kernel with zeros
                kernel = zeros(1, num_bins);
                % Define weights for the stimulus window
                stimulus_mask = (time_vector >= stimulus_onset) & (time_vector <= stimulus_offset-50);
                kernel(stimulus_mask) = 1 / sum(stimulus_mask); % Normalize weights within stimulus window
                % Convolve to calculate firing rate (specific to stimulus)
               stim_kernel(t,:) = conv2(fr_single(t,:), flip(kernel), 'same');

               if ismember (t, [hit_easy_idx; miss_easy_idx])
                   easy_go_stim_kernel(t,:) = smoothdata(circshift (stim_kernel(t,:),-(dist(baseline_onset,stimulus_onset-50)))) ;
                   easy_ngo_stim_kernel(t,:) = zeros(1, num_bins);
                   hard_go_stim_kernel(t,:) = zeros(1, num_bins);
                   hard_ngo_stim_kernel(t,:) = zeros(1, num_bins);
               elseif ismember (t, [fa_easy_idx; cr_easy_idx])
                   easy_ngo_stim_kernel(t,:) =  smoothdata(circshift (stim_kernel(t,:),-(dist(baseline_onset,stimulus_onset-50)))) ;
                   easy_go_stim_kernel(t,:) = zeros(1, num_bins);
                   hard_go_stim_kernel(t,:) = zeros(1, num_bins);
                   hard_ngo_stim_kernel(t,:) = zeros(1, num_bins);
               elseif ismember (t, [hit_hard_idx; miss_hard_idx])
                   hard_go_stim_kernel(t,:) =  smoothdata(circshift (stim_kernel(t,:),-(dist(baseline_onset,stimulus_onset-50)))) ;
                   easy_go_stim_kernel(t,:) = zeros(1, num_bins);
                   easy_ngo_stim_kernel(t,:) = zeros(1, num_bins);
                   hard_ngo_stim_kernel(t,:) = zeros(1, num_bins);
               elseif ismember (t, [fa_hard_idx; cr_hard_idx])
                   hard_ngo_stim_kernel(t,:) =  smoothdata(circshift (stim_kernel(t,:),-(dist(baseline_onset,stimulus_onset-50)))) ;
                   easy_go_stim_kernel(t,:) = zeros(1, num_bins);
                   easy_ngo_stim_kernel(t,:) = zeros(1, num_bins);
                   hard_go_stim_kernel(t,:) = zeros(1, num_bins);
               end

               %lick_kernel
               kernel = zeros(1, num_bins) ;
               if ~isempty(licks_aligned{t})
                   idx_lick = find(licks_aligned{t} + 0.2 <= (length(fr_single(t,:))/1000) ) ;
                   if~isempty(idx_lick)
                       for l = 1:length(idx_lick)
                           lick_onset = ((licks_aligned{t}(1,l) +0.2 )*1000)- 25 ;  % Stimulus starts at 200 ms
                           lick_offset = ((licks_aligned{t}(1,l) +0.2 )*1000) + 25 ;%licks_aligned{t}(end) ; % Stimulus ends at 400 ms
                           if lick_onset > length(fr_single(t,:))
                               lick_offset = length(fr_single(t,:)) - 50 ;
                           end
                           if lick_offset > length(fr_single(t,:))
                               lick_offset = length(fr_single(t,:))  ;
                           end
                           % Initialize kernel with zeros
                           kernel = zeros(1, num_bins);
                           % Define weights for the stimulus window
                           lick_mask = (time_vector >= lick_onset) & (time_vector <= lick_offset);
                           kernel(lick_mask) = 1 / sum(lick_mask); % Normalize weights within stimulus window
                           % Convolve to calculate firing rate (specific to stimulus)
                           lick_kernel_all(l,:) = conv(fr_single(t,:), flip(kernel), 'same');
                           lick_kernel_all(l,:) = smoothdata(circshift (lick_kernel_all(l,:),+ stimulus_onset/2)) ;
                       end
                       lick_kernel(t,:) = zeros(1, num_bins) ;
                       lick_kernel(t,stimulus_offset:length(fr_single(t,:))) = flip(nanmean(lick_kernel_all(1:idx_lick,stimulus_offset:length(fr_single(t,:))),1)) ;
                   end
               elseif isempty(licks_aligned{t})
                   lick_kernel(t,:) = zeros(1, num_bins) ;
               end
                      
                    bias_Go_kernel(t,:) = zeros(1, num_bins);
                    bias_NoGo_kernel(t,:) = zeros(1, num_bins);

                    if t > 1
                    % response bias Go
                    if ismember (t-1, [fa_easy_idx; fa_hard_idx; ])
                        bias_Go_kernel(t,200:end) =   1 ;
                    elseif ismember (t-1, [cr_easy_idx; cr_hard_idx])
                        bias_Go_kernel(t,200:end) =  -1 ;
                    end

                    % response bias No Go
                    if ismember (t-1, [hit_easy_idx; hit_hard_idx])
                        bias_NoGo_kernel(t,200:end) =  1 ;
                    elseif ismember (t-1, [miss_easy_idx; miss_hard_idx ])
                        bias_NoGo_kernel(t,200:end) =  -1 ;
                    end

                    end 
            end
            
            % % Create predictors for the GLM
            X_neuron =[reshape(easy_go_stim_kernel',[],1),...
                reshape(easy_ngo_stim_kernel',[],1),...
                reshape(hard_go_stim_kernel',[],1),...
                reshape(hard_ngo_stim_kernel',[],1),...
                reshape(lick_kernel',[],1),...
                reshape(bias_Go_kernel',[],1),...
                reshape(bias_NoGo_kernel',[],1)];

            [b_full, ~, stats_full] = glmfit(X_neuron, reshape(fr_single',[],1), 'normal');
            y_pred_full = glmval(b_full, X_neuron, 'identity');

            num_predictors = size(X_neuron,2); % Number of predictors

            % Firing rate per neuron
            y_FR = nanmean(fr_single(:,:),1);

            y_pred_all = reshape(y_pred_full, length(fr_single), n_trials)' ;

            y_pred_all_cell{g,i}{c,1} = y_pred_all ; 

            % Explained variance for full model
            SS_total = sum((fr_single - y_FR).^2);
            SS_residual_full = sum((fr_single - y_pred_all).^2);
            explained_var_full{g,i}(c) = 1 - (SS_residual_full / SS_total);

            % Step 2: Relative contribution of predictors (GLM coefficients normalized)
            contributions_per_predictor{g,i}(:,c) = abs(b_full(2:end)) ./ sum(abs(b_full(2:end)));

            y_pred_partial = [] ; 
            % % Step 3: Explained variance per predictor
            for p = 1:num_predictors-1
                X_partial = X_neuron(:,p:num_predictors);
                [b_partial, ~] = glmfit(X_partial, reshape(fr_single',[],1), 'normal');

                y_pred_partial(:,p) = glmval(b_partial, X_partial, 'identity') ;

                y_pred_partial_r = reshape(y_pred_partial(:,p), length(fr_single), n_trials)' ;

                y_pred_partial_r_cell{g,i}{c,p} = y_pred_partial_r ;

                SS_residual_partial = sum(( fr_single - y_pred_partial_r).^2) ;
                explained_var_per_predictor{g,i}(c,p) = 1 - (SS_residual_partial / SS_total);
            end

            y_pred_single = [] ;
            % % Step 4: Explained variance per predictor
            for p = 1:num_predictors
                X_single = X_neuron(:,p);
                [b_single, ~] = glmfit(X_single, reshape(fr_single',[],1), 'normal');

                y_pred_single(:,p) = glmval(b_single, X_single, 'identity') ;

                y_pred_single_r = reshape(y_pred_single(:,p), length(fr_single), n_trials)' ;

                y_pred_single_r_cell{g,i}{c,p} = y_pred_single_r ;

                SS_residual_single = sum(( fr_single - y_pred_single_r).^2) ;
                explained_var_single_predictor{g,i}(c,p) = 1 - (SS_residual_single / SS_total);
            end


             % Define time lags in seconds (-0.6s to 0.6s)
            bin_size = 0.001; % Example bin size in seconds
            lag_range = -0.6:0.1:0.6; % Lag range in seconds
            num_lags = length(lag_range);

             % Total bins in your time series
            N = length(fr_single(:,:)); % Assuming firing rate has N time bins

            % % % Step 5: Cross-validated explained variance for different time lags
             fr_single_2 = fr_single ;

            for lag_idx = 1:num_lags
                time_lag = lag_range(lag_idx);
                shift_bins = round(time_lag / bin_size); % Convert time lag to bin shifts
                %shift = randperm(abs(shift_bins)) ;
               % for t = 1:n_trials
                   % fr_single_2(t,1:abs(shift_bins)) = fr_single_2(t,shift) ;
                    y_FR_shifted = circshift(reshape(fr_single',[],1), shift_bins) ;
               % end
                  
                % Fit GLM for circularly shifted firing rates
                [b_time, ~] = glmfit(X_neuron, reshape(y_FR_shifted',[],1), 'normal');
                y_pred_time = glmval(b_time, X_neuron, 'identity');
                % Compute explained variance
                  y_FR_shifted =  reshape(y_FR_shifted, length(fr_single), n_trials)' ;
                  SS_total = sum((y_FR_shifted - nanmean(y_FR_shifted(:,:),1)).^2);
                 y_pred_time = reshape(y_pred_time, length(fr_single), n_trials)' ;

                y_pred_time_cell{g,i}{c,lag_idx} = y_pred_time ;

                SS_residual_time = sum((y_FR_shifted - y_pred_time).^2);
                cv_explained_time{g,i}(c, lag_idx) = 1 - (SS_residual_time / SS_total);
            end
  
        end
    end
end
    toc
   
    %% plot example neurons
    close all
    figure
    g = 1
    i = 1
    for c = 1:30
        figure
        subplot (2,1,1)
        plot(rasterX{g,i}{c,1},rasterY{g,i}{c,1},'k')
        xticks([0 200 400 600 800])
        xticklabels([-200 0 200 400 600])
        xlabel('time(ms)')
        ylabel('# trials')
        box off
        subplot(2,1,2)
        plot(mean(frArray_Smoothed{g,i}{c,1}),'k')
        hold on
        plot(mean(y_pred_all_cell{g,i}{c,1}),'r')
        hold on

        box off
        xlabel('time(ms)')
        xticks([0 200 400 600 800])
        xticklabels([-200 0 200 400 600])
        xline(200,'--k')
        xline(300,'--k')
        ylabel('FR(Hz)')
        figure
        bar (contributions_per_predictor{g,i}(:,c))
        ylim([0 0.25])
    end

    %% plot averag explained varianc,e lag and relative contribution 
    close all
    explained_var_all = nan(size(et,1),250) ;

    for g   = 1:size(et,1)
        % run per group
        p_string = {'easy go','easy no go', 'hard go', 'hard no go', 'lick','Go  bias','NoGo bias'} ;
        p_string_without = {'- easy go',' - easy no go', '- hard go', '- hard no go', '- lick','- Go  bias'} ;


        frArray_Smoothed{g, i};

        explained_var_per_group = [] ;
        explained_var_per_group = abs(cat(2,explained_var_full{g,:})') ;
        explained_var_per_group_all{g} = explained_var_per_group ;

        contributions_per_predictor_per_group = [] ;
        contributions_per_predictor_per_group = horzcat(contributions_per_predictor{g,:});
        contributions_per_predictor_all{g} = horzcat(contributions_per_predictor{g,:});

        explained_var_per_predictor_per_group = [];
        explained_var_per_predictor_per_group = abs(vertcat(explained_var_per_predictor{g,:})) ;

        explained_var_single_predictor_per_group = [];
        explained_var_single_predictor_per_group = abs(vertcat(explained_var_single_predictor{g,:})) ;
        explained_var_single_predictor_all{g} = abs(vertcat(explained_var_single_predictor{g,:})) ;

        cv_explained_time_per_group = [] ;
        cv_explained_time_per_group = vertcat(cv_explained_time{g,:}) ;

        mean_contributions = [] ;
        mean_explained_var = [] ;
        ste_contributions = [] ;
        ste_explained_var = [] ;
        mean_cv_explained = [] ;
        ste_cv_explained = [] ;

        % Population-level analysis
        mean_contributions = mean(contributions_per_predictor_per_group, 2);
        for p = 1:num_predictors
            ste_contributions(p,1) = std(contributions_per_predictor_per_group(p,:))/sqrt (size(contributions_per_predictor_per_group,2));
        end

        mean_explained_var = mean(explained_var_per_predictor_per_group, 1);
        for p = 1:num_predictors-1
            ste_explained_var(p,1) = std(explained_var_per_predictor_per_group(p,:))/sqrt (size(explained_var_per_predictor_per_group,2));
        end


        mean_cv_explained = mean(cv_explained_time_per_group, 1);
        for l = 1:length(lag_range)
            ste_cv_explained(l,1) = std(cv_explained_time_per_group(:,l))/sqrt (size(cv_explained_time_per_group,2));
        end

        % 4. performance with time lag
        % Plot mean explained variance with error bars
        figure(1);
        subplot(1,2,g)
        errorbar(lag_range, mean_cv_explained,ste_cv_explained,'linestyle', '-','Color',[.3 .3 .3]);
        hold on
        box off
        xlabel('Time Lag (s)');
        ylabel('Explained Variance');
        title('Explained Variance vs Time Lag');


        % 5. Scatter: Explained variance per neuron for full model
        figure(2);
        boxplot(explained_var_per_group,'Positions',g,'Colors','k','Symbol','.')
        hold on
        scatter(ones(size(explained_var_per_group)).*(g+(rand(size(explained_var_per_group))-0.5)/10),explained_var_per_group,'Marker','.','MarkerEdgeColor',...
            [.3 .3 .3],'MarkerFaceColor','none')
        title('Full Model Explained Variance per Neuron');
        xlabel('Neuron Index'); ylabel('Explained Variance');
        xticks([1,2])
        xticklabels({'adolescent','adult'})
        box off

        figure
        imagesc (contributions_per_predictor_per_group(:,:)')
        colormap("parula")
        colorbar;
        caxis([0 1]);
        box off
        xlabel('predictors');
        ylabel('neuron');
        title('Relative Contribution')

    end


    %% LME full explained variance 

    for g   = 1:size(et,1) % run per group

        for i =1:size(et,2) % run per recording

            for c = 1:height(frArray_Smoothed{g,i})

                var_all{g,i}(c,1) = g ;
                var_all{g,i}(c,2) = GNG_rec_all_cell{1,g}(i).Mouse ;
                var_all{g,i}(c,3) = GNG_rec_all_cell{1,g}(i).Rec ;
            end
        end
            explained_var_per_group = [] ;
        explained_var_per_group = abs(cat(2,explained_var_full{g,:})') ;
        explained_var_per_group_all{g} = explained_var_per_group ; 

            var_per_group{g} = vertcat(var_all{g,:});

    end 

    var_all = [var_per_group{1}', var_per_group{2}']';
    explained_all = [explained_var_per_group_all{1}', explained_var_per_group_all{2}']';
T = table ; 
T.Group = var_all(:,1) ;
T.Mouse = var_all(:,2) ;
T.Recording = var_all(:,3) ;
T.Condition = explained_all;

T.Group = categorical(T.Group);
T.Mouse = categorical(T.Mouse);
T.Recording = categorical(T.Recording);



% Fit the linear mixed-effects model
lme = fitlme(T, 'Condition ~ Group + (1|Mouse) + (1|Recording)')
% Display model summary
disp(lme)

% --- Extract fixed effects ---
fixedTbl = anova(lme);  % Includes F, p, etc.



%% relative contribution - overall
 

clear contributions_per_predictor_all
for g = 1:2

    contributions_per_predictor_per_group = [] ;
        contributions_per_predictor_per_group = horzcat(contributions_per_predictor{g,:});

        contributions_per_predictor_all{g} = contributions_per_predictor_per_group ; 

 % Population-level analysis
        mean_contributions_all{g} = mean(contributions_per_predictor_per_group, 2);
        for p = 1:num_predictors
            ste_explained_var_cell{g}(p,1) = std(contributions_per_predictor_per_group(p,:))/sqrt (size(contributions_per_predictor_per_group,2));
        end
end 

figure 
hold on 
          % go easy
        errorbar(1,mean_contributions_all{1}(1,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(1,mean_contributions_all{1}(1,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(2,mean_contributions_all{2}(1,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(2,mean_contributions_all{2}(1,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');
        % no go easy
        errorbar(3,mean_contributions_all{1}(2,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(3,mean_contributions_all{1}(2,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(4,mean_contributions_all{2}(2,:), ste_contributions_all{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(4,mean_contributions_all{2}(2,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');

                 % go hard
        errorbar(5,mean_contributions_all{1}(3,:), ste_contributions_all{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(5,mean_contributions_all{1}(3,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(6,mean_contributions_all{2}(3,:), ste_contributions_all{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(6,mean_contributions_all{2}(3,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');
        % no go hard 
        errorbar(7,mean_contributions_all{1}(4,:), ste_contributions_all{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(7,mean_contributions_all{1}(4,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(8,mean_contributions_all{2}(4,:), ste_contributions_all{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(8,mean_contributions_all{2}(4,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');


        % LIck
        errorbar(9,mean_contributions_all{1}(5,:), ste_contributions_all{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(9,mean_contributions_all{1}(5,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(10,mean_contributions_all{2}(5,:), ste_contributions_all{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(10,mean_contributions_all{2}(5,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');

        % go bias
        errorbar(11,mean_contributions_all{1}(6,:), ste_contributions_all{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(11,mean_contributions_all{1}(6,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(12,mean_contributions_all{2}(6,:), ste_contributions_all{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(12,mean_contributions_all{2}(6,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');

        % no go Bias
        errorbar(13,mean_contributions_all{1}(7,:), ste_contributions_all{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(13,mean_contributions_all{1}(7,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(14,mean_contributions_all{2}(7,:), ste_contributions_all{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(14,mean_contributions_all{2}(7,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');

    
        box off
        title('Average Relative Contribution of Predictors');
        xlabel('Predictors'); ylabel('Mean Contribution'); 
        ylim([0 0.25])

    for n = 1:num_predictors
        [p_contrib(n),~,stats] = ranksum( contributions_per_predictor_all{1}(n,:),contributions_per_predictor_all{2}(n,:));
        z_contrib(n) = stats.zval ; 
            [~,p_contrib,stats]= bonferroni_holm( p_contrib,0.01) ;

    end

  

        %% explained variance per predictor

        clear explained_var_per_predictor_all
        for g = 1:2
            contributions_per_predictor_per_group = [] ;
            explained_var_per_predictor_per_group = abs(vertcat(explained_var_per_predictor{g,:})) ;
            explained_var_per_predictor_all{g} = explained_var_per_predictor_per_group ;

            explained_var_per_predictor_all{g}(:,4) = mean(explained_var_per_predictor_all{g}(:,1:4),2) ;
            mean_explained_var_cell{g} = mean(explained_var_per_predictor_per_group, 1)' ;
            for p = 1:num_predictors-1
                ste_explained_var_cell{g}(p,1) = std(explained_var_per_predictor_per_group(p,:))/sqrt (size(explained_var_per_predictor_per_group,1));
            end
        end


     figure
hold on 
          % go easy
        errorbar(1,mean_explained_var_cell{1}(1,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(1,mean_explained_var_cell{1}(1,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(2,mean_explained_var_cell{2}(1,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(2,mean_explained_var_cell{2}(1,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');
        % no go easy
        errorbar(3,mean_explained_var_cell{1}(2,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(3,mean_explained_var_cell{1}(2,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(4,mean_explained_var_cell{2}(2,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(4,mean_explained_var_cell{2}(2,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');

                 % go hard
        errorbar(5,mean_explained_var_cell{1}(3,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(5,mean_explained_var_cell{1}(3,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(6,mean_explained_var_cell{2}(3,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(6,mean_explained_var_cell{2}(3,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');
        % no go hard 
        errorbar(7,mean_explained_var_cell{1}(4,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(7,mean_explained_var_cell{1}(4,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(8,mean_explained_var_cell{2}(4,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(8,mean_explained_var_cell{2}(4,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');


        % LIck
        errorbar(9,mean_explained_var_cell{1}(5,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(9,mean_explained_var_cell{1}(5,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(10,mean_explained_var_cell{2}(5,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(10,mean_explained_var_cell{2}(5,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');

        % go bias
        errorbar(11,mean_explained_var_cell{1}(6,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(11,mean_explained_var_cell{1}(6,:),'FaceColor',[.5 .5 .5],'EdgeColor','none');
        errorbar(12,mean_explained_var_cell{2}(6,:), ste_explained_var_cell{1}(5,:),'linestyle','none','Color',[0 0 0])
            bar(12,mean_explained_var_cell{2}(6,:),'FaceColor',[.2 .2 .2],'EdgeColor','none');

 
        box off
        title('Explained Variance per predictor');
        xlabel('Predictors'); ylabel('Explained Variance'); 

    for n = 1:num_predictors-1
        [P_var(n),~,stats] = ranksum( explained_var_per_predictor_all{1}(:,n),explained_var_per_predictor_all{2}(:,n));
        z_var(n) = stats.zval ; 
    
            [~,P_var,stats]= bonferroni_holm( P_var,0.05) ;

    end
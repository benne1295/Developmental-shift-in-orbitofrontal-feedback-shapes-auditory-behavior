clearvars -except GNG_rec_all_cell

startRange = -0.2
stopRange = 0.6
binSize = 0.001
smoothSize = 20
min_spikes = 100 ;
rasterScale = 5
colors_ado_adu = {[.5 .7 .2],[ .2 .4 .2]} ;

%% exrtact event_times and spike_times of all units
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

% test if neurons are modulated by the tone
for g   = 1:size(et,1) % run per group

    for i = 1:size(et,2) % run per recording

        spikeTimes = st{g,i}  ; %extract spike times

        eventTimes_all = et{g, i}.all ;

        for c = 1:length(st{g,i})  %run per cluster = per neuron

            %sort out spike times outside of the event time range
            spikeTime = spikeTimes{c,1}(spikeTimes{c,1}>min(eventTimes_all+startRange) & spikeTimes{c,1}<max(eventTimes_all+stopRange));
            if length(spikeTime) >=  min_spikes

                %extract the basics
                [~, psth,GNG_analysis.binArray{g,i}{c,1},binCenters] = GNG_binning_PSTH(startRange,stopRange, binSize, eventTimes_all, spikeTime);

                % [rasterX{g,i}{c,1}, rasterY{g,i}{c,1}] = GNG_raster (GNG_analysis.binArray{g,i}{c,1}, binSize, binCenters, rasterScale);

                % extract the smoothed basics
                [~,psth_Smoothed, GNG_analysis.frArray_Smoothed{g,i}{c,1},spikeCounts] ...
                    = GNG_smoothed_PSTH (startRange, stopRange, binCenters,GNG_analysis.binArray{g,i}{c,1}, binSize, smoothSize);

                baseline_psth = psth (1,1:150) ;
                onset_psth = psth (1,101:350) ;

                %excitated neurons
                % compare the psth of tone and baseline
                [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','right');
                if  p < 0.05 % sort out non-excited units
                    clusterIDs_exc{g,i}(c,:) = c ;
                end

                %supressed neurons
                [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','left');
                if  p < 0.05 % sort out non-excited units
                    clusterIDs_inh{g,i}(c,:) = c ;

                end

                %look at both excited and inhibited neurons or population decoding
                [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','both');
                if  p < 0.05 % sort out non-excited units
                    clusterIDs_modh{g,i}(c,:) = c ;

                end

            end
        end

        if ~isempty(clusterIDs_modh{g,i})
            clusterIDs_modh{g,i}(clusterIDs_modh{g,i}==0) = [] ;
        end

        if ~isempty(clusterIDs_inh{g,i})
            clusterIDs_inh{g,i}(clusterIDs_inh{g,i}==0) = [] ;
        end

        if ~isempty(clusterIDs_exc{g,i})
            clusterIDs_exc{g,i}(clusterIDs_exc{g,i}==0) = [] ;
        end


        for cn =  1:length(clusterIDs_modh{g,i})
            st_mod{g,i}{cn,1} =  st{g,i}{clusterIDs_modh{g,i}(cn,1),1}
        end

    end
end
%% create FR_array_matrix
clc
tic
% test if neurons are modulated by the tone
for g   = 1:size(et,1) % run per group

    for i = 1:size(et,2) % run per recording
        spikeTimes = st_mod{g,i}  ; %extract spike times

        et_single = et{g, i};
        % Extract indices for different trial types
                hit_easy_idx = et_single.tr_resp_idx(1, ~isnan(et_single.tr_resp_idx(1, :)))';
        fa_easy_idx = et_single.tr_resp_idx(2, ~isnan(et_single.tr_resp_idx(2, :)))';
        cr_easy_idx = et_single.tr_resp_idx(4, ~isnan(et_single.tr_resp_idx(4, :)))';

        fa_hard_idx = et_single.tr_resp_idx(8, ~isnan(et_single.tr_resp_idx(8, :)))';
        cr_hard_idx = et_single.tr_resp_idx(10, ~isnan(et_single.tr_resp_idx(10, :)))';

        %------------- Lick bias -----------
        %after licks
        licks = sort([  fa_easy_idx' fa_hard_idx']) ;
        after_licks = licks + 1 ;
        after_licks = after_licks(1:end-2) ;
        % for specific stim
        after_licks = intersect(after_licks,hit_easy_idx)' ;
        events_after_licks = et{g, i}.all(after_licks) ;
        %after no licks
        nlicks = sort([cr_easy_idx'  cr_hard_idx']) ;
        after_nlicks = nlicks + 1 ;
        after_nlicks = after_nlicks(1:end-2) ;
        % for specific stim
        after_nlicks = intersect(after_nlicks,hit_easy_idx)' ;
        events_after_nlicks = et{g, i}.all(after_nlicks) ;


        for c = 1:length(st_mod{g,i})  %run per cluster = per neuron

            %------------- Lick bias -----------

            %sort out spike times outside of the event time range
            spikeTimes_after_licks = spikeTimes{c,1}(spikeTimes{c,1} > min(events_after_licks + startRange) & spikeTimes{c,1} < max(events_after_licks + stopRange));
            spikeTimes_after_nlicks = spikeTimes{c,1}(spikeTimes{c,1} > min(events_after_nlicks + startRange) & spikeTimes{c,1} < max(events_after_nlicks + stopRange));

            if length(spikeTimes_after_licks) >=  min_spikes
                if length(spikeTimes_after_nlicks) >=  min_spikes
                    % Compute binned and smoothed firing rates
                    [~, ~, binArray_alicks{g, i}{c,1}, binCenters] = GNG_binning_PSTH(startRange, stopRange, binSize, events_after_licks, spikeTimes_after_licks);
                    [~, ~, binArray_anlicks{g, i}{c,1}, ~] = GNG_binning_PSTH(startRange, stopRange, binSize, events_after_nlicks, spikeTimes_after_nlicks);

                    [~, psth_Smoothed_alicks{g,i}(c,:), frArray_Smoothed_alicks{g, i}{c,1}, ~] = GNG_smoothed_PSTH(startRange, stopRange, binCenters, binArray_alicks{g, i}{c,1}, binSize, smoothSize);
                    [~,psth_Smoothed_anlicks{g,i}(c,:), frArray_Smoothed_anlicks{g, i}{c,1}, ~] = GNG_smoothed_PSTH(startRange, stopRange, binCenters, binArray_anlicks{g, i}{c,1}, binSize, smoothSize);
                end
            end


        end
    end
end

% take out the cells that didn't pass the thresholds of minimal spikes
for g   = 1:size(et,1) % run per group
    for i = 1:size(et,2) % run per recording
        if ~isempty( psth_Smoothed_alicks{g,i})
            psth_Smoothed_alicks{g,i}(sum(psth_Smoothed_alicks{g,i}, 2) == 0, :) = [];

        end
        if ~isempty( psth_Smoothed_anlicks{g,i})
            psth_Smoothed_anlicks{g,i}(sum(psth_Smoothed_anlicks{g,i}, 2) == 0, :) = [];
        end


    end
end

%% Perform PCA on the full 800 time points for each trial
clc
close all

for g = 1:size(et, 1)
    distances = nan (7,800) ;
    for i = 1:size(et, 2)

        et_single = et{g, i};
        g
        i

        events_A = psth_Smoothed_alicks{g, i};
        events_B = psth_Smoothed_anlicks{g, i};


        if ~isempty( events_A)
            if ~isempty( events_B)

                if height(events_A) > 10
                    if height(events_B) > 10

                        % Combine data for all trials and create labels
                        data = [events_A; events_B]; % Rows = trials, Columns = 800 time bins
                        data_ = [events_A; events_B]; % Rows = trials, Columns = 800 time bins

                        labels = [ones(size(events_A, 1), 1); ones(size(events_B, 1), 1) * 2]; % 1 = hit, 2 = CR

                        % Standardize the data (z-score normalization per time bin)
                        data_standardized_A = zscore(events_A, 0, 1); % Mean = 0, Std Dev = 1 for each time bin (feature)
                        data_standardized_B = zscore(events_B, 0, 1); % Mean = 0, Std Dev = 1 for each time bin (feature)

                        % Perform PCA on the full dataset
                        [coeff_A, score_A, ~, ~, explained_A] = pca(data_standardized_A);

                        [coeff_B, score_B, ~, ~, explained_B] = pca(data_standardized_B);

                        % Calculate distances
                        mean_dist_A(g,i) = nanmean(vecnorm(diff(score_A(:, 1:3), 1, 1), 2, 2));
                        cum_dist_A(g,i) = sum(vecnorm(diff(score_B(:, 1:3), 1, 1), 2, 2));

                        mean_dist_B(g,i) = nanmean(vecnorm(diff(score_B(:, 1:3), 1, 1), 2, 2));
                        cum_dist_B(g,i) = sum(vecnorm(diff(score_B(:, 1:3), 1, 1), 2, 2));

                        % Plot clusters in the first three PCs
                        % 
                        % figure;
                        % %subplot(3,3,i)
                        % initialColor = [0.5, 0, 0.5];  
                        % for j = 1:height(coeff_A)
                        %     colorShift = 0.0005;
                        %     currentColor = initialColor + colorShift * j;
                        %     currentColor = min(max(currentColor, 0), 1);
                        %     scatter3(coeff_A(j, 1), coeff_A(j, 2), coeff_A(j, 3), 10, currentColor, 'filled');
                        %     hold on
                        % end
                        % scatter3(coeff_A(1, 1), coeff_A(1, 2), coeff_A(1, 3), 10, 'k', 'Marker','*');
                        % scatter3(coeff_A(height(coeff_A), 1), coeff_A(height(coeff_A), 2), coeff_A(height(coeff_A), 3), 10, 'r', 'Marker','*');

                        % 
                        % hold on;
                        % initialColor = [0.5, 0.5, 0]	;  
                        % for j = 1:height(coeff_B)
                        %     colorShift = 0.0005;
                        %     currentColor = initialColor + colorShift * j;
                        %     currentColor = min(max(currentColor, 0), 1);
                        %     scatter3(coeff_B(j, 1), coeff_B(j, 2), coeff_B(j, 3), 10, currentColor, 'filled');
                        %     hold on
                        % end
                        % scatter3(coeff_B(1, 1), coeff_B(1, 2), coeff_B(1, 3), 10, 'k', 'Marker','*');
                        % scatter3(coeff_B(height(coeff_B), 1), coeff_B(height(coeff_B), 2), coeff_B(height(coeff_B), 3), 10, 'r', 'Marker','*');
                        % 
                        % xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
                        % hold off;


                        coeff_A_reduced = nanmean(coeff_A(:, 1:3),2); % Using the first 2 components of USV coefficients
                        coeff_B_reduced = nanmean(coeff_B(:, 1:3),2); % Using the first 2 components of NBN coefficients

            
                            % Calculate the distance between the coefficients of the corresponding components
                            distances(i,:) = sqrt(sum((coeff_A_reduced - coeff_B_reduced).^2, 2))'; % Euclidean distance between coefficients
         
                    end
                end
            end
        end

    end
    distances_cell{1,g} = distances;



    % Plot SEM as a patch
    x = [1:800, fliplr(1:800)];


    cum_distance(1,g,:) = cumsum(nanmean(distances(:,:))) ;
    sem_cum_distance(1,g,:) = cumsum(nanstd(distances(:,:)))/sqrt(7) ;

    figure(1+4)
    plot(1:800, squeeze(cum_distance(1,g,:)), 'Color', colors_ado_adu{g}, 'LineWidth', 3)
    hold on
    % Plot SEM as a patch
    y_cum = [squeeze(cum_distance(1,g,:))' + squeeze(sem_cum_distance(1,g,:))', ...
        fliplr(squeeze(cum_distance(1,g,:))' - squeeze(sem_cum_distance(1,g,:))')];
    patch('XData', x, 'YData', y_cum, 'FaceColor', colors_ado_adu{g}, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
    xticks([0 200 400 600 800])
    xticklabels([-200 0 200 400 600])
    yticks([0 5 10 15 20])
    ylim([0  20])

    xline(200,'--k')
    xline(300,'--k')
    box off
    xlabel('time(ms)')
    ylabel('cumulative distance')




end


h = kstest( distances_cell{1,1}(:,:) )
for bias = 1
    for bin = 1:size(distances_cell{1,1}(:,:),2)
        [p(bin,:),~,~] = ranksum(distances_cell{bias,1}(:,bin),...
            distances_cell{bias,2}(:,bin),'alpha',0.05,'tail','both') ;
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


end

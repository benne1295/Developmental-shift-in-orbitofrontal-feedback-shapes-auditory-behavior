clc 
close all
clearvars -except GNG_rec_all_cell

%%
addpath(genpath('Z:\Shared\Benne\Praegel_et_al_2025\analysis'))

[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2025\data\optogenetics'...
    , ['Select GNG_rec_all_cell ']);
addpath(path)
load (file)

%% Parameters
clc
tone_dur = 0.1 ; % sec
response_window = 2 ; % sec
stim = [1:4] ; % ID
freqs = [7.07 9.17 10.95 14.14] ; % frequencies of stim IDs

color_eh= {[0.3010 0.7450 0.9330], [0 0.4470 0.7410]} ;
color_opto = {[0.6350 0.0780 0.1840]} ;
color_opto_fade = {[0.6350 0.0780 0.1840 0.5]} ;
color_mean_marker = {[0.5 0.5 0.5 0.5]} ;

color_plot_bias = {[0.4660 0.6740 0.1880 0.5] [0.8500 0.3250 0.0980 0.5],[0.9290 0.6940 0.1250 0.5], [0.4940 0.1840 0.5560 0.5]};
color_plot_bias_opto = {[0.4660 0.6740 0.3880 0.5], [0.8500 0.3250 0.0680 0.5],[0.9290 0.6940 0.3250 0.5], [0.4940 0.1840 0.8560 0.5]};
color_plot = {[0.4660 0.6740 0.1880] [0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]};

L = {['-'], ['-']} ;
L2 = {['--'], ['-']} ;

tone_dur = 0.1 ; % sec
response_window = 2 ; % sec
dprime_threshold = 1 ;
min_binsize = 30;
startRange = 0; %beginning of baseline
stopRange = 0.1; %end of the offset
startReinf = 0.6 ; % start reinforcement
length_base = 200 ; %ms
reinforcement = 500 + length_base ; %ms
max_trial = 1000 ; % n trials
trial_cut_off = 100 ; % n tirals
startRange_resp = - 0.1 ;
stopRange_resp = 2 ;
binSize = 0.001 ;
smoothSize = 5 ;
length_stim = (abs(startRange_resp) + stopRange_resp)*1000 ;
trial_criterion = 5 ;

stim = [1 2 3 4 5 6 7];

binBorders = startRange_resp:binSize:stopRange_resp ; %define the time points
numBins = length(binBorders)-1 ;

colors_ado_adu = {[.5 .7 .2],[ .2 .4 .2]} ;
directory = 'Z:\Shared\Benne\Praegel_et_al_2025\figures'

%% caclulate discrimination and lick performance 
close all
clc

for g = 1:size(GNG_rec_all_cell,2)
    for i   = 1:size(GNG_rec_all_cell{1,g},1) % mouse

        stim_types = GNG_rec_all_cell{1,g}(i).Behavior.stim_types;
        stim_ids =GNG_rec_all_cell{1,g}(i).Behavior.stim_ids;
        trial_responses = GNG_rec_all_cell{1,g}(i).Behavior.trial_responses;
        lick_times = GNG_rec_all_cell{1,g}(i).Behavior.lick_times;
        stim_times = GNG_rec_all_cell{1,g}(i).Behavior.stim_times;

          % retrieve trial identities
        hit_ID{g,i} = find (stim_types == 1 & trial_responses == 1);
        FA_ID{g,i} = find (stim_types == -1 & trial_responses == 1);
        miss_ID{g,i} = find (stim_types == 1 & trial_responses == 0);
        CR_ID{g,i} = find (stim_types == -1 & trial_responses == 0);

        % easy
        index_stim = find (stim_ids == stim(1) | stim_ids == stim(4));
        [go_licks_all_easy(g,i),ngo_licks_all_easy(g,i),dprimes_easy(g,i),~,~,~,cbias_easy(g,i)] =....
            GNG_lick_rt_dprime ...
            (index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

        % hard
        index_stim = find (stim_ids == stim(2) | stim_ids == stim(3));
        [go_licks_all_hard(g,i),ngo_licks_all_hard(g,i),dprimes_hard(g,i),~,~,~,cbias_hard(g,i)] =...
            GNG_lick_rt_dprime ...
            (index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

        %mean
        index_stim = find (stim_ids == stim(2) | stim_ids == stim(3)| stim_ids == stim(1) | stim_ids == stim(4));
        [go_licks_all(g,i),ngo_licks_all(g,i),dprimes(g,i),~,~,~,cbias(g,i)] =...
            GNG_lick_rt_dprime ...
            (index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

        mouse_ID(g,i) = GNG_rec_all_cell{1,g}(i).Mouse ;
        session_ID(g,i) = GNG_rec_all_cell{1,g}(i).Rec ;

        % lick raster of recording session
        eventTimes_all = GNG_rec_all_cell{1, g}(i).eventTimes.all;
        eventTimes_licks = GNG_rec_all_cell{1, g}(i).eventTimes.times_vec_licks ;

        % creater raster matrix of all lick times across all trials
        for it = 1:length(eventTimes_all) ;
            licks_after_stim = lick_times (lick_times > stim_times(it)) ;
            licks_per_trial = licks_after_stim(licks_after_stim < stim_times(it)+tone_dur+response_window) ;
            licks_from_stim_onset = licks_per_trial - stim_times(it) ;
            licks_raster{g,i}(it,1:length(licks_from_stim_onset)) = licks_from_stim_onset ;
        end

        eventTimes_all = GNG_rec_all_cell{1, g}(i).eventTimes.all;
        eventTimes_licks = GNG_rec_all_cell{1, g}(i).eventTimes.times_vec_licks ;

        lick_time{g,i} = cell (length(eventTimes_all),1) ;
        for E = 1:length(eventTimes_all)

            idx_lick{E,1} =  find(eventTimes_licks > (eventTimes_all(1,E) + startRange)...
                & eventTimes_licks < (eventTimes_all(1,E) + stopRange)) ;
            if ~isempty(idx_lick{E,1})
                lick_time{g,i}{E,1} =  eventTimes_licks(idx_lick{E,1}) - eventTimes_all(1,E) ;
            end
        end

 
    end
end
%% Lick Raster - still make the CRs purple 

close all
tic
% choose example mouse
g = 2 % adult
i = 5 % recording

size_E =  size(GNG_rec_all_cell{1,g}(i).eventTimes.all,2) ;
licks_raster{g,i} (licks_raster{g,i} == 0) = nan ;
figure
for it = 1:size_E
    if ismember (it, hit_ID{g,i})
        rectangle( 'Position' , [0 it tone_dur it+1],'Facecolor',[0.4660 0.6740 0.1880],'EdgeColor','none')
        hold on
        scatter(licks_raster{g,i}(it,:) , it + 0.5, 'MarkerFaceColor' ,[0.5 0.5 0.5],'Marker','.','MarkerEdgeColor',[0.5 0.5 0.5])
        hold on
        if  ~isempty(lick_time{g,i}{it,1})
            scatter(lick_time{g,i}{it,1}(:,:), it + 0.5, 'MarkerFaceColor' ,[0.5 0.5 0.5],'Marker','.','MarkerEdgeColor',[0.5 0.5 0.5])
        end
    elseif ismember (it, FA_ID{g,i})
        rectangle( 'Position' , [0 it tone_dur it+1],'Facecolor',[0.9290 0.6940 0.1250],'EdgeColor','none')
        scatter(licks_raster{g,i}(it,:) , it + 0.5 , 'MarkerFaceColor',[0.5 0.5 0.5],'Marker','.','MarkerEdgeColor',[0.5 0.5 0.5])

        if  ~isempty(lick_time{g,i}{it,1})
            scatter(lick_time{g,i}{it,1}(:,:), it + 0.5, 'MarkerFaceColor' ,[0.5 0.5 0.5],'Marker','.','MarkerEdgeColor',[0.5 0.5 0.5])
        end

    elseif ismember (it, miss_ID{g,i})
        rectangle( 'Position' , [0 it tone_dur it+1],'Facecolor',[0.8500 0.3250 0.0980],'EdgeColor','none')
    elseif ismember (it, CR_ID{g,i})
        rectangle( 'Position' , [0 it tone_dur it+1],'Facecolor',[0 0.4470 0.7410],'EdgeColor','none')
    end
end
rectangle( 'Position' , [0 (size_E+1) tone_dur (size_E+1)+1],'Facecolor',[1 1 1],'EdgeColor','none')

clear xlim
xlim ([0 tone_dur+response_window])
xticks ([0 0.1 0.6 2])
xlabel ('time(sec)')
ylabel ('trials')
xline(0.6,'-.k')

ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
movegui('east');
 
%% plot psychometric curve

L = {['-'],['--'],['-'],['--']}
go_licks_all_easy(go_licks_all_easy == 0) = nan ; 
go_licks_all_hard(go_licks_all_hard == 0) = nan ; 
ngo_licks_all_hard(ngo_licks_all_hard == 0) = nan ; 
ngo_licks_all_easy(ngo_licks_all_easy == 0) = nan ; 

clc
close all
figure
for g = 1:size(GNG_rec_all_cell,2)

    mean_licks = nanmean([go_licks_all_easy(g,:)' go_licks_all_hard(g,:)' ngo_licks_all_hard(g,:)' ngo_licks_all_easy(g,:)']) ;
    ste_licks = nanstd([go_licks_all_easy(g,:)' go_licks_all_hard(g,:)' ngo_licks_all_hard(g,:)' ngo_licks_all_easy(g,:)'],1,1)

    errorbar(freqs, mean_licks ,ste_licks,'Color',[.5 .5 .5],'linestyle',L{g},'linewidth',3)
    hold on

    ylim([0 1]);
yticks([0:0.2:1])
yline(0.5,'--','Color','k','linewidth',2);
xline(10,'--','Color','k','linewidth',2);
xlim([7 14.1])
xticks([7 9.2 10.9 14.1])
xlabel('freq (kHz)');
ylabel('lick rate');
set(gca, 'XDir','reverse','XScale','log');
box off;
hold on    
ax = gca;
ax.XAxis.FontSize = 20 ;
ax.YAxis.FontSize = 20 ;
ax.Title.FontSize = 20 ;
ax.Title.FontSize = 20 ;

end


%% dprime 

 clc
 close all
clear stats
 figure

 violinplot([dprimes(1,:) ;dprimes(2,:)]',...
     {'adolescent ','adult'},"ViolinColor",{[0.5 0.5 0.5; 0.5 0.5 0.5]}) ;
 xticks ([1 2])
 xlim([0 3]);
 ylim([0 4]);
 yticks([ 0 1 2 3 4])
 box off;
 hold on
 ylabel("discrimination (d') ")
 yline(1,'linestyle','--','linewidth',2)
 ax = gca;
 ax.XAxis.FontSize = 20;
 ax.YAxis.FontSize = 20;
 movegui('east');
 box off;
  h = kstest(dprimes_easy) ;
  [p(1,1),h,stats{1,1}] = ranksum(dprimes(1,:),dprimes(2,:),'tail','left') ;

%% cbias 
 clc
 close all
clear stats
 figure

 violinplot([cbias(1,:) ;cbias(2,:)]',...
     {'adolescent ','adult'},"ViolinColor",{[0.5 0.5 0.5; 0.5 0.5 0.5]}) ;
 xticks ([1 2])
 xlim([0 3]);
 ylim([-2 1.5]);
 yticks([ -2 -1 0 1])
 box off;
 hold on
 ylabel("c-bias")
 yline(0,'linestyle','--','linewidth',2)
 ax = gca;
 ax.XAxis.FontSize = 20;
 ax.YAxis.FontSize = 20;
 movegui('east');
 box off;
  h = kstest(cbias(1,:)) ;
  [p(1,1),h,stats{1,1}] = ranksum(cbias(1,:),cbias(2,:),'tail','left') ;

  
%% caclulate the bias per stimulus  difficulty
close all
clc
stim = [1 2 3 4 5 6 7];

for g = 1:size(GNG_rec_all_cell,2)
    for i   = 1:size(GNG_rec_all_cell{1,g},1) % mouse

        stim_types = GNG_rec_all_cell{1,g}(i).Behavior.stim_types;
        stim_ids =GNG_rec_all_cell{1,g}(i).Behavior.stim_ids;
        trial_responses = GNG_rec_all_cell{1,g}(i).Behavior.trial_responses;
        lick_times = GNG_rec_all_cell{1,g}(i).Behavior.lick_times;
        stim_times = GNG_rec_all_cell{1,g}(i).Behavior.stim_times;

        %hit
        idx_bias{1} = find(stim_types == 1 & trial_responses == 1 );
        idx_bias{1} = idx_bias{1}(1:end-1) +1 ;
        %miss
        idx_bias{2} =  find(stim_types == 1 & trial_responses == 0);
        idx_bias{2} = idx_bias{2}(1:end-1) +1 ;
        %fa
        idx_bias{3} = find (stim_types == -1 & trial_responses == 1);
        idx_bias{3} = idx_bias{3}(1:end-1) +1 ;
        %cr
        idx_bias{4} = find (stim_types == -1 & trial_responses == 0);
        idx_bias{4} = idx_bias{4}(1:end-1) +1 ;

        % all biases
        for b = 1:size(idx_bias,2)

            % easy bias easy
            index_stim = find (stim_ids == stim(1) | stim_ids == stim(4));
            if length(idx_bias{b}) > 1
                idx_intersect = intersect (index_stim,idx_bias{b}) ;
                if length(idx_intersect) > 1

                    [go_licks_all_bias_easy(g,i,b),ngo_licks_all_bias_easy(g,i,b),~,~,~,~] =....
                        GNG_lick_rt_dprime ...
                        (idx_intersect, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

                elseif length(idx_intersect) <= 1
                    go_licks_all_bias_easy(g,i,b) = nan ;
                    ngo_licks_all_bias_easy(g,i,b) = nan ;
                end
            elseif length(idx_bias{b})  <= 1
                go_licks_all_bias_easy(g,i,b) = nan ;
                ngo_licks_all_bias_easy(g,i,b) = nan ;
            end
            % hard bias
            index_stim = find (stim_ids == stim(2) | stim_ids == stim(3));
            if length(idx_bias{b}) > 1
                idx_intersect = intersect (index_stim,idx_bias{b}) ;
                if length(idx_intersect) > 1

                    [go_licks_all_bias_hard(g,i,b),ngo_licks_all_bias_hard(g,i,b),~,~,~,~] =....
                        GNG_lick_rt_dprime ...
                        (idx_intersect, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

                elseif length(idx_intersect) <= 1
                    go_licks_all_bias_hard(g,i,b) = nan ;
                    ngo_licks_all_bias_hard(g,i,b) = nan ;
                end
            elseif length(idx_bias{b})  <= 1
                go_licks_all_bias_hard(g,i,b) = nan ;
                ngo_licks_all_bias_hard(g,i,b) =nan ;
            end
        end
    end
end
%% create mean difficulty lick curve

for g = 1:size(GNG_rec_all_cell,2)

    for b = 1:size(idx_bias,2)
        for i   = 1:numel(GNG_rec_all_cell{1,g}) % mouse

            licks_mice_bias(i,:) = [go_licks_all_bias_easy(g,i,b) go_licks_all_bias_hard(g,i,b) ngo_licks_all_bias_hard(g,i,b) ngo_licks_all_bias_easy(g,i,b)] ;

        end

        mean_licks_bias = nanmean(licks_mice_bias) ;
        stderr_licks_bias = nanstd(licks_mice_bias)/sqrt(size(licks_mice_bias,1)) ;

        figure(1)
        subplot(1,2,g)
        errorbar(freqs, mean_licks_bias ,stderr_licks_bias,'Color',color_plot_bias{b},'linestyle',L2{2},'linewidth',3)
        hold on
        ylim([0 1]);
        yticks([0:0.2:1])
        yline(0.5,'--','Color','k','linewidth',2);
        xline(10,'--','Color','k','linewidth',2);
        xlim([6 16])
        xticks([7 9.2 10.9 14.1])

        xlabel('freq (kHz)');
        ylabel('lick rate');
        set(gca, 'XDir','reverse','XScale','log');
        box off;
        hold on
        ax = gca;
        ax.XAxis.FontSize = 20 ;
        ax.YAxis.FontSize = 20 ;
        ax.Title.FontSize = 20 ;
        ax.Title.FontSize = 20 ;

    end
end


%% plot the  bias of Go and No-Go responses 
close all
clear licks_mice_hit licks_mice_fa licks_mice_cr licks_mice_miss  delta_diff_all
clear stats
for g = 1:size(GNG_rec_all_cell,2)

    session_ID_g = reshape(session_ID(g,:,:),[],1) ;
    session_ID_g = session_ID_g(~isnan(session_ID_g)) ;

    for b = 1:size(idx_bias,2)
        for i   = 1:numel(GNG_rec_all_cell{1,g}) % mouse

            licks_mice_hit(i,:) =  [go_licks_all_bias_easy(g,i,1) go_licks_all_bias_hard(g,i,1) ngo_licks_all_bias_hard(g,i,1) ngo_licks_all_bias_easy(g,i,1)] ;
            licks_mice_fa(i,:) =   [go_licks_all_bias_easy(g,i,3) go_licks_all_bias_hard(g,i,3) ngo_licks_all_bias_hard(g,i,3) ngo_licks_all_bias_easy(g,i,3)] ;
            licks_mice_cr(i,:) =  [go_licks_all_bias_easy(g,i,4) go_licks_all_bias_hard(g,i,4) ngo_licks_all_bias_hard(g,i,4) ngo_licks_all_bias_easy(g,i,4)] ;
            licks_mice_miss(i,:) =  [go_licks_all_bias_easy(g,i,2) go_licks_all_bias_hard(g,i,2) ngo_licks_all_bias_hard(g,i,2) ngo_licks_all_bias_easy(g,i,2)] ;

        end
    end

    licks_mice_hit(licks_mice_hit == 0) = nan ;
    licks_mice_fa(licks_mice_fa == 0) = nan ;
    licks_mice_cr(licks_mice_cr == 0) = nan ;
    licks_mice_miss(licks_mice_miss == 0) = nan ;
    

    delta_diff_all(g,1,1:numel(GNG_rec_all_cell{1,g})) =  nanmean((licks_mice_hit) - (licks_mice_miss),2)' ;
    delta_diff_all(g,2,1:numel(GNG_rec_all_cell{1,g})) =  nanmean((licks_mice_fa) - (licks_mice_cr),2)' ;
end



figure
 violinplot([ squeeze(delta_diff_all(1,1,:))....
     squeeze(delta_diff_all(2,1,:))],...
     {'adolescent ','adult '},"ViolinColor",{[.5 .5 .5; .5 .5 .5]}) ;

 xlim([0 3])
 xticks([1 2])
 ylim([-0.5 0.75])
 yticks([-0.5:0.25:0.75])
 ylabel("delta lick rate ")
 yline(1,'linestyle','--','linewidth',2)
 xticklabels({'adolescent','adult'})
 ax = gca;
 ax.XAxis.FontSize = 20;
 ax.YAxis.FontSize = 20;
 movegui('east');
 box off;
 yline(0,'--','Color','k','linewidth',2);


 [p(1,1),h,stats{1,1}] = ranksum(squeeze(delta_diff_all(1,1,:)),squeeze(delta_diff_all(2,1,:)),'tail','right') ;
 [p(2,1),h,stats{2,1}] = signrank(squeeze(delta_diff_all(1,1,:))) ;
 [p(3,1),h,stats{3,1}] = signrank(squeeze(delta_diff_all(2,1,:))) ;
p

 figure
 violinplot([ squeeze(delta_diff_all(1,2,:))....
     squeeze(delta_diff_all(2,2,:))],...
     {'adolescent ','adult '},"ViolinColor",{[.5 .5 .5; .5 .5 .5]}) ;

 xlim([0 3])
 xticks([1 2])
 ylim([-0.5 0.75])
 yticks([-0.5:0.25:0.75])
 ylabel("delta lick rate ")
 yline(1,'linestyle','--','linewidth',2)
 xticklabels({'adolescent','adult'})
 ax = gca;
 ax.XAxis.FontSize = 20;
 ax.YAxis.FontSize = 20;
 movegui('east');
 box off;
 yline(0,'--','Color','k','linewidth',2);
 [p(1,1),h,stats{1,1}] = ranksum(squeeze(delta_diff_all(1,2,:)),squeeze(delta_diff_all(2,2,:)),'tail','right') ;
 [p(2,1),h,stats{2,1}] = signrank(squeeze(delta_diff_all(1,2,:))) ;
 [p(3,1),h,stats{3,1}] = signrank(squeeze(delta_diff_all(2,2,:))) ;
p

%% plot the correlation between delta bias GO and dprime 

close all
clc

for g = 1:size(GNG_rec_all_cell,2)

    for diff = [1 2]

            dprime = dprimes(g,1:numel(GNG_rec_all_cell{1,g})) ;
            delta_lick = squeeze(delta_diff_all(g,diff,1:numel(GNG_rec_all_cell{1,g})))' ;
        
          figure(diff)
        if g == 1
            scatter(dprime,delta_lick,'Marker','o','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor','none','LineWidth',2);
        elseif g == 2
            scatter(dprime,delta_lick,'Marker','o','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'LineWidth',2);
        end

        hold on
        mdl = fitlm(dprime,delta_lick)
        h =  plot(mdl)
        hold on
        delete(h(1))
          if g == 2
    h(2).Color = [.5 .5 .5] ;
        h(3).Color = [.5 .5 .5] ;
                elseif g == 1

        h(2).Color = [.5 .5 .5] ;
        h(3).Color = [.5 .5 .5] ;
                end

        h(2).LineWidth = 2 ;
        h(3).LineWidth =  2 ;
        h(2).LineStyle = L2{g} ;
        h(3).LineStyle =  L2{g} ;
        legend off
        hold on
        xlim([0 4])
        ylim([-0.5 0.75])
        xlabel("d'");
        ylabel('Delta lick rate');
        if diff ==1
            title({'R', 'Hit vs Miss'})
        elseif diff ==2
            title({'R', 'CR vs FA'})
   
        end
        box off;
        hold on
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');

        [R, P] = corr(dprime,delta_lick,'Type','Pearson','Tail','left') ;
        R_hard(g,diff) = R(2,1) ; 
        P_hard(g,diff) =  P(2,1) ;

    end
end
for diff = 1:2
    R_txt_adolescent = num2str(R_hard(1,diff))
    R_txt_adult = num2str(R_hard(2,diff))

    P_txt_adolescent = num2str(P_hard(1,diff))
    P_txt_adult = num2str(P_hard(2,diff))
    figure(diff)
    text(3.5, 0.4,R_txt_adolescent)
    text(3.5, 0.3,P_txt_adolescent)
    text(3.5, 0.2,R_txt_adult)
    text(3.5, 0.1,P_txt_adult)
end
%% plot the correlation to the c-bias 


close all
clc

for g = 1:size(GNG_rec_all_cell,2)

    for diff = [1 2]

            c_bias = cbias(g,1:numel(GNG_rec_all_cell{1,g})) ;
            delta_lick = squeeze(delta_diff_all(g,diff,1:numel(GNG_rec_all_cell{1,g})))' ;
        
          figure(diff)
        if g == 1
            scatter(c_bias,delta_lick,'Marker','o','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor','none','LineWidth',2);
        elseif g == 2
            scatter(c_bias,delta_lick,'Marker','o','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'LineWidth',2);
        end

        hold on
        mdl = fitlm(c_bias,delta_lick)
        h =  plot(mdl)
        hold on
        delete(h(1))
          if g == 2
    h(2).Color = [.5 .5 .5] ;
        h(3).Color = [.5 .5 .5] ;
                elseif g == 1

        h(2).Color = [.5 .5 .5] ;
        h(3).Color = [.5 .5 .5] ;
                end

        h(2).LineWidth = 2 ;
        h(3).LineWidth =  2 ;
        h(2).LineStyle = L2{g} ;
        h(3).LineStyle =  L2{g} ;
        legend off
        hold on
        xlim([-1.5 1.5])
        ylim([-0.5 0.75])
        xlabel("c-bias");
        ylabel('Delta lick rate');
        if diff ==1
            title({'R', 'Hit vs Miss'})
        elseif diff ==2
            title({'R', 'CR vs FA'})
   
        end
        box off;
        hold on
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');

        [R, P] = corr(c_bias,delta_lick,'Type','Pearson','Tail','left') ;
        R_hard(g,diff) = R(2,1) ; 
        P_hard(g,diff) =  P(2,1) ;

    end
end
for diff = 1:2
    R_txt_adolescent = num2str(R_hard(1,diff))
    R_txt_adult = num2str(R_hard(2,diff))

    P_txt_adolescent = num2str(P_hard(1,diff))
    P_txt_adult = num2str(P_hard(2,diff))
    figure(diff)
    text(3.5, 0.4,R_txt_adolescent)
    text(3.5, 0.3,P_txt_adolescent)
    text(3.5, 0.2,R_txt_adult)
    text(3.5, 0.1,P_txt_adult)
end
%% bias n-2
close all
clc
stim = [1 2 3 4 5 6 7];

for g = 1:size(GNG_rec_all_cell,2)
    for i   = 1:size(GNG_rec_all_cell{1,g},1) % mouse

        stim_types = GNG_rec_all_cell{1,g}(i).Behavior.stim_types;
        stim_ids =GNG_rec_all_cell{1,g}(i).Behavior.stim_ids;
        trial_responses = GNG_rec_all_cell{1,g}(i).Behavior.trial_responses;
        lick_times = GNG_rec_all_cell{1,g}(i).Behavior.lick_times;
        stim_times = GNG_rec_all_cell{1,g}(i).Behavior.stim_times;

        %hit
        idx_bias{1} = find(stim_types == 1 & trial_responses == 1 );
        idx_bias{1} = idx_bias{1}(1:end-2) +2 ;
        %miss
        idx_bias{2} =  find(stim_types == 1 & trial_responses == 0);
        idx_bias{2} = idx_bias{2}(1:end-2) +2 ;
        %fa
        idx_bias{3} = find (stim_types == -1 & trial_responses == 1);
        idx_bias{3} = idx_bias{3}(1:end-2) +2;
        %cr
        idx_bias{4} = find (stim_types == -1 & trial_responses == 0);
        idx_bias{4} = idx_bias{4}(1:end-2) +2 ;

        % all biases
        for b = 1:size(idx_bias,2)

            % easy bias easy
            index_stim = find (stim_ids == stim(1) | stim_ids == stim(4));
            if length(idx_bias{b}) > 1
                idx_intersect = intersect (index_stim,idx_bias{b}) ;
                if length(idx_intersect) > 1

                    [go_licks_all_bias_easy(g,i,b),ngo_licks_all_bias_easy(g,i,b),~,~,~,~] =....
                        GNG_lick_rt_dprime ...
                        (idx_intersect, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

                elseif length(idx_intersect) <= 1
                    go_licks_all_bias_easy(g,i,b) = nan ;
                    ngo_licks_all_bias_easy(g,i,b) = nan ;
                end
            elseif length(idx_bias{b})  <= 1
                go_licks_all_bias_easy(g,i,b) = nan ;
                ngo_licks_all_bias_easy(g,i,b) = nan ;
            end
            % hard bias
            index_stim = find (stim_ids == stim(2) | stim_ids == stim(3));
            if length(idx_bias{b}) > 1
                idx_intersect = intersect (index_stim,idx_bias{b}) ;
                if length(idx_intersect) > 1

                    [go_licks_all_bias_hard(g,i,b),ngo_licks_all_bias_hard(g,i,b),~,~,~,~] =....
                        GNG_lick_rt_dprime ...
                        (idx_intersect, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

                elseif length(idx_intersect) <= 1
                    go_licks_all_bias_hard(g,i,b) = nan ;
                    ngo_licks_all_bias_hard(g,i,b) = nan ;
                end
            elseif length(idx_bias{b})  <= 1
                go_licks_all_bias_hard(g,i,b) = nan ;
                ngo_licks_all_bias_hard(g,i,b) =nan ;
            end
        end
    end
end
%% plot the bias n-2
close all
clear licks_mice_hit licks_mice_fa licks_mice_cr licks_mice_miss  delta_diff_all
clear stats
for g = 1:size(GNG_rec_all_cell,2)

    session_ID_g = reshape(session_ID(g,:,:),[],1) ;
    session_ID_g = session_ID_g(~isnan(session_ID_g)) ;

    for b = 1:size(idx_bias,2)
        for i   = 1:numel(GNG_rec_all_cell{1,g}) % mouse

            licks_mice_hit(i,:) =  [go_licks_all_bias_easy(g,i,1) go_licks_all_bias_hard(g,i,1) ngo_licks_all_bias_hard(g,i,1) ngo_licks_all_bias_easy(g,i,1)] ;
            licks_mice_fa(i,:) =   [go_licks_all_bias_easy(g,i,3) go_licks_all_bias_hard(g,i,3) ngo_licks_all_bias_hard(g,i,3) ngo_licks_all_bias_easy(g,i,3)] ;
            licks_mice_cr(i,:) =  [go_licks_all_bias_easy(g,i,4) go_licks_all_bias_hard(g,i,4) ngo_licks_all_bias_hard(g,i,4) ngo_licks_all_bias_easy(g,i,4)] ;
            licks_mice_miss(i,:) =  [go_licks_all_bias_easy(g,i,2) go_licks_all_bias_hard(g,i,2) ngo_licks_all_bias_hard(g,i,2) ngo_licks_all_bias_easy(g,i,2)] ;

        end
    end

    licks_mice_hit(licks_mice_hit == 0) = nan ;
    licks_mice_fa(licks_mice_fa == 0) = nan ;
    licks_mice_cr(licks_mice_cr == 0) = nan ;
    licks_mice_miss(licks_mice_miss == 0) = nan ;
    

    delta_diff_all(g,1,1:numel(GNG_rec_all_cell{1,g})) =  nanmean((licks_mice_hit) - (licks_mice_miss),2)' ;
    delta_diff_all(g,2,1:numel(GNG_rec_all_cell{1,g})) =  nanmean((licks_mice_fa) - (licks_mice_cr),2)' ;
end



figure
 violinplot([ squeeze(delta_diff_all(1,1,:))....
     squeeze(delta_diff_all(2,1,:))],...
     {'adolescent ','adult '},"ViolinColor",{[.5 .5 .5; .5 .5 .5]}) ;

 xlim([0 3])
 xticks([1 2])
 ylim([-0.5 0.75])
 yticks([-0.5:0.25:0.75])
 ylabel("delta lick rate ")
 yline(1,'linestyle','--','linewidth',2)
 xticklabels({'adolescent','adult'})
 ax = gca;
 ax.XAxis.FontSize = 20;
 ax.YAxis.FontSize = 20;
 movegui('east');
 box off;
 yline(0,'--','Color','k','linewidth',2);


 [p(1,1),h,stats{1,1}] = ranksum(squeeze(delta_diff_all(1,1,:)),squeeze(delta_diff_all(2,1,:)),'tail','right') ;
 [p(2,1),h,stats{2,1}] = signrank(squeeze(delta_diff_all(1,1,:))) ;
 [p(3,1),h,stats{3,1}] = signrank(squeeze(delta_diff_all(2,1,:))) ;
p

 figure
 violinplot([ squeeze(delta_diff_all(1,2,:))....
     squeeze(delta_diff_all(2,2,:))],...
     {'adolescent ','adult '},"ViolinColor",{[.5 .5 .5; .5 .5 .5]}) ;

 xlim([0 3])
 xticks([1 2])
 ylim([-0.5 0.75])
 yticks([-0.5:0.25:0.75])
 ylabel("delta lick rate ")
 yline(1,'linestyle','--','linewidth',2)
 xticklabels({'adolescent','adult'})
 ax = gca;
 ax.XAxis.FontSize = 20;
 ax.YAxis.FontSize = 20;
 movegui('east');
 box off;
 yline(0,'--','Color','k','linewidth',2);
 [p(1,1),h,stats{1,1}] = ranksum(squeeze(delta_diff_all(1,2,:)),squeeze(delta_diff_all(2,2,:)),'tail','right') ;
 [p(2,1),h,stats{2,1}] = signrank(squeeze(delta_diff_all(1,2,:))) ;
 [p(3,1),h,stats{3,1}] = signrank(squeeze(delta_diff_all(2,2,:))) ;
p

%% auto-correlation
G = size(GNG_rec_all_cell,2);
edges = -1:0.05:1;

% Gather per-group vectors (drop NaNs)
for g = 1:G
    nM = size(GNG_rec_all_cell{1,g},1);
    vFA{g,1}  = delta_FA(g,1:nM);  vFA{g}  = vFA{g}(~isnan(vFA{g}));
    vHIT{g,1} = delta_HIT(g,1:nM); vHIT{g} = vHIT{g}(~isnan(vHIT{g}));
end

% -------- Histograms (Δ_FA)
figure; hold on
for g = 1:G
    if ~isempty(vFA{g})
        h(g) = histogram(vFA{g}, 'BinEdges', edges, 'Normalization','probability'); %#ok<SAGROW>
        h(g).DisplayStyle = 'stairs'; h(g).LineWidth = 2;
        xline(median(vFA{g}), '--', sprintf('G%d med=%.2f', g, median(vFA{g})), 'LabelVerticalAlignment','bottom');
    end
end
xline(0,'k:'); xlabel('\Delta_{FA} = P(FA|FA) - P(FA|CR)'); ylabel('Fraction of mice')
legend(arrayfun(@(x) sprintf('Group %d',x),1:G,'UniformOutput',false),'Location','best')
title('No-Go outcome autocorrelation (lag 1)')

% -------- Histograms (Δ_HIT)
figure; hold on
for g = 1:G
    if ~isempty(vHIT{g})
        hh(g) = histogram(vHIT{g}, 'BinEdges', edges, 'Normalization','probability'); %#ok<SAGROW>
        hh(g).DisplayStyle = 'stairs'; hh(g).LineWidth = 2;
        xline(median(vHIT{g}), '--', sprintf('G%d med=%.2f', g, median(vHIT{g})), 'LabelVerticalAlignment','bottom');
    end
end
xline(0,'k:'); xlabel('\Delta_{HIT} = P(Hit|Hit) - P(Hit|Miss)'); ylabel('Fraction of mice')
legend(arrayfun(@(x) sprintf('Group %d',x),1:G,'UniformOutput',false),'Location','best')
title('Go outcome autocorrelation (lag 1)')

% -------- Statistical comparisons (simple + robust)
if G == 2
    % Between groups (Δ_FA)
    [~,p_t_FA] = ttest2(vFA{1}, vFA{2});
    p_rs_FA = ranksum(vFA{1}, vFA{2});
    d_FA = (mean(vFA{1})-mean(vFA{2})) / sqrt(((numel(vFA{1})-1)*var(vFA{1}) + (numel(vFA{2})-1)*var(vFA{2}))/(numel(vFA{1})+numel(vFA{2})-2));
    % Within group vs 0 (Δ_FA)
    [~,p_t_FA_g1] = ttest(vFA{1}, 0); p_sr_FA_g1 = signrank(vFA{1}, 0);
    [~,p_t_FA_g2] = ttest(vFA{2}, 0); p_sr_FA_g2 = signrank(vFA{2}, 0);

    % Between groups (Δ_HIT)
    [~,p_t_H] = ttest2(vHIT{1}, vHIT{2});
    p_rs_H = ranksum(vHIT{1}, vHIT{2});
    d_H = (mean(vHIT{1})-mean(vHIT{2})) / sqrt(((numel(vHIT{1})-1)*var(vHIT{1}) + (numel(vHIT{2})-1)*var(vHIT{2}))/(numel(vHIT{1})+numel(vHIT{2})-2));
    % Within group vs 0 (Δ_HIT)
    [~,p_t_H_g1] = ttest(vHIT{1}, 0); p_sr_H_g1 = signrank(vHIT{1}, 0);
    [~,p_t_H_g2] = ttest(vHIT{2}, 0); p_sr_H_g2 = signrank(vHIT{2}, 0);

    % Print concise summary
    fprintf('\n=== Δ_FA (No-Go) ===\n');
    fprintf('G1 median=%.3f (n=%d), G2 median=%.3f (n=%d)\n', median(vFA{1}), numel(vFA{1}), median(vFA{2}), numel(vFA{2}));
    fprintf('Between groups: ttest2 p=%.4g, ranksum p=%.4g, Cohen d=%.2f\n', p_t_FA, p_rs_FA, d_FA);
    fprintf('Vs 0:  G1 t=%.4g / signrank=%.4g | G2 t=%.4g / signrank=%.4g\n', p_t_FA_g1, p_sr_FA_g1, p_t_FA_g2, p_sr_FA_g2);

    fprintf('\n=== Δ_HIT (Go) ===\n');
    fprintf('G1 median=%.3f (n=%d), G2 median=%.3f (n=%d)\n', median(vHIT{1}), numel(vHIT{1}), median(vHIT{2}), numel(vHIT{2}));
    fprintf('Between groups: ttest2 p=%.4g, ranksum p=%.4g, Cohen d=%.2f\n', p_t_H, p_rs_H, d_H);
    fprintf('Vs 0:  G1 t=%.4g / signrank=%.4g | G2 t=%.4g / signrank=%.4g\n', p_t_H_g1, p_sr_H_g1, p_t_H_g2, p_sr_H_g2);
end


%% binned lick rate and autocorrelation
bin_its = 10 ; 
for g = 1:size(GNG_rec_all_cell,2)
    nMice = size(GNG_rec_all_cell{1,g},1);
    FA_autocorr_delta(g,1:nMice) = NaN;     % preallocate per group

    for i = 1:nMice % mouse
        stim_types = GNG_rec_all_cell{1,g}(i).Behavior.stim_types;
        stim_ids = GNG_rec_all_cell{1,g}(i).Behavior.stim_ids;
        trial_responses = GNG_rec_all_cell{1,g}(i).Behavior.trial_responses;
        lick_times = GNG_rec_all_cell{1,g}(i).Behavior.lick_times;
        stim_times = GNG_rec_all_cell{1,g}(i).Behavior.stim_times;

        index_stim = find(stim_ids == stim(2) | stim_ids == stim(3) | stim_ids == stim(1) | stim_ids == stim(4));
        index_stim = index_stim(100:end);
        [~,~,dprimes(g,i),~,~,~,cbias(g,i)] = ...
            GNG_lick_rt_dprime(index_stim, lick_times, stim_times, stim_types, trial_responses, stim_ids, tone_dur, response_window, stim);

        % trial outcome fluctuations
        index_stim = 1:length(stim_times);
        bin_step = length(stim_times) / bin_its;
        bins = round(1:bin_step:length(stim_times)-bin_its);
        for bin = 1:length(bins)-1
            [go_licks_all(g,i,bin), ngo_licks_all(g,i,bin), ~, ~, ~, ~, ~] = ...
                GNG_lick_rt_dprime(index_stim(bins(bin):bins(bin+1)), lick_times, stim_times, ...
                stim_types, trial_responses, stim_ids, tone_dur, response_window, stim);
        end

        % easy
        index_stim = find(stim_ids == stim(1) | stim_ids == stim(4));
        [go_licks_all_easy(g,i), ngo_licks_all_easy(g,i), ~, ~, ~, ~, ~] = ...
            GNG_lick_rt_dprime(index_stim, lick_times, stim_times, stim_types, trial_responses, stim_ids, tone_dur, response_window, stim);

        % hard
        index_stim = find(stim_ids == stim(2) | stim_ids == stim(3));
        [go_licks_all_hard(g,i), ngo_licks_all_hard(g,i), ~, ~, ~, ~, ~] = ...
            GNG_lick_rt_dprime(index_stim, lick_times, stim_times, stim_types, trial_responses, stim_ids, tone_dur, response_window, stim);

        % lick raster of recording session
        eventTimes_all = GNG_rec_all_cell{1,g}(i).eventTimes.all;
        eventTimes_licks = GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_licks;

        n_trials(g,i) = length(eventTimes_all);
        % retrieve trial identities
        hit_ID{g,i}  = find(stim_types ==  1 & trial_responses == 1);
        FA_ID{g,i}   = find(stim_types == -1 & trial_responses == 1);
        miss_ID{g,i} = find(stim_types ==  1 & trial_responses == 0);
        CR_ID{g,i}   = find(stim_types == -1 & trial_responses == 0);

        % create raster matrix of all lick times across all trials
        for it = 1:length(eventTimes_all)
            licks_after_stim = lick_times(lick_times > stim_times(it));
            licks_per_trial = licks_after_stim(licks_after_stim < stim_times(it)+tone_dur+response_window);
            licks_from_stim_onset = licks_per_trial - stim_times(it);
            licks_from_stim_onset = licks_from_stim_onset(licks_from_stim_onset > 0.1);
            licks_raster{g,i}(it,1:length(licks_from_stim_onset)) = licks_from_stim_onset;
            if ismember(it, hit_ID{g,i}) == 1
                reward_raster{g,i}(it,1) = licks_from_stim_onset(1);
                if licks_from_stim_onset(1) < 0.6, reward_raster{g,i}(it,1) = 0.6; end
            end
            if ismember(it, FA_ID{g,i}) == 1
                pun_raster{g,i}(it,1) = licks_from_stim_onset(1);
                if licks_from_stim_onset(1) < 0.6, pun_raster{g,i}(it,1) = 0.6; end
            end
        end

        eventTimes_all = GNG_rec_all_cell{1,g}(i).eventTimes.all;
        eventTimes_licks = GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_licks;

        lick_time{g,i} = cell(length(eventTimes_all),1);
        for E = 1:length(eventTimes_all)
            idx_lick{E,1} = find(eventTimes_licks > (eventTimes_all(1,E) + startRange) & ...
                                 eventTimes_licks < (eventTimes_all(1,E) + stopRange));
            if ~isempty(idx_lick{E,1})
                lick_time{g,i}{E,1} = eventTimes_licks(idx_lick{E,1}) - eventTimes_all(1,E);
            end
        end

% No-Go (FA vs CR): Δ_FA = P(FA|FA) - P(FA|CR)
ng_idx = find(stim_types == -1);
x_ng = double(trial_responses(ng_idx) == 1); % 1=FA, 0=CR
if numel(x_ng) > 1
    prev = x_ng(1:end-1); nextv = x_ng(2:end);
    nFA = sum(prev==1); nCR = sum(prev==0);
    if nFA>0 && nCR>0
        pFA_FA = sum(prev==1 & nextv==1) / nFA;
        pFA_CR = sum(prev==0 & nextv==1) / nCR;
        delta_FA(g,i) = pFA_FA - pFA_CR;
        PFA_after_FA(g,i)  = pFA_FA;
        PFA_after_CR(g,i)  = pFA_CR;
    else
        delta_FA(g,i) = NaN; PFA_after_FA(g,i)=NaN; PFA_after_CR(g,i)=NaN;
    end
else
    delta_FA(g,i) = NaN; PFA_after_FA(g,i)=NaN; PFA_after_CR(g,i)=NaN;
end

% Go (Hit vs Miss): Δ_HIT = P(Hit|Hit) - P(Hit|Miss)
go_idx = find(stim_types == 1);
x_go = double(trial_responses(go_idx) == 1); % 1=Hit, 0=Miss
if numel(x_go) > 1
    prev = x_go(1:end-1); nextv = x_go(2:end);
    nHit  = sum(prev==1); nMiss = sum(prev==0);
    if nHit>0 && nMiss>0
        pHit_Hit  = sum(prev==1 & nextv==1) / nHit;
        pHit_Miss = sum(prev==0 & nextv==1) / nMiss;
        delta_HIT(g,i) = pHit_Hit - pHit_Miss;
        PHIT_after_HIT(g,i)  = pHit_Hit;
        PHIT_after_MISS(g,i) = pHit_Miss;
    else
        delta_HIT(g,i) = NaN; PHIT_after_HIT(g,i)=NaN; PHIT_after_MISS(g,i)=NaN;
    end
else
    delta_HIT(g,i) = NaN; PHIT_after_HIT(g,i)=NaN; PHIT_after_MISS(g,i)=NaN;
end

    end
end

%% plot the autocorrelaiton of FAs and Hits 
clear h 
G = size(GNG_rec_all_cell,2);
edges = -1:0.05:1;

% Gather per-group vectors (drop NaNs)
for g = 1:G
    nM = size(GNG_rec_all_cell{1,g},1);
    vFA{g,1}  = delta_FA(g,1:nM);  vFA{g}  = vFA{g}(~isnan(vFA{g}));
    vHIT{g,1} = delta_HIT(g,1:nM); vHIT{g} = vHIT{g}(~isnan(vHIT{g}));
end

figure; 
hold on
for g = 1:G
    if ~isempty(vFA{g})
        h(g) = histogram(vFA{g}, 'BinEdges', edges, 'Normalization','probability')
  
    end
end
xline(0,'k:'); 
xlabel('\Delta_{FA} = P(FA|FA) - P(FA|CR)'); 
ylabel('Fraction of mice')
title('No-Go outcome autocorrelation (lag 1)')

figure; 
hold on
for g = 1:G
    if ~isempty(vHIT{g})
        hh(g) = histogram(vHIT{g}, 'BinEdges', edges, 'Normalization','probability'); %#ok<SAGROW>
   
    end
end
xline(0,'k:');
xlabel('\Delta_{HIT} = P(Hit|Hit) - P(Hit|Miss)');
ylabel('Fraction of mice')
title('Go outcome autocorrelation (lag 1)')

if G == 2
    % Between groups (Δ_FA)
    [~,p_t_FA] = ttest2(vFA{1}, vFA{2});
    p_rs_FA = ranksum(vFA{1}, vFA{2});
    d_FA = (mean(vFA{1})-mean(vFA{2})) / sqrt(((numel(vFA{1})-1)*var(vFA{1}) + (numel(vFA{2})-1)*var(vFA{2}))/(numel(vFA{1})+numel(vFA{2})-2));
    % Within group vs 0 (Δ_FA)
    [~,p_t_FA_g1] = ttest(vFA{1}, 0); p_sr_FA_g1 = signrank(vFA{1}, 0);
    [~,p_t_FA_g2] = ttest(vFA{2}, 0); p_sr_FA_g2 = signrank(vFA{2}, 0);

    % Between groups (Δ_HIT)
    [~,p_t_H] = ttest2(vHIT{1}, vHIT{2});
    p_rs_H = ranksum(vHIT{1}, vHIT{2});
    d_H = (mean(vHIT{1})-mean(vHIT{2})) / sqrt(((numel(vHIT{1})-1)*var(vHIT{1}) + (numel(vHIT{2})-1)*var(vHIT{2}))/(numel(vHIT{1})+numel(vHIT{2})-2));
    % Within group vs 0 (Δ_HIT)
    [~,p_t_H_g1] = ttest(vHIT{1}, 0); p_sr_H_g1 = signrank(vHIT{1}, 0);
    [~,p_t_H_g2] = ttest(vHIT{2}, 0); p_sr_H_g2 = signrank(vHIT{2}, 0);

    % Print concise summary
    fprintf('\n=== Δ_FA (No-Go) ===\n');
    fprintf('G1 median=%.3f (n=%d), G2 median=%.3f (n=%d)\n', median(vFA{1}), numel(vFA{1}), median(vFA{2}), numel(vFA{2}));
    fprintf('Between groups: ttest2 p=%.4g, ranksum p=%.4g, Cohen d=%.2f\n', p_t_FA, p_rs_FA, d_FA);
    fprintf('Vs 0:  G1 t=%.4g / signrank=%.4g | G2 t=%.4g / signrank=%.4g\n', p_t_FA_g1, p_sr_FA_g1, p_t_FA_g2, p_sr_FA_g2);

    fprintf('\n=== Δ_HIT (Go) ===\n');
    fprintf('G1 median=%.3f (n=%d), G2 median=%.3f (n=%d)\n', median(vHIT{1}), numel(vHIT{1}), median(vHIT{2}), numel(vHIT{2}));
    fprintf('Between groups: ttest2 p=%.4g, ranksum p=%.4g, Cohen d=%.2f\n', p_t_H, p_rs_H, d_H);
    fprintf('Vs 0:  G1 t=%.4g / signrank=%.4g | G2 t=%.4g / signrank=%.4g\n', p_t_H_g1, p_sr_H_g1, p_t_H_g2, p_sr_H_g2);
end

%% binned lick rate over time 

Gplot = min(2, size(GNG_rec_all_cell,2));           % plot up to 2 groups
figure; t = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Precompute max #bins to define x-axis
maxBins = 0;
for g = 1:Gplot
    nM = size(GNG_rec_all_cell{1,g},1);
    if nM>0
        maxBins = max([maxBins, size(go_licks_all,3), size(ngo_licks_all,3)]);
    end
end
x = 1:maxBins;
co = lines(Gplot);

% Containers for stats
GO_mean_per_mouse  = cell(Gplot,1);
NGO_mean_per_mouse = cell(Gplot,1);
GO_slope_per_mouse  = cell(Gplot,1);
NGO_slope_per_mouse = cell(Gplot,1);

%  Go trials
nexttile; hold on
for g = 1:Gplot
    nM = size(GNG_rec_all_cell{1,g},1);
    go_mat = squeeze(go_licks_all(g,1:nM,1:maxBins));   % [mice x bins]
    % mean & SEM across mice per bin (ignore NaNs)
    mu = nanmean(go_mat,1);
    n_eff = sum(~isnan(go_mat),1);
    se = nanstd(go_mat,0,1) ./ max(1,sqrt(n_eff));
    errorbar(x, mu, se, 'LineWidth', 1.6, 'CapSize', 0, 'Color', co(g,:));

    % store per-mouse means & slopes for stats
    GO_mean_per_mouse{g} = nanmean(go_mat,2);
    % per-mouse slope across bins (simple linear trend)
    s = NaN(nM,1);
    for m = 1:nM
        yy = go_mat(m,:);
        idx = isfinite(yy);
        if sum(idx) >= 3
            p = polyfit(x(idx), yy(idx), 1);   % p(1) = slope
            s(m) = p(1);
        end
    end
    GO_slope_per_mouse{g} = s;
end
xlabel('Bin (early \rightarrow late in session)');
ylabel('Go lick rate');
title('Go trials');
ylim([0 1]);
box off

% No-Go trials
nexttile; hold on
for g = 1:Gplot
    nM = size(GNG_rec_all_cell{1,g},1);
    ngo_mat = squeeze(ngo_licks_all(g,1:nM,1:maxBins)); % [mice x bins]
    mu = nanmean(ngo_mat,1);
    n_eff = sum(~isnan(ngo_mat),1);
    se = nanstd(ngo_mat,0,1) ./ max(1,sqrt(n_eff));
    errorbar(x, mu, se, 'LineWidth', 1.6, 'CapSize', 0, 'Color', co(g,:));

    NGO_mean_per_mouse{g} = nanmean(ngo_mat,2);
    s = NaN(size(ngo_mat,1),1);
    for m = 1:size(ngo_mat,1)
        yy = ngo_mat(m,:);
        idx = isfinite(yy);
        if sum(idx) >= 3
            p = polyfit(x(idx), yy(idx), 1);
            s(m) = p(1);
        end
    end
    NGO_slope_per_mouse{g} = s;
end
xlabel('Bin (early \rightarrow late in session)');
ylabel('No-Go lick rate');
title('No-Go trials');
ylim([0 1]);
box off

for g = 1:Gplot
    s_go  = GO_slope_per_mouse{g};   s_go  = s_go(isfinite(s_go));
    s_ngo = NGO_slope_per_mouse{g};  s_ngo = s_ngo(isfinite(s_ngo));
    if ~isempty(s_go)
        [~,p_go_trend]  = ttest(s_go, 0);
    else
        p_go_trend = NaN;
    end
    if ~isempty(s_ngo)
        [~,p_ngo_trend] = ttest(s_ngo, 0);
    else
        p_ngo_trend = NaN;
    end
    fprintf('Group %d trend over bins: Go slope vs 0 p=%.4g, No-Go slope vs 0 p=%.4g\n', g, p_go_trend, p_ngo_trend);
end

if Gplot == 2
    % Go
    a = GO_mean_per_mouse{1}; b = GO_mean_per_mouse{2};
    a = a(isfinite(a)); b = b(isfinite(b));
    if ~isempty(a) && ~isempty(b)
        [~,p_go_grp] = ttest2(a,b);
        d_go = (mean(a)-mean(b)) / sqrt(((numel(a)-1)*var(a) + (numel(b)-1)*var(b)) / (numel(a)+numel(b)-2));
    else
        p_go_grp = NaN; d_go = NaN;
    end
    % No-Go
    a2 = NGO_mean_per_mouse{1}; b2 = NGO_mean_per_mouse{2};
    a2 = a2(isfinite(a2)); b2 = b2(isfinite(b2));
    if ~isempty(a2) && ~isempty(b2)
        [~,p_ngo_grp] = ttest2(a2,b2);
        d_ngo = (mean(a2)-mean(b2)) / sqrt(((numel(a2)-1)*var(a2) + (numel(b2)-1)*var(b2)) / (numel(a2)+numel(b2)-2));
    else
        p_ngo_grp = NaN; d_ngo = NaN;
    end
    fprintf('Overall Go lick rate: G1 mean=%.3f, G2 mean=%.3f, ttest2 p=%.4g, d=%.2f\n', ...
        mean(a,'omitnan'), mean(b,'omitnan'), p_go_grp, d_go);
    fprintf('Overall No-Go lick rate: G1 mean=%.3f, G2 mean=%.3f, ttest2 p=%.4g, d=%.2f\n', ...
        mean(a2,'omitnan'), mean(b2,'omitnan'), p_ngo_grp, d_ngo);
end

clearvars -except GNG_opto_all_cell
close all
clc
%% load the the opto file
addpath(genpath('Z:\Shared\Benne\Praegel_et_al_2025\analysis'))

[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2025\data\optogenetics'...
    , ['Select opto sessions ']);
addpath(path)
 load (file)

%% Parameters
tone_dur = 0.1 ; % sec
response_window = 2 ; % sec
stim = [1:16] ; % ID
freqs = [7.07 9.17 10.95 14.14] % frequencies of stim IDs

color_eh= {[0.3010 0.7450 0.9330], [0 0.4470 0.7410]} ;
color_opto = {[0.6350 0.0780 0.1840]} ;
color_opto_fade = {[0.6350 0.0780 0.1840 0.5]} ; 
color_mean_marker = {[0.5 0.5 0.5 0.5]} ;

color_plot_bias = {[0.4660 0.6740 0.1880 0.5], [0.8500 0.3250 0.0980 0.5],[0.9290 0.6940 0.1250 0.5],[0.4940 0.1840 0.5560 0.5]};
color_plot_bias_opto = {[0.4660 0.6740 0.3880 0.5], [0.8500 0.3250 0.0680 0.5],[0.9290 0.6940 0.3250 0.5], [0.4940 0.1840 0.8560 0.5]};

L = {['-'], ['-']} ;

%% analyse dprime and  lick rates 

clc
post_bias  = {'light on', 'light off'} ;
opto_go_easy = [8 1] ;
opto_go_hard = [10 3] ;
opto_ngo_easy = [9 2] ;
opto_ngo_hard = [11 4] ;


for g = 1:size(GNG_opto_all_cell,1) % transition or expert
    for i   = 1:size(GNG_opto_all_cell,2) % mouse
        for s = 1:numel(GNG_opto_all_cell{g,i}) % session

            stim_types = GNG_opto_all_cell{g,i}(s).Behavior.stim_types;
            stim_ids =GNG_opto_all_cell{g,i}(s).Behavior.stim_ids;
            trial_responses = GNG_opto_all_cell{g,i}(s).Behavior.trial_responses;
            lick_times = GNG_opto_all_cell{g,i}(s).Behavior.lick_times;
            stim_times = GNG_opto_all_cell{g,i}(s).Behavior.stim_times;

            mouse_ID(g,i,s) = i ;
            session_ID(g,i,s) = s ;

            %nature of stimulus
            for o = 1:length(post_bias)

                % easy
                index_stim = find (stim_ids == stim(opto_go_easy(o)) | stim_ids == stim(opto_ngo_easy(o)));
                [go_licks_all_easy( o,g,i,s),ngo_licks_all_easy( o,g,i,s),dprimes_easy( o,g,i,s),~,~,~] =....
                    GNG_lick_rt_dprime ...
                    (index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

                % hard
                index_stim = find (stim_ids == stim(opto_go_hard(o)) | stim_ids == stim(opto_ngo_hard(o)));
                [go_licks_all_hard(o,g,i,s),ngo_licks_all_hard( o,g,i,s),dprimes_hard(o,g,i,s),~,~,~] =....
                    GNG_lick_rt_dprime ...
                    (index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

                %biases 
                %hit
                idx_bias{1} = find(stim_types == 1 & trial_responses == 1 & stim_ids == stim(opto_go_hard(o)) | stim_ids == stim(opto_go_easy(o)));
                idx_bias{1} = idx_bias{1}(1:end-1) +1 ;
                %miss
                idx_bias{2} =  find(stim_types == 1 & trial_responses == 0 & stim_ids == stim(opto_go_hard(o)) | stim_ids == stim(opto_go_easy(o)));
                idx_bias{2} = idx_bias{2}(1:end-1) +1 ;
                %fa
                idx_bias{3} = find (stim_types == -1 & trial_responses == 1 & stim_ids == stim(opto_ngo_hard(o)) | stim_ids == stim(opto_ngo_easy(o)));
                idx_bias{3} = idx_bias{3}(1:end-1) +1 ;
                %cr
                idx_bias{4} = find (stim_types == -1 & trial_responses == 0 & stim_ids == stim(opto_ngo_hard(o)) | stim_ids == stim(opto_ngo_easy(o)));
                idx_bias{4} = idx_bias{4}(1:end-1) +1 ;

                for b = 1:size(idx_bias,2)
                    for p = 1:length(post_bias)

                        if length(idx_bias{b}) > 1

                            index_stim = find (stim_ids == stim(opto_go_easy(p)) | stim_ids == stim(opto_ngo_easy(p)));
                            idx_intersect = intersect (index_stim,idx_bias{b}) ;

                            if length(idx_intersect) > 1
                                [go_licks_all_easy_bias(p,o,g,b,i,s),ngo_licks_all_easy_bias(p,o,g,b,i,s),dprimes_easy_bias(p,o,g,b,i,s),~,~,~] =....
                                    GNG_lick_rt_dprime ...
                                    (idx_intersect, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

                            elseif length(idx_intersect) <=  1
                                go_licks_all_easy_bias(p,o,g,b,i,s) = nan ;
                                ngo_licks_all_easy_bias(p,o,g,b,i,s) = nan ;
                                dprimes_easy_bias(p,o,g,b,i,s) = nan ;
                            end

                            index_stim = find (stim_ids == stim(opto_go_hard(p)) | stim_ids == stim(opto_ngo_hard(p)));
                            idx_intersect = intersect (index_stim,idx_bias{b}) ;

                            if length(idx_intersect) > 1 

                                [go_licks_all_hard_bias(p,o,g,b,i,s),ngo_licks_all_hard_bias(p,o,g,b,i,s),dprimes_hard_bias(p,o,g,b,i,s),~,~,~] =....
                                    GNG_lick_rt_dprime ...
                                    (idx_intersect, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

                            elseif length(idx_intersect) <=  1
                                go_licks_all_hard_bias(p,o,g,b,i,s) = nan ;
                                ngo_licks_all_hard_bias(p,o,g,b,i,s) = nan ;
                                dprimes_hard_bias(p,o,g,b,i,s) = nan ;
                            end


                        elseif length(idx_bias{b}) <= 1
                            go_licks_all_easy_bias(p,o,g,b,i,s) = nan ;
                            ngo_licks_all_easy_bias(p,o,g,b,i,s) = nan ;
                            dprimes_easy_bias(p,o,g,b,i,s) = nan ;
                            go_licks_all_hard_bias(p,o,g,b,i,s) = nan ;
                            ngo_licks_all_hard_bias(p,o,g,b,i,s) = nan ;
                            dprimes_hard_bias(p,o,g,b,i,s) = nan ;
                        end
                    end
                end
            end
        end
    end
end
%% dprime 
clc
close all

for g = 1:size(GNG_opto_all_cell,1)

    dprimes_easy_g = reshape(squeeze(dprimes_easy(2,g,:,:)),[],1) ;
    dprimes_easy_g = dprimes_easy_g(~isnan(dprimes_easy_g));

    dprimes_easy_opto_g = reshape(squeeze(dprimes_easy(1,g,:,:)),[],1) ;
    dprimes_easy_opto_g = dprimes_easy_opto_g(~isnan(dprimes_easy_opto_g));

    dprimes_hard_g = reshape(squeeze(dprimes_hard(2,g,:,:)),[],1) ;
    dprimes_hard_g = dprimes_hard_g(~isnan(dprimes_hard_g));

    dprimes_hard_opto_g = reshape(squeeze(dprimes_hard(1,g,:,:)),[],1) ;
    dprimes_hard_opto_g = dprimes_hard_opto_g(~isnan(dprimes_hard_opto_g));

    dprimes_easy_g(dprimes_easy_g == 0) = nan ;
    dprimes_easy_opto_g(dprimes_easy_opto_g == 0) = nan ;

    dprimes_hard_g(dprimes_hard_g == 0) = nan ;
    dprimes_hard_opto_g(dprimes_hard_opto_g == 0) = nan ;

    figure
    plot([1 2],[mean([dprimes_easy_g dprimes_hard_g],2)  mean([dprimes_easy_opto_g dprimes_hard_opto_g],2)]...
        ,'Color',[0.5 0.5 0.5 0.2],'Marker','o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',2)
    hold on

    scatter(ones(size(dprimes_easy_opto_g)).*(1 + ones(size(dprimes_easy_opto_g))),mean([dprimes_easy_opto_g dprimes_hard_opto_g],2)...
        ,'Marker','o','MarkerFaceColor',color_opto{1},'MarkerEdgeColor',color_opto{1})
    hold on
    box off
    ylim([-0.5 3])
    ylabel("d'")
    xlim([0 3])
    xticks([1 2])
    xticklabels({'light off','light on'})
    [p,~,stats] = signrank(mean([dprimes_easy_g dprimes_hard_g],2),mean([dprimes_easy_opto_g dprimes_hard_opto_g],2),'tail','right')
            ax = gca;
        ax.XAxis.FontSize = 20 ;
        ax.YAxis.FontSize = 20 ;
        ax.Title.FontSize = 20 ;
        ax.Title.FontSize = 20 ;
end


%% psychometric curve 
clc
close all

go_licks_all_easy(go_licks_all_easy == 0) = nan ;
go_licks_all_hard(go_licks_all_hard == 0) = nan ;
ngo_licks_all_hard(ngo_licks_all_hard == 0) = nan ;
ngo_licks_all_easy(ngo_licks_all_easy == 0) = nan ;


for g =1:size(GNG_opto_all_cell,1)
    figure

    for o = 1:length(post_bias)

        for i   = 1:size(GNG_opto_all_cell,2) % mouse
            for s = 1:numel(GNG_opto_all_cell{g,i}) %ession

                if s == 1
                    go_licks_all_easy(isnan(go_licks_all_easy)) = 0 ;
                    go_licks_all_hard(isnan(go_licks_all_hard)) = 0 ;
                    ngo_licks_all_hard(isnan(ngo_licks_all_hard)) = 0 ;
                    ngo_licks_all_easy(isnan(ngo_licks_all_easy)) = 0 ;
                end

                
                if o == 1
                    plot(freqs, [go_licks_all_easy(o,g,i,s) go_licks_all_hard(o,g,i,s) ngo_licks_all_hard(o,g,i,s) ngo_licks_all_easy(o,g,i,s)]...
                        ,'Color',color_opto_fade{1},'linewidth',2,'linestyle','-')
                    hold on
                     xticks([7 9.2 10.9 14.1])
                    xlabel('freq (kHz)');
                    ylabel('lick rate');
                    set(gca, 'XDir','reverse','XScale','log');
                            ylim([0 1]);
                elseif o == 2

                    plot(freqs, [go_licks_all_easy(o,g,i,s) go_licks_all_hard(o,g,i,s) ngo_licks_all_hard(o,g,i,s) ngo_licks_all_easy(o,g,i,s)]...
                        ,'Color',color_mean_marker{1},'linewidth',2,'linestyle','-')
                    hold on
                    xticks([7 9.2 10.9 14.1])
                    xlabel('freq (kHz)');
                    ylabel('lick rate');
                    set(gca, 'XDir','reverse','XScale','log');
                            ylim([0 1]);


                end
            end

go_licks_all_easy(go_licks_all_easy == 0) = nan ;
go_licks_all_hard(go_licks_all_hard == 0) = nan ;
ngo_licks_all_hard(ngo_licks_all_hard == 0) = nan ;
ngo_licks_all_easy(ngo_licks_all_easy == 0) = nan ;



            mean_licks_mice(i,:) =  [nanmean(go_licks_all_easy(o,g,i,:)) nanmean(go_licks_all_hard(o,g,i,:)) nanmean(ngo_licks_all_hard(o,g,i,:)) nanmean(ngo_licks_all_easy(o,g,i,:))] ;
            std_licks_mice(i,:) =  [nanstd(go_licks_all_easy(o,g,i,:)) nanstd(go_licks_all_hard(o,g,i,:)) nanstd(ngo_licks_all_hard(o,g,i,:)) nanstd(ngo_licks_all_easy(o,g,i,:))] ;



        end

        mean_licks = mean(mean_licks_mice) ;
        std_licks = mean(std_licks_mice) ;

        session_ID_g = reshape(session_ID(g,:,:),[],1) ;

        if o == 1

            errorbar(freqs, mean_licks ,std_licks/...
                sqrt(length(session_ID_g)),'Color',color_opto_fade{1},'linestyle',L{1},'linewidth',5)
            hold on
        elseif o == 2

            errorbar(freqs, mean_licks ,std_licks/...
                sqrt(length(session_ID_g)),'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)
            hold on
        end

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
% stats corrected for go and no go 
[p,~,stats] = signrank(reshape( [reshape(go_licks_all_easy(1,1,:,:),[],1) reshape(go_licks_all_hard(1,1,:,:),[],1)],[],1),...
    reshape( [reshape(go_licks_all_easy(2,1,:,:),[],1) reshape(go_licks_all_hard(2,1,:,:),[],1)],[],1))
p = p*2 


[p,~,stats] = signrank(reshape( [reshape(ngo_licks_all_easy(1,1,:,:),[],1) reshape(ngo_licks_all_hard(1,1,:,:),[],1)],[],1),...
    reshape( [reshape(ngo_licks_all_easy(2,1,:,:),[],1) reshape(ngo_licks_all_hard(2,1,:,:),[],1)],[],1))

p = p*2  




%% lick rate after no go biases
clc
close all

go_licks_all_easy_bias(go_licks_all_easy_bias == 0) = nan ;
go_licks_all_hard_bias(go_licks_all_hard_bias == 0) = nan ;
ngo_licks_all_hard_bias(ngo_licks_all_hard_bias == 0) = nan ;
ngo_licks_all_easy_bias(ngo_licks_all_easy_bias == 0) = nan ;


post_bias  = {'light on', 'light off'} ;
%easy
for g =1:size(GNG_opto_all_cell,1)

    for o = 1:length(post_bias)
figure
        for p = 1:length(post_bias)

            for b = [ 3 4]%1:size(idx_bias,2)

                licks_mice_opto_bias(p,o,g,b,:,:) = [reshape(squeeze(go_licks_all_easy_bias(p,o,g,b,:,:)),[],1), reshape(squeeze(go_licks_all_hard_bias(p,o,g,b,:,:)),[],1),...
                    reshape(squeeze(ngo_licks_all_hard_bias(p,o,g,b,:,:)),[],1), reshape(squeeze(ngo_licks_all_easy_bias(p,o,g,b,:,:)),[],1)] ;

                mean_licks_mice_opto_bias(p,o,g,b,:) = nanmean([reshape(squeeze(go_licks_all_easy_bias(p,o,g,b,:,:)),[],1), reshape(squeeze(go_licks_all_hard_bias(p,o,g,b,:,:)),[],1),...
                    reshape(squeeze(ngo_licks_all_hard_bias(p,o,g,b,:,:)),[],1), reshape(squeeze(ngo_licks_all_easy_bias(p,o,g,b,:,:)),[],1)]) ;

                stderr_licks_mice_opto_bias(p,o,g,b,:) = nanstd([reshape(squeeze(go_licks_all_easy_bias(p,o,g,b,:,:)),[],1), reshape(squeeze(go_licks_all_hard_bias(p,o,g,b,:,:)),[],1),...
                    reshape(squeeze(ngo_licks_all_hard_bias(p,o,g,b,:,:)),[],1), reshape(squeeze(ngo_licks_all_easy_bias(p,o,g,b,:,:)),[],1)])/sqrt(5) ;

                subplot(1,2,p)
                errorbar(freqs, reshape(mean_licks_mice_opto_bias(p,o,g,b,:),[],4) ,reshape(stderr_licks_mice_opto_bias(p,o,g,b,:),[],4)...
                    ,'Color',color_plot_bias{b},'linestyle',L{1},'linewidth',2)
                hold on
                ylim([0 1]);
                yticks([0:0.2:1])
                yline(0.5,'--','Color','k','linewidth',2);
                xline(10,'--','Color','k','linewidth',2);
                xlim([6 16])
                xticks([7 9.2 10.9 14.1])
                xlabel('freq (kHz)');
                ylabel('lick rate');
                title({post_bias{p},'after' post_bias{o}})
                set(gca, 'XDir','reverse','XScale','log');
                box off;
                hold on
                ax = gca;
                ax.XAxis.FontSize = 20 ;
                ax.YAxis.FontSize = 20 ;
                ax.Title.FontSize = 20 ;
                ax.Title.FontSize = 20 ;



                %% Repeated-measures ANOVA for b = 3 and b = 4 (same code block)
% Within-subject factor: Condition (4 levels)
% Conditions: GoEasy, GoHard, NoGoHard, NoGoEasy
% Runs two separate RM-ANOVAs: one for b=3 and one for b=4


    % --- 1) Build subjects x conditions matrix ---
    Y = [ ...
        reshape(squeeze(go_licks_all_easy_bias(p,o,g,b,:,:)), [], 1), ...
        reshape(squeeze(go_licks_all_hard_bias(p,o,g,b,:,:)), [], 1), ...
        reshape(squeeze(ngo_licks_all_hard_bias(p,o,g,b,:,:)), [], 1), ...
        reshape(squeeze(ngo_licks_all_easy_bias(p,o,g,b,:,:)), [], 1) ...
        ];

    % Sanity checks
    assert(size(Y,2) == 4, 'Expected 4 within-subject conditions.');
    if any(~isfinite(Y(:)))
        warning('b=%d: Found NaN/Inf. Removing subjects (rows) with any non-finite values.', b);
        bad = any(~isfinite(Y), 2);
        Y = Y(~bad, :);
    end
    assert(size(Y,1) >= 2, 'b=%d: Not enough subjects after cleanup.', b);

    % --- 2) Create wide table (required by fitrm) ---
    T = array2table(Y, 'VariableNames', {'GoEasy','GoHard','NoGoHard','NoGoEasy'});
    T.Subject = (1:height(T))';
    T = movevars(T, 'Subject', 'Before', 1);

    % --- 3) Within-subject design ---
    WithinDesign = table( ...
        categorical({'GoEasy'; 'GoHard'; 'NoGoHard'; 'NoGoEasy'}), ...
        'VariableNames', {'Condition'});

    % --- 4) Fit RM model and run ANOVA ---
    rm = fitrm(T, 'GoEasy-NoGoEasy ~ 1', 'WithinDesign', WithinDesign);

    ranovatbl = ranova(rm, 'WithinModel', 'Condition');



            end
        end
    end
end


%%  delta bias after light on and light off - clean
clc
close all
clear p change_lick_all
go_licks_all_easy_bias(go_licks_all_easy_bias == 0) = nan ;
go_licks_all_hard_bias(go_licks_all_hard_bias == 0) = nan ;
ngo_licks_all_hard_bias(ngo_licks_all_hard_bias == 0) = nan ;
ngo_licks_all_easy_bias(ngo_licks_all_easy_bias == 0) = nan ;



for g = 1:size(GNG_opto_all_cell,1)
    figure
    for p = 1:length(post_bias)

if p == 1
    o = 1 ; 
elseif p == 2
    o = 2; 
end 

        licks_bias_fa =   reshape(licks_mice_opto_bias(p,o,g,3,:,:),[],4) ;
        licks_bias_cr = reshape(licks_mice_opto_bias(p,o,g,4,:,:),[],4) ;


        change_lick = nanmean((licks_bias_fa - licks_bias_cr),2) ;
        change_lick_all(:,p) = change_lick ; 

    end 

    figure
    plot([1 2],change_lick_all(:,:)...
        ,'Color',[0.5 0.5 0.5 0.2],'Marker','o','MarkerFaceColor',color_opto{1},'MarkerEdgeColor',color_opto{1},'LineWidth',2)
    hold on

    scatter(ones(size(change_lick_all(:,2))).*(1 + ones(size(change_lick_all(:,1)))),change_lick_all(:,2)...
        ,'Marker','o','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor',[.5 .5 .5])
    hold on


        ylim([-0.25 0.5]);
        yticks([-0.25:0.25:0.5])
        yline(0,'--','Color','k','linewidth',2);
        xline(10,'--','Color','k','linewidth',2);
        xlim([0 3])
        xticks([1 2])
        xticklabels({'light on','light off'});
        ylabel('delta lick rate');
        box off;
        hold on
        ax = gca;
        ax.XAxis.FontSize = 20 ;
        ax.YAxis.FontSize = 20 ;
        ax.Title.FontSize = 20 ;
        ax.Title.FontSize = 20 ;
    end
                [p_avg,~,stats] = signrank( change_lick_all(:,1), change_lick_all(:,2),'tail','right')



%%  delta bias after light on and light off - unclean - light on after light off
clc
close all
clear p change_lick_all
go_licks_all_easy_bias(go_licks_all_easy_bias == 0) = nan ;
go_licks_all_hard_bias(go_licks_all_hard_bias == 0) = nan ;
ngo_licks_all_hard_bias(ngo_licks_all_hard_bias == 0) = nan ;
ngo_licks_all_easy_bias(ngo_licks_all_easy_bias == 0) = nan ;



for g = 1:size(GNG_opto_all_cell,1)
    figure
    for p = 1:length(post_bias)

if p == 1
    o = 2 ; 
elseif p == 2
    o = 2; 
end 

        licks_bias_fa =   reshape(licks_mice_opto_bias(p,o,g,3,:,:),[],4) ;
        licks_bias_cr = reshape(licks_mice_opto_bias(p,o,g,4,:,:),[],4) ;


        change_lick = nanmean((licks_bias_fa - licks_bias_cr),2) ;
        change_lick_all(:,p) = change_lick ; 

    end 

    figure
    plot([1 2],change_lick_all(:,:)...
        ,'Color',[0.5 0.5 0.5 0.2],'Marker','o','MarkerFaceColor',color_opto{1},'MarkerEdgeColor',color_opto{1},'LineWidth',2)
    hold on

    scatter(ones(size(change_lick_all(:,2))).*(1 + ones(size(change_lick_all(:,1)))),change_lick_all(:,2)...
        ,'Marker','o','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor',[.5 .5 .5])
    hold on


        ylim([-0.25 0.5]);
        yticks([-0.25:0.25:0.5])
        yline(0,'--','Color','k','linewidth',2);
        xline(10,'--','Color','k','linewidth',2);
        xlim([0 3])
        xticks([1 2])
        xticklabels({'light on','light off'});
        ylabel('delta lick rate');
        box off;
        hold on
        ax = gca;
        ax.XAxis.FontSize = 20 ;
        ax.YAxis.FontSize = 20 ;
        ax.Title.FontSize = 20 ;
        ax.Title.FontSize = 20 ;
    end
                [p_avg,~,stats] = signrank( change_lick_all(:,1), change_lick_all(:,2),'tail','right')

%%  delta bias after light on and light off - unclean - light off after light on
clc
close all
clear p change_lick_all
go_licks_all_easy_bias(go_licks_all_easy_bias == 0) = nan ;
go_licks_all_hard_bias(go_licks_all_hard_bias == 0) = nan ;
ngo_licks_all_hard_bias(ngo_licks_all_hard_bias == 0) = nan ;
ngo_licks_all_easy_bias(ngo_licks_all_easy_bias == 0) = nan ;



for g = 1:size(GNG_opto_all_cell,1)
    figure
    for o = 1:length(post_bias)

if o == 1
    p = 2 ; 
elseif o == 2
    p = 2; 
end 

        licks_bias_fa =   reshape(licks_mice_opto_bias(p,o,g,3,:,:),[],4) ;
        licks_bias_cr = reshape(licks_mice_opto_bias(p,o,g,4,:,:),[],4) ;


        change_lick = nanmean((licks_bias_fa - licks_bias_cr),2) ;
        change_lick_all(:,o) = change_lick ; 

    end 

    figure
    plot([1 2],change_lick_all(:,:)...
        ,'Color',[0.5 0.5 0.5 0.2],'Marker','o','MarkerFaceColor',color_opto{1},'MarkerEdgeColor',color_opto{1},'LineWidth',2)
    hold on

    scatter(ones(size(change_lick_all(:,2))).*(1 + ones(size(change_lick_all(:,1)))),change_lick_all(:,2)...
        ,'Marker','o','MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor',[.5 .5 .5])
    hold on


        ylim([-0.25 0.5]);
        yticks([-0.25:0.25:0.5])
        yline(0,'--','Color','k','linewidth',2);
        xline(10,'--','Color','k','linewidth',2);
        xlim([0 3])
        xticks([1 2])
        xticklabels({'light on','light off'});
        ylabel('delta lick rate');
        box off;
        hold on
        ax = gca;
        ax.XAxis.FontSize = 20 ;
        ax.YAxis.FontSize = 20 ;
        ax.Title.FontSize = 20 ;
        ax.Title.FontSize = 20 ;
    end
                [p_avg,~,stats] = signrank( change_lick_all(:,1), change_lick_all(:,2),'tail','right')


%% Lick  parameters
startRange_resp = - 0.1 ;
stopRange_resp = 2 ;
binSize = 0.001 ;
smoothSize = 20 ;

binBorders = startRange_resp:binSize:stopRange_resp ; %define the time points
numBins = length(binBorders)-1 ;

stim_IDs = [1 2 3 4 8 9 10 11]
%% analyse lick times 

mouse_ID = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;
session_ID = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;

for g = 1:size(GNG_opto_all_cell,1)
    for i   = 1:size(GNG_opto_all_cell,2) % mouse
        for s = 1:numel(GNG_opto_all_cell{g,i}) % session
            
            stim_types = GNG_opto_all_cell{g,i}(s).Behavior.stim_types;
            stim_ids =GNG_opto_all_cell{g,i}(s).Behavior.stim_ids;
            trial_responses = GNG_opto_all_cell{g,i}(s).Behavior.trial_responses;
            lick_times = GNG_opto_all_cell{g,i}(s).Behavior.lick_times;
            stim_times = GNG_opto_all_cell{g,i}(s).Behavior.stim_times;
            
            % easy go
            for stim_ID = 1:length(stim_IDs)
                index_stim = find (stim_ids == stim(stim_IDs(stim_ID)));
                binArray = zeros(length(index_stim), numBins) ;
                
                for r = 1:length(index_stim)
                    [n,binCenters] = histdiff(lick_times(), stim_times(index_stim(r)), binBorders); %get the spike per binsize
                    binArray(r,:) = n;


                 reaction_times(g,i,s,stim_ID,r) =  lick_times(index_stim(r)) ;
                end
                inclRange = binCenters > startRange_resp & binCenters<= stopRange_resp;
                lickCounts = sum(binArray(:,inclRange),2)./(stopRange_resp - startRange_resp); % sum the spikecounts
                
                gaussian_window = gausswin(round(smoothSize*6),3); %apply a gaussian window
                summed_window = gaussian_window./sum(gaussian_window); %sum all winows
                
                binArray_smoothed = conv2(summed_window,1,binArray', 'same')'./binSize;  %conv the summed window
                binArray_lick_all{g,i,s,stim_ID}  = binArray ;
                psth_lick(g,i,s,stim_ID,:) = mean(binArray./binSize) ; % normalize to Hz
                frArray_lick_all{g,i,s,stim_ID} = binArray_smoothed ;
                psth_lick_Smoothed(g,i,s,stim_ID,:) = mean(binArray_smoothed);
            end
        end
    end
end
%% lick plot parameters
tone_dur = 0.1 ; % sec
response_window = 2 ; % sec
stim = [1:16] ; % ID
freqs = [7.07 9.17 10.95 14.14] % frequencies of stim IDs

color_eh= {[0.3010 0.7450 0.9330], [0 0.4470 0.7410]} ;
color_opto_eh = {[0.8500 0.3250 0.0980],[0.6350 0.0780 0.1840]} ;
color_opto = {[0.6350 0.0780 0.1840]} ; 
color_opto_fade = {[0.6350 0.0780 0.1840 0.5]} ;
color_mean_marker = {[0.5 0.5 0.5 0.5]} ;

color_plot_bias = {[0.4660 0.6740 0.1880 0.5],[0.9290 0.6940 0.1250 0.5], [0.8500 0.3250 0.0980 0.5], [0.4940 0.1840 0.5560 0.5]};
color_plot_bias_opto = {[0.4660 0.6740 0.3880 0.5],[0.9290 0.6940 0.3250 0.5], [0.8500 0.3250 0.0680 0.5], [0.4940 0.1840 0.8560 0.5]};

L = {['-'], ['-']} ;

%% lick psth for go stimuli


    figure(g)
    % go easy
    mean_psth_smoothed_stimID_easy = squeeze(mean(mean(squeeze(psth_lick_Smoothed(g,:,:,1,:)))))'  ;
    stderr_psth_smoothed_stimID_easy = squeeze(std(std(squeeze(psth_lick_Smoothed(g,:,:,1,:)))))' / sqrt(7) ;
    % go hard
    mean_psth_smoothed_stimID_hard = squeeze(mean(mean(squeeze(psth_lick_Smoothed(g,:,:,3,:)))))'  ;
    stderr_psth_smoothed_stimID_hard = squeeze(std(std(squeeze(psth_lick_Smoothed(g,:,:,3,:)))))' / sqrt(7) ;
    % go easy opto
    mean_psth_smoothed_stimID_easy_opto = squeeze(mean(mean(squeeze(psth_lick_Smoothed(g,:,:,5,:)))))'  ;
    stderr_psth_smoothed_stimID_easy_opto = squeeze(std(std(squeeze(psth_lick_Smoothed(g,:,:,5,:)))))' / sqrt(7) ;
    % go hard opto
    mean_psth_smoothed_stimID_hard_opto = squeeze(mean(mean(squeeze(psth_lick_Smoothed(g,:,:,7,:)))))'  ;
    stderr_psth_smoothed_stimID_hard_opto = squeeze(std(std(squeeze(psth_lick_Smoothed(g,:,:,7,:)))))' / sqrt(7) ;


    mean_psth_smoothed_stimID_non = mean([mean_psth_smoothed_stimID_easy; mean_psth_smoothed_stimID_hard],1)
    stderr_psth_smoothed_stimID = mean([stderr_psth_smoothed_stimID_easy; stderr_psth_smoothed_stimID_hard],1)
    plot(mean_psth_smoothed_stimID,'Color' ,[.3 .3 .3])
    hold on
    patch([[1,1:2100-1] flip([1,1:2100-1])] , [mean_psth_smoothed_stimID_non + ....
        stderr_psth_smoothed_stimID flip(mean_psth_smoothed_stimID_non - stderr_psth_smoothed_stimID)],[.3 .3 .3] ,...
        'facealpha' , 0.2, 'EdgeColor','none')
    hold on


    mean_psth_smoothed_stimID_opto = mean([mean_psth_smoothed_stimID_easy_opto; mean_psth_smoothed_stimID_hard_opto],1)
    stderr_psth_smoothed_stimID = mean([stderr_psth_smoothed_stimID_easy_opto; stderr_psth_smoothed_stimID_hard_opto],1)
    plot(mean_psth_smoothed_stimID,'Color' ,color_opto_eh{2})
    hold on
    patch([[1,1:2100-1] flip([1,1:2100-1])] , [mean_psth_smoothed_stimID_opto + ....
        stderr_psth_smoothed_stimID flip(mean_psth_smoothed_stimID_opto - stderr_psth_smoothed_stimID)],color_opto_eh{2} ,...
        'facealpha' , 0.2, 'EdgeColor','none')
    hold on

    xlim([0 2100])
    xticks([100  700 2100])
    xticklabels([0  600 2000])
    xline(100,'k--')
    xline(200,'k--')
    xline(700,'k--')
    xlabel('time(ms)')
    ylabel('Lick (Hz)')
    ylim([0 10])
    box off
    ax = gca;
    ax.XAxis.FontSize = 20 ;
    ax.YAxis.FontSize = 20 ;
    ax.Title.FontSize = 20 ;
    ax.Title.FontSize = 20 ;


[p,~,stats] = signrank(mean_psth_smoothed_stimID_opto, mean_psth_smoothed_stimID_non)



%% reation time 
clc
close all

for g = 1:size(GNG_opto_all_cell,1)

    for i   = 1:size(GNG_opto_all_cell,2) % mouse
        for s = 1:numel(GNG_opto_all_cell{g,i}) % session

            % easy go
            for stim_ID = 1:length(stim_IDs)
                mean_r(g,i,s,stim_ID,:) = mean(reaction_times(g,i,s,stim_ID,:))


            end
        end
    end


    mean_r_non = (reshape(squeeze(squeeze(mean_r(g,:,:,5,:))),[],1) ) ; 
    mean_r_non(mean_r_non == 0) = nan ; 

    mean_r_opto = (reshape(squeeze(squeeze(mean_r(g,:,:,2,:))),[],1) ) ; 
    mean_r_opto(mean_r_opto == 0) = nan ; 


    figure
    plot([1 2],[mean_r_non  mean_r_opto]...
        ,'Color',[0.5 0.5 0.5 0.2],'Marker','o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',2)
    hold on

    scatter(ones(size(mean_r_opto)).*(1 + ones(size(mean_r_opto))),mean_r_opto...
        ,'Marker','o','MarkerFaceColor',color_opto{1},'MarkerEdgeColor',color_opto{1})
    hold on
    box off
    ylabel("time(ms)")
    xlim([0 3])
    xticks([1 2])
    xticklabels({'light off','light on'})
    [p,~,stats] = signrank(mean_r_non,mean_r_opto,'tail','right')
            ax = gca;
        ax.XAxis.FontSize = 20 ;
        ax.YAxis.FontSize = 20 ;
        ax.Title.FontSize = 20 ;
        ax.Title.FontSize = 20 ;
end


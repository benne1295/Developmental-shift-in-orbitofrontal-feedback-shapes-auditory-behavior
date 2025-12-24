%% mean fluorescence MGreen Lantern ACx
clear all
close all
clc
n_groups = 2 ;
colors_ado_adu = {[.5 .7 .2],[ .2 .4 .2]} ; 
%% mean fluorescence MGreen Lantern ACx

all_norm_mean = nan(5,n_groups) ;
all_norm_intDen = nan(5,n_groups) ;

for group = 1:n_groups
    if group == 1
        n_brains = 5 ;
    elseif group == 2
        n_brains = 4 ;
    end
    for brain = 1:n_brains
        dname = uigetdir('Z:\Shared\Benne\Feedback_Anatomy\MGreen_Lantern');

        dirname=fullfile(dname);
        cd(dirname)
        saving_folder = dirname;
        data = readtable('Results.xlsx') ;

        group
        brain + 1
        cell_fluoresc_table_acx{brain,group} = data ;
    end
end

%% mean fluorescence MGreen Lantern OFC
clc
n_groups = 2 ;
all_norm_mean_ofc = nan(5,n_groups) ;
all_norm_intDen_ofc = nan(5,n_groups) ;
for group = 1:n_groups
    if group == 1
        n_brains = 5 ;
    elseif group == 2
        n_brains = 4 ;
    end
    for brain = 1:n_brains
        dname = uigetdir('Z:\Shared\Benne\Feedback_Anatomy\MGreen_Lantern');

        dirname=fullfile(dname);
        cd(dirname)
        saving_folder = dirname;
        data = readtable('Results.xlsx') ;

        group
        brain + 1
        cell_fluoresc_table_ofc{brain,group} = data ;
    end
end


%% mean fluorescence and density of fluorescent signal in ACx
close all
all_norm_mean_acx = nan(5,n_groups) ;
all_norm_intDen_acx = nan(5,n_groups) ;
all_norm_mean = nan(5,n_groups) ;
all_norm_intDen = nan(5,n_groups) ;

for group = 1:n_groups
    if group == 1
        n_brains = 5 ;
    elseif group == 2
        n_brains = 4 ;
    end
    for brain = 1:n_brains



          meanColumn =  cell_fluoresc_table_ofc{brain,group}.Mean;
        intDenColumn =  cell_fluoresc_table_ofc{brain,group}.IntDen;

        % Separate odd and even rows
        oddRows = 1:2:height(cell_fluoresc_table_ofc{brain,group}); % Odd-indexed rows (injection side)
        evenRows = 2:2:height(cell_fluoresc_table_ofc{brain,group}); % Even-indexed rows (contralateral side)

        norm_mean =  abs(meanColumn(oddRows) - meanColumn(evenRows)) ;
        norm_intDen =  abs(intDenColumn(oddRows) - intDenColumn(evenRows)) ;

        all_norm_mean_ofc (brain,group) = nanmean(norm_mean) ;
        all_norm_intDen_ofc(brain,group) = nanmean(norm_intDen) ;

        meanColumn =  cell_fluoresc_table_acx{brain,group}.Mean;
        intDenColumn =  cell_fluoresc_table_acx{brain,group}.IntDen;

        % Separate odd and even rows
        oddRows = 1:2:height(cell_fluoresc_table_acx{brain,group}); % Odd-indexed rows (injection side)
        evenRows = 2:2:height(cell_fluoresc_table_acx{brain,group}); % Even-indexed rows (contralateral side)

        norm_mean =  abs(meanColumn(oddRows) - meanColumn(evenRows)) ;
        norm_intDen =  abs(intDenColumn(oddRows) - intDenColumn(evenRows)) ;

        all_norm_mean_acx (brain,group) = nanmean(norm_mean) ;
        all_norm_intDen_acx(brain,group) = nanmean(norm_intDen) ;

        

 all_norm_mean(brain,group) = all_norm_mean_acx(brain,group) / all_norm_mean_ofc(brain,group) ;

 all_norm_intDen(brain,group) = all_norm_intDen_acx(brain,group) / all_norm_intDen_ofc(brain,group) ;

    end
end
%%
  figure
violinplot([all_norm_mean(:,1), all_norm_mean(:,2)],...
    {'adolescent ','adult '},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;
box off
ylabel('ratio of mean fluorescence')
% ylim([0 .3])
% yticks([0 .1 .2 .3])
box off;
title('ACx ')
hold on;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

[p,h,stats] = ranksum(all_norm_mean(:,1),all_norm_mean(:,2),'tail','right')


  figure
violinplot([all_norm_mean_acx(:,1), all_norm_mean_acx(:,2)],...
    {'adolescent ','adult '},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;
box off
ylabel(' mean fluorescence')
% ylim([0 .3])
% yticks([0 .1 .2 .3])
box off;
title('ACx ')
hold on;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

[p,h,stats] = ranksum(all_norm_mean_acx(:,1),all_norm_mean_acx(:,2))


  figure
violinplot([all_norm_mean_ofc(:,1), all_norm_mean_ofc(:,2)],...
    {'adolescent ','adult '},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;
box off
ylabel('mean fluorescence')
% ylim([0 .3])
% yticks([0 .1 .2 .3])
box off;
title('OFC ')
hold on;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

[p,h,stats] = ranksum(all_norm_mean_ofc(:,1),all_norm_mean_ofc(:,2))
%% feedback cell count 
clearvars -except cell_table_all area_table_all sum_area sum_brain
close all
clc
colors_ado_adu = {[.5 .7 .2],[ .2 .4 .2]} ; 
n_groups = 2 ;
colors_ado_adu = {[.5 .7 .2],[ .2 .4 .2]} ; 

%% Parameters
for group = 1:n_groups
    if group == 1
    n_brains = 6 ;
    elseif group == 2
        n_brains = 6 ;
    end 
    for brain = 1:n_brains
        %% define directory
        dname = uigetdir('Z:\Shared\Benne_Adolescence\Reto-AAV- feedback tracing');
        
        % making a folder to save images 
        dirname=fullfile(dname);
        cd(dirname)
        saving_folder = dirname;
        cell_table = readtable('cell_count.txt') ;
        %% reorder table
        cell_table.slice_number = cell_table.Var2 ;
        cell_table.area = cell_table.Var3 ;
        cell_table.cell_count = cell_table.Var5 ;
        
        cell_table.Var1 = [] ;
        cell_table.Var2 = [] ;
        cell_table.Var3 = [] ;
        cell_table.Var4 = [] ;
        cell_table.Var5 = [] ;
        
        cell_table_all{1,group}{1,brain} = cell_table ;
        
        for area = 1:length(cortex)
            index = contains(cell_table.area, cortex(area)) ;
            area_table = cell_table(index == 1,:) ;
            
            area_table.slice_number = char(area_table.slice_number) ;
            area_table.area = char(area_table.area) ;
            
            area_table_all{1,group}{area,brain} = area_table ;
            sum_area{1,group}(area,brain) = sum(area_table.cell_count) ;
        end
        sum_brain{1,group}(:,brain) = sum(sum_area{1,group}(:,brain));
        
         
    end
end
%% caclulate 
      
        for group = 1:n_groups
         n_brains =  size(area_table_all{1,group},2) ;
            
            for brain = 1:n_brains
                for area = 1:length(cortex)
                    
                    %percent area to overall A1
                    perc_area_A1{1,group}(area,brain) =  sum_area{1,group}(area,brain)/sum_area{1,group}(1,brain);
                    
                    perc_area_mean_A1{1,group}(area,brain) =  sum_area{1,group}(area,brain)/mean(sum_area{1,group}(1,:));
                    
                    %percent area to all counted cells
                    perc_total_area{1,group}(area,brain) =  sum_area{1,group}(area,brain)/sum_brain{1,group}(:,brain) ;
                   
                    perc_total_area_mean{1,group}(area,brain) =  sum_area{1,group}(area,brain)/mean(sum_brain{1,group}(:,:)) ;

                end
            end
        end



%% ratio of mgb and ofc cells to ACx 

close all

for g = 1:n_groups

    MGB{1,g} = perc_area_A1{1,g}(5,:) ;
    OFC{1,g} = mean(perc_area_A1{1,g}(6:8,:)) ;
end

figure
violinplot([MGB{1,1}', MGB{1,2}'],...
    {'adolescent ','adult '},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;
box off
ylabel('#frac to A1')
ylim([0 .3])
yticks([0 .1 .2 .3])
box off;
title('MGB')
hold on;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

[p_mgb,~,stats_mgb] = ranksum( MGB{1,1} ,  MGB{1,2},'tail', 'both', 'alpha', 0.05)
ylabel('#frac to A1')
title('MGB')
box off;

figure
violinplot([OFC{1,1}', OFC{1,2}'],...
    {'adolescent ','adult '},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;
ylim([0 .3])
yticks([0 .1 .2 .3])
box off
ylabel('#frac to A1')
box off;
title('OFC')
hold on;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

[p_ofc,~,stats_ofc] = ranksum( OFC{1,1} ,  OFC{1,2},'tail', 'both', 'alpha', 0.05)

%% feedback ratio 

figure
for g = 1:n_groups

   feedback_ratio{1,g} =  (MGB{1,g} ./ OFC{1,g}) ;
end 

violinplot([feedback_ratio{1,1}', feedback_ratio{1,2}'],...
    {'adolescent ','adult '},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;
ylim([0 5])
yticks([0 2 4])

ylabel('#ratio of MGb to OFC fraction')
box off;
title('ratio of MGB - OFC')
 ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');
[p_feedback,~,~] = ranksum( feedback_ratio{1,1} ,  feedback_ratio{1,2},'tail', 'both', 'alpha', 0.05) 


%% total cell number mgb and ofc 

close all
for g = 1:n_groups
    MGB{1,g} = sum_area{1,g}(5,:) ;
    OFC{1,g} = mean(sum_area{1,g}(6:8,:)) ;
    ACx{1,g} =  sum_area{1,g}(1,:) ;

end 

figure
violinplot([MGB{1,1}', MGB{1,2}'],...
    {'adolescent ','adult '},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;

[p_mgb,~,stats_mgb] = ranksum( MGB{1,1} ,  MGB{1,2},'tail', 'both', 'alpha', 0.05)
ylabel('number of cells')
ylim([0 1000])
yticks([0:250:1000])
xticks([1 2])
xticklabels({'adolescent','adult'})
box off;
title('MGB')
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

figure
violinplot([OFC{1,1}', OFC{1,2}'],...
    {'adolescent ','adult '},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;

[p_ofc,~,stats_ofc] = ranksum( OFC{1,1} ,  OFC{1,2},'tail', 'both', 'alpha', 0.05)
ylabel('number of cells')
ylim([0 1000])
yticks([0:250:1000])
xticks([1 2])
xticklabels({'adolescent','adult'})
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;
title('OFC')

figure
violinplot([ACx{1,1}', ACx{1,2}'],...
    {'adolescent ','adult '},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;

[p_acx,~,stats_acx] = ranksum( ACx{1,1} ,  ACx{1,2},'tail', 'both', 'alpha', 0.05)
ylabel('number of cells')
xticks([1 2])
xticklabels({'adolescent','adult'})
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;
title('ACx')

        
    
%% feedback cell number 

close all
figure
for g = 1:n_groups

    feedback_ratio{1,g} =  (MGB{1,g} ./ OFC{1,g})  ;
end
violinplot([feedback_ratio{1,1}', feedback_ratio{1,2}'],...
    {'adolescent ','adult '},"ViolinColor",{[colors_ado_adu{1}; colors_ado_adu{2}]} ) ;
ylabel('ratio of cells')
ylim([0 5])
yticks([0:1:4])
xticks([1 2])
xticklabels({'adolescent','adult'})
title('n MGB / n OFC')

box off;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
[p_feedback,~,~] = ranksum( feedback_ratio{1,1} ,  feedback_ratio{1,2},'tail', 'both', 'alpha', 0.05)
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');



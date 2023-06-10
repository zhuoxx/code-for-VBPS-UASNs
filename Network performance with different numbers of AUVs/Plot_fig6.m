clc;clear;close all

path1 = 'C:\Users\zxx\Desktop\upload1\Network performance with different numbers of AUVs\vbpsi_AUVnumber\plot_data';

path2 = 'C:\Users\zxx\Desktop\upload1\Network performance with different numbers of AUVs\vbpsc_AUVnumber\plot_data';

path3 = 'C:\Users\zxx\Desktop\upload1\Network performance with different numbers of AUVs\random_AUVnumber\plot_data';

path4 = 'C:\Users\zxx\Desktop\upload1\Network performance with different numbers of AUVs\laccm_AUVnumber\plot_data';

%% vbps-I

vuvoi_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vuvoi_vbps_average'))));
vuvoi_collision_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vuvoi_collision_vbps_average'))));
vuvoi_plot_vbps_average  = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vuvoi_plot'))));

vuvoitotaldelay_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vuvoitotaldelay_vbps_average'))));
vuvoisenddelay_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vuvoisenddelay_vbps_average'))));

vuvoidiscardnodecount_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vuvoidiscardnodecount_vbps_average'))));
vuvoicollisioncount_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vuvoicollisioncount_vbps_average'))));
vuvoiauvcollisioncount_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vuvoiauvcollisioncount_vbps_average'))));

vuvoithroughput_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vuvoithroughput_vbps_average'))));
vuvoitotalthroughput_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vuvoitotalthroughput_vbps_average'))));

vmvoitotaldelay_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vmvoitotaldelay_vbps_average'))));
vmvoisenddelay_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vmvoisenddelay_vbps_average'))));

vmvoidiscardnodecount_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vmvoidiscardnodecount_vbps_average'))));
vmvoicollisioncount_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vmvoicollisioncount_vbps_average'))));
vmvoiauvcollisioncount_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vmvoiauvcollisioncount_vbps_average'))));

vmvoithroughput_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vmvoithroughput_vbps_average'))));
vmvoitotalthroughput_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vmvoitotalthroughput_vbps_average'))));

vmvoi_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vmvoi_vbps_average'))));
vmvoi_collision_vbps_average  =  cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vmvoi_collision_vbps_average'))));
vmvoi_plot_vbps_average  = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpsauvnodenumbers_vmvoi_plot'))));

%% vbps-c
vuvoi_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vuvoi_vbpsc_average'))));
vuvoi_collision_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vuvoi_collision_vbpsc_average'))));
vuvoi_plot_vbpsc_average  = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vuvoi_plot'))));


vuvoitotaldelay_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vuvoitotaldelay_vbpsc_average'))));
vuvoisenddelay_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vuvoisenddelay_vbpsc_average'))));

vuvoidiscardnodecount_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vuvoidiscardnodecount_vbpsc_average'))));
vuvoicollisioncount_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vuvoicollisioncount_vbpsc_average'))));
vuvoiauvcollisioncount_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vuvoiauvcollisioncount_vbpsc_average'))));

vuvoithroughput_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vuvoithroughput_vbpsc_average'))));
vuvoitotalthroughput_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vuvoitotalthroughput_vbpsc_average'))));

vmvoitotaldelay_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vmvoitotaldelay_vbpsc_average'))));
vmvoisenddelay_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vmvoisenddelay_vbpsc_average'))));

vmvoidiscardnodecount_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vmvoidiscardnodecount_vbpsc_average'))));
vmvoicollisioncount_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vmvoicollisioncount_vbpsc_average'))));
vmvoiauvcollisioncount_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vmvoiauvcollisioncount_vbpsc_average'))));

vmvoithroughput_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vmvoithroughput_vbpsc_average'))));
vmvoitotalthroughput_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vmvoitotalthroughput_vbpsc_average'))));

vmvoi_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vmvoi_vbpsc_average'))));
vmvoi_collision_vbpsc_average  =  cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vmvoi_collision_vbpsc_average'))));
vmvoi_plot_vbpsc_average  = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscauvnodenumbers_vmvoi_plot'))));


%% random
vutotaldelay_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vutotaldelay_random_average'))));
vusenddelay_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vusenddelay_random_average'))));

vudiscardnodecount_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vudiscardnodecount_random_average'))));
vucollisioncount_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vucollisioncount_random_average'))));
vuauvcollisioncount_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vuauvcollisioncount_random_average'))));

vuthroughput_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vuthroughput_random_average'))));
vutotalthroughput_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vutotalthroughput_random_average'))));

vuvoi_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vuvoi_random_average'))));
vuvoi_collision_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vuvoi_collision_random_average'))));
vuvoi_plot_random_average  = cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vuvoi_plot'))));

vmtotaldelay_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vmtotaldelay_random_average'))));
vmsenddelay_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vmsenddelay_random_average'))));

vmdiscardnodecount_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vmdiscardnodecount_random_average'))));
vmcollisioncount_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vmcollisioncount_random_average'))));
vmauvcollisioncount_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vmauvcollisioncount_random_average'))));

vmthroughput_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vmthroughput_random_average'))));
vmtotalthroughput_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vmtotalthroughput_random_average'))));

vmvoi_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vmvoi_random_average'))));
vmvoi_collision_random_average  =  cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vmvoi_collision_random_average'))));
vmvoi_plot_random_average  = cell2mat(struct2cell(load(fullfile(path3,'20211220randomauvnodenumbers_vmvoi_plot'))));


%% laccm
vutotaldelay_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vutotaldelay_laccm_average'))));
vusenddelay_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vusenddelay_laccm_average'))));

vudiscardnodecount_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vudiscardnodecount_laccm_average'))));
vucollisioncount_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vucollisioncount_laccm_average'))));
vuauvcollisioncount_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vuauvcollisioncount_laccm_average'))));

vuthroughput_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vuthroughput_laccm_average'))));
vutotalthroughput_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vutotalthroughput_laccm_average'))));

vuvoi_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vuvoi_laccm_average'))));
vuvoi_collision_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vuvoi_collision_laccm_average'))));
vuvoi_plot_laccm_average  = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vuvoi_plot'))));

vmtotaldelay_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vmtotaldelay_laccm_average'))));
vmsenddelay_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vmsenddelay_laccm_average'))));

vmdiscardnodecount_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vmdiscardnodecount_laccm_average'))));
vmcollisioncount_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vmcollisioncount_laccm_average'))));
vmauvcollisioncount_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vmauvcollisioncount_laccm_average'))));

vmthroughput_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vmthroughput_laccm_average'))));
vmtotalthroughput_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vmtotalthroughput_laccm_average'))));

vmvoi_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vmvoi_laccm_average'))));
vmvoi_collision_laccm_average  =  cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vmvoi_collision_laccm_average'))));
vmvoi_plot_laccm_average  = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmauvnodenumbers_vmvoi_plot'))));


%% Plot
clr(1,:) = [0 113 188]/255;
clr(2,:) = [216 82 24]/255;
clr(3,:) = [236 176 31]/255;
clr(4,:) = [125 46 141]/255;
clr(5,:) = [118 171 47]/255;
clr(6,:) = [76 189 237]/255;
clr(7,:) = [161 19 46]/255;
clr(8,:) = [255 255 255]/255;
clr(9,:) = [25 35 45]/255;
clr(10,:) = [100 100 100]/255;
clr(11,:) = [150 150 150]/255;
clr(12,:) = [0 0 0]/255;




%% Collision Probability
figure;
 

for iNA_number = 1:5%1:length(vmtrafficload)
    
    A5(iNA_number) = sum(vmvoiauvcollisioncount_vbps_average(iNA_number,:));
end
plot(1:5,A5,'*-b','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));
hold on

 

for iNA_number = 1:5%1:length(vmtrafficload)
    
    A6(iNA_number) = sum(vmvoiauvcollisioncount_vbpsc_average(iNA_number,:));
end
plot(1:5,A6,'o-b','LineWidth',1.8,'MarkerSize',8,'color',clr(2,:));
hold on




for iNA_number = 1:5%1:length(vmtrafficload)
    
    A7(iNA_number) = sum(vmauvcollisioncount_random_average(iNA_number,:));
end
plot(1:5,A7,'x-b','LineWidth',1.8,'MarkerSize',8,'color',clr(3,:));
hold on


for iNA_number = 1:5%1:length(vmtrafficload)
    
    A8(iNA_number) = sum(vmauvcollisioncount_laccm_average(iNA_number,:));
end
plot(1:5,A8,'s-b','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
hold on


backColor  =  [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',1:1:5);
xlabel('\fontsize{16}Numbers of AUVs ');
ylabel('\fontsize{16}Collision Probability among AUVs');
legend1 = legend('\fontsize{16}\fontname{Times New Roman}VBPS-I','\fontsize{16}\fontname{Times New Roman}VBPS-C','\fontsize{16}\fontname{Times New Roman}RA','\fontsize{16}\fontname{Times New Roman}LACCM','Location','NorthEast');
set(legend1,...
    'Position',[0.142261904761905 0.661507968423228 0.226785714285714 0.248015873015873],'Color',[1 1 1]);
 
%% VoI
figure1 = figure;
backColor  =  [245 249 253]/255;
clr(1,:) = [0 113 188]/255;
clr(2,:) = [216 82 24]/255;
clr(3,:) = [236 176 31]/255;
clr(4,:) = [125 46 141]/255;
clr(5,:) = [118 171 47]/255;
clr(6,:) = [76 189 237]/255;
clr(7,:) = [161 19 46]/255;
clr(8,:) = [255 255 255]/255;
clr(9,:) = [25 35 45]/255;
clr(10,:) = [100 100 100]/255;
clr(11,:) = [150 150 150]/255;
clr(12,:) = [0 0 0]/255;


clr_b(1,:) = [0 113 188]/255;
clr_b(2,:) = [77 141 188]/255;
clr_b(3,:) = [114 157 188]/255;
clr_b(4,:) = [143 171 188]/255;
clr_b(5,:) = [167 186 188]/255;


clr_b(6,:) = [216 82 24]/255;
clr_b(7,:) = [216 118 76]/255;
clr_b(8,:) = [216 147 115]/255;
clr_b(9,:) = [216 165 146]/255;
clr_b(10,:) = [216 191 180]/255;

clr_b(11,:) = [236 176 31]/255;
clr_b(12,:) = [236 186 79]/255;
clr_b(13,:) = [236 199 117]/255;
clr_b(14,:) = [236 215 156]/255;
clr_b(15,:) = [236 227 204]/255;

clr_b(16,:) = [125 46 141]/255;
clr_b(17,:) = [132 67 145]/255;
clr_b(18,:) = [140 87 150]/255;
clr_b(19,:) = [155 129 161]/255;
clr_b(20,:) = [176 156 181]/255;


for iNA_number = 1:5
   vuvoi_vbps(iNA_number,:) = vuvoi_collision_vbps_average(iNA_number,:)+vmvoi_collision_vbps_average(iNA_number,:);
   vuvoi_vbpsc(iNA_number,:) = vuvoi_collision_vbpsc_average(iNA_number,:)+vmvoi_collision_vbpsc_average(iNA_number,:);
   vuvoi_random(iNA_number,:) = vuvoi_collision_random_average(iNA_number,:)+vmvoi_collision_random_average(iNA_number,:);
   vuvoi_laccm(iNA_number,:) = vuvoi_collision_laccm_average(iNA_number,:)+vmvoi_collision_laccm_average(iNA_number,:);
end


for iNA_number = 1:5
    for i = 1:5
        if i == 1
            b = bar((iNA_number-1)*5+i,vuvoi_vbps(iNA_number,:),'stacked','EdgeColor',clr(12,:));
            set(b(1),'facecolor',clr_b(1,:));
            set(b(2),'facecolor',clr_b(2,:));
            set(b(3),'facecolor',clr_b(3,:));
            set(b(4),'facecolor',clr_b(4,:));
            set(b(5),'facecolor',clr_b(5,:));
            hold on
        end
        if i == 2
            c = bar((iNA_number-1)*5+i,vuvoi_vbpsc(iNA_number,:),'stacked','EdgeColor',clr(12,:));
            set(c(1),'facecolor',clr_b(6,:));
            set(c(2),'facecolor',clr_b(7,:));
            set(c(3),'facecolor',clr_b(8,:));
            set(c(4),'facecolor',clr_b(9,:));
            set(c(5),'facecolor',clr_b(10,:));
            hold on
        end
        if i == 3
            d = bar((iNA_number-1)*5+i,vuvoi_random(iNA_number,:),'stacked','EdgeColor',clr(12,:));
            set(d(1),'facecolor',clr_b(11,:));
            set(d(2),'facecolor',clr_b(12,:));
            set(d(3),'facecolor',clr_b(13,:));
            set(d(4),'facecolor',clr_b(14,:));
            set(d(5),'facecolor',clr_b(15,:));
            hold on
        end
        if i == 4
           e = bar((iNA_number-1)*5+i,vuvoi_laccm(iNA_number,:),'stacked','EdgeColor',clr(12,:));
            set(e(1),'facecolor',clr_b(16,:));
            set(e(2),'facecolor',clr_b(17,:));
            set(e(3),'facecolor',clr_b(18,:));
            set(e(4),'facecolor',clr_b(19,:));
            set(e(5),'facecolor',clr_b(20,:));
            hold on
        end
        if i == 5
            f = bar((iNA_number-1)*5+i,[0,0,0,0,0],'stacked','EdgeColor',clr(12,:));
            set(f(1),'facecolor',clr_b(16,:));
            set(f(2),'facecolor',clr_b(17,:));
            set(f(3),'facecolor',clr_b(18,:));
            set(f(4),'facecolor',clr_b(19,:));
            set(f(5),'facecolor',clr_b(20,:));
            hold on
            
        end
       
    end
     
end

backColor  =  [255 255 255]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
  set(gca,'xtick',5:5:25);
set(gca,'ytick',0:1000:7000);
set(gca,'Color','none',...,
    'XTickLabel',...
    {' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '});
xlabel('\fontsize{16} Numbers of AUVs');
ylabel('\fontsize{16} Cumulative VoI');



axes2  =  axes('Parent',figure1);
hold(axes2,'on');
voi = zeros(5,5);
voi(1,5) =  0.1;
b = bar(1:5,voi);
ch  =  get(b,'children');
set(ch{1},'facecolor',clr(8,:));
set(ch{2},'facecolor',clr(8,:));
set(ch{3},'facecolor',clr(8,:));
set(ch{4},'facecolor',clr(8,:));
 set(ch{4},'facecolor',backColor,'EdgeColor',backColor);



backColor  =  [245 249 253]/255;
set(gca,'Color','none');

set(gca,'xtick',1:5);
 set(gca,'ytick',[]);set(gca,'ytick',0:1000:7000);
legend12 = legend([ch{1} ch{2} ch{3} ch{4}],'\fontsize{16}\fontname{Times New Roman}VBPS-I','\fontsize{16}\fontname{Times New Roman}VBPS-C','\fontsize{16}\fontname{Times New Roman}RA','\fontsize{16}\fontname{Times New Roman}LACCM','Location','NorthEast');
set(legend12, 'color', [1,1,1]);
 %%%Note: To avoid displaying extraneous data in the figure, it may be desirable to alter the blue color utilized in constructing the legend and x-y coordinates to a white color.



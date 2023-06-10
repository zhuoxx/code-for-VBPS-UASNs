clc;clear;close all

path1 = 'C:\Users\zxx\Desktop\upload1\Collision probability with vary packet sizes\vbpsi_packetlength\plot_data';

path2 = 'C:\Users\zxx\Desktop\upload1\Collision probability with vary packet sizes\vbpsc_packetlength\plot_data';

path3 = 'C:\Users\zxx\Desktop\upload1\Collision probability with vary packet sizes\random_packetlength\plot_data';

path4 = 'C:\Users\zxx\Desktop\upload1\Collision probability with vary packet sizes\laccm_packetlength\plot_data';
%% vbps
 
vuvoi_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vuvoi_vbps_average'))));
vuvoi_collision_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vuvoi_collision_vbps_average'))));


vuvoitotaldelay_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vuvoitotaldelay_vbps_average'))));
vuvoisenddelay_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vuvoisenddelay_vbps_average'))));

vuvoidiscardnodecount_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vuvoidiscardnodecount_vbps_average'))));
vuvoicollisioncount_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vuvoicollisioncount_vbps_average'))));
vuvoiauvcollisioncount_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vuvoiauvcollisioncount_vbps_average'))));

vuvoithroughput_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vuvoithroughput_vbps_average'))));
vuvoitotalthroughput_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vuvoitotalthroughput_vbps_average'))));

vmvoitotaldelay_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vmvoitotaldelay_vbps_average'))));
vmvoisenddelay_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vmvoisenddelay_vbps_average'))));

vmvoidiscardnodecount_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vmvoidiscardnodecount_vbps_average'))));
vmvoicollisioncount_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vmvoicollisioncount_vbps_average'))));
vmvoiauvcollisioncount_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vmvoiauvcollisioncount_vbps_average'))));

vmvoithroughput_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vmvoithroughput_vbps_average'))));
vmvoitotalthroughput_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vmvoitotalthroughput_vbps_average'))));

vmvoi_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vmvoi_vbps_average'))));
vmvoi_collision_vbps_average = cell2mat(struct2cell(load(fullfile(path1,'20211220vbpspacketlength_vmvoi_collision_vbps_average'))));

 

%% vbpsc
 
 
vuvoi_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vuvoi_vbpsc_average'))));
vuvoi_collision_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vuvoi_collision_vbpsc_average'))));


vuvoitotaldelay_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vuvoitotaldelay_vbpsc_average'))));
vuvoisenddelay_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vuvoisenddelay_vbpsc_average'))));

vuvoidiscardnodecount_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vuvoidiscardnodecount_vbpsc_average'))));
vuvoicollisioncount_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vuvoicollisioncount_vbpsc_average'))));
vuvoiauvcollisioncount_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vuvoiauvcollisioncount_vbpsc_average'))));

vuvoithroughput_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vuvoithroughput_vbpsc_average'))));
vuvoitotalthroughput_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vuvoitotalthroughput_vbpsc_average'))));

vmvoitotaldelay_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vmvoitotaldelay_vbpsc_average'))));
vmvoisenddelay_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vmvoisenddelay_vbpsc_average'))));

vmvoidiscardnodecount_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vmvoidiscardnodecount_vbpsc_average'))));
vmvoicollisioncount_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vmvoicollisioncount_vbpsc_average'))));
vmvoiauvcollisioncount_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vmvoiauvcollisioncount_vbpsc_average'))));

vmvoithroughput_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vmvoithroughput_vbpsc_average'))));
vmvoitotalthroughput_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vmvoitotalthroughput_vbpsc_average'))));

vmvoi_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vmvoi_vbpsc_average'))));
vmvoi_collision_vbpsc_average = cell2mat(struct2cell(load(fullfile(path2,'20211220vbpscpacketlength_vmvoi_collision_vbpsc_average'))));

 
%% random
%queue
vutotaldelay_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vutotaldelay_random_average'))));
vusenddelay_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vusenddelay_random_average'))));

vudiscardnodecount_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vudiscardnodecount_random_average'))));
vucollisioncount_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vucollisioncount_random_average'))));
vuauvcollisioncount_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vuauvcollisioncount_random_average'))));

vuthroughput_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vuthroughput_random_average'))));
vutotalthroughput_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vutotalthroughput_random_average'))));

vuvoi_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vuvoi_random_average'))));
vuvoi_collision_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vuvoi_collision_random_average'))));

vmtotaldelay_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vmtotaldelay_random_average'))));
vmsenddelay_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vmsenddelay_random_average'))));

vmdiscardnodecount_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vmdiscardnodecount_random_average'))));
vmcollisioncount_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vmcollisioncount_random_average'))));
vmauvcollisioncount_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vmauvcollisioncount_random_average'))));

vmthroughput_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vmthroughput_random_average'))));
vmtotalthroughput_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vmtotalthroughput_random_average'))));

vmvoi_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vmvoi_random_average'))));
vmvoi_collision_random_average = cell2mat(struct2cell(load(fullfile(path3,'20211220randompacketlength_vmvoi_collision_random_average'))));

 %% laccm
vutotaldelay_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vutotaldelay_laccm_average'))));
vusenddelay_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vusenddelay_laccm_average'))));

vudiscardnodecount_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vudiscardnodecount_laccm_average'))));
vucollisioncount_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vucollisioncount_laccm_average'))));
vuauvcollisioncount_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vuauvcollisioncount_laccm_average'))));

vuthroughput_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vuthroughput_laccm_average'))));
vutotalthroughput_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vutotalthroughput_laccm_average'))));

vuvoi_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vuvoi_laccm_average'))));
vuvoi_collision_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vuvoi_collision_laccm_average'))));

vmtotaldelay_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vmtotaldelay_laccm_average'))));
vmsenddelay_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vmsenddelay_laccm_average'))));

vmdiscardnodecount_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vmdiscardnodecount_laccm_average'))));
vmcollisioncount_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vmcollisioncount_laccm_average'))));
vmauvcollisioncount_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vmauvcollisioncount_laccm_average'))));

vmthroughput_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vmthroughput_laccm_average'))));
vmtotalthroughput_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vmtotalthroughput_laccm_average'))));

vmvoi_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vmvoi_laccm_average'))));
vmvoi_collision_laccm_average = cell2mat(struct2cell(load(fullfile(path4,'20211220laccmpacketlength_vmvoi_collision_laccm_average'))));

 
%%
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

L_packetrange = 500:500:5000;

 
  
 
 

for iL_packet = 1:length(L_packetrange)%1:length(vmtrafficload)
    A1(iL_packet) = sum(vuvoicollisioncount_vbps_average(iL_packet,:));
    
end

for iL_packet = 1:length(L_packetrange)%1:length(vmtrafficload)
    A2(iL_packet) = sum(vuvoicollisioncount_vbpsc_average(iL_packet,:));
    
end

for iL_packet = 1:length(L_packetrange)%1:length(vmtrafficload)
    A3(iL_packet) = sum(vucollisioncount_random_average(iL_packet,:));
    
end

for iL_packet = 1:length(L_packetrange)%1:length(vmtrafficload)
    A4(iL_packet) = sum(vucollisioncount_laccm_average(iL_packet,:));
    
end

for iL_packet = 1:length(L_packetrange)%1:length(vmtrafficload)
    
    A5(iL_packet) = sum(vmvoicollisioncount_vbps_average(iL_packet,:));
end


 
for iL_packet = 1:length(L_packetrange)%1:length(vmtrafficload)
    
    A6(iL_packet) = sum(vmvoicollisioncount_vbpsc_average(iL_packet,:));
end


 

for iL_packet = 1:length(L_packetrange)%1:length(vmtrafficload)
    
    A7(iL_packet) = sum(vmcollisioncount_random_average(iL_packet,:));
end



for iL_packet = 1:length(L_packetrange)%1:length(vmtrafficload)
    
    A8(iL_packet) = sum(vmcollisioncount_laccm_average(iL_packet,:));
end


 

figure1 = figure;
collisionrate = [A1+A5;A2+A6;A3+A7;A4+A8];
b = bar(500:500:5000,collisionrate);
ch  =  get(b,'children');
set(ch{1},'facecolor',clr(1,:));
set(ch{2},'facecolor',clr(6,:));
set(ch{3},'facecolor',clr(3,:));
set(ch{4},'facecolor',clr(4,:));
 

 
plot(500:500:5000,A1+A5,'*-b','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));
hold on
plot(500:500:5000,A2+A6,'o-b','LineWidth',1.8,'MarkerSize',8,'color',clr(2,:));
hold on
plot(500:500:5000,A3+A7,'x-b','LineWidth',1.8,'MarkerSize',8,'color',clr(3,:));
hold on
plot(500:500:5000,A4+A8,'s-b','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));

backColor  =  [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',500:500:5000);
xlabel('\fontsize{16} Packet Sizes (Bits) ');
ylabel('\fontsize{16} Collision Probability for Static nodes');
 
 h = legend([ch{1} ch{2} ch{3} ch{4}],'\fontsize{16}\fontname{Times New Roman}VBPS-I','\fontsize{16}\fontname{Times New Roman}VBPS-C','\fontsize{16}\fontname{Times New Roman}RA','\fontsize{16}\fontname{Times New Roman}LACCM','Location','NorthEast');

 set(h, 'color', [1,1,1]);

position = [0.45 0.2 0.42 0.2];
axes3  =  axes('Parent',figure1,...
    'Position',position);
hold(axes3,'on');

plot(500:500:5000,A1+A5,'*-b','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));
hold on
plot(500:500:5000,A2+A6,'o-b','LineWidth',1.8,'MarkerSize',8,'color',clr(2,:));
hold on
plot(500:500:5000,A4+A8,'s-b','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
set(gca,'xtick',500:500:5000);
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([500 ,5000,0,0.02])
 
annotation(figure1,'arrow',[0.901785714285714 0.876785714285714],...
    [0.134714285714286 0.221428571428571]);

 
annotation(figure1,'arrow',[0.1375 0.432142857142857],...
    [0.125190476190476 0.195238095238095]);

 

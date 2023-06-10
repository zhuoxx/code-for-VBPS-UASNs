%%%%%%%%%%%%%%%%%%%%%% Main file for VBPS-I with different numbers of AUVs

clc;clear;close all;

%%-------------------------Parameters----------------------------------------
dc  =  3500;%m
f  =  26000;%Hz
R  =  13900;%bps
L_packet  =  1000;%bit
t_data  =  L_packet/R;
vs  =  1500;%m/s
n  =  ceil((dc/vs+t_data)/t_data);

v_A  =  10;%m/s

iteration_total = 1000;
NArange = 1:5;

%%%%%%%%%%%%%%%%%%%%%%%%%simulation network in this paper, it is better not
%%%%%%%%%%%%%%%%%%%%%%%%%to modifiy this part%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate locations for predeployed static nodes
N0 = 3;        %Number of sink nodes
N1 = N0*5;     %Number of relay nodes
N2 = N0*5*3;   %Number of sensor nodes
N = N0+N1+N2;  %Total static node number

path1 = 'C:\Users\zxx\Desktop\upload1\Network performance with different numbers of AUVs\vbpsi_AUVnumber\networkdata'; %networkdata file path
path2 = 'C:\Users\zxx\Desktop\upload1\Network performance with different numbers of AUVs\vbpsi_AUVnumber\plot_data'; %plotdata file path

 
[Nx,Ny,Nz]  =  generate_staticnode_location(N0,N1,N2,dc,path1); %Generate locations for predeployed static nodes
 

%% Generate the network topology and sending time for each static nodes based on TDMA protocol
Round  =  50;  %Communication round
max_pathcount  =  120;
totaltimeslot  =  Round*max_pathcount;
[delay_table,dis_table,TM,IM,timeslot,sequence,node,path,Ts]  =  generate_TDMA(N,Nx,Ny,Nz,dc,vs,N0,N1,N2,path1,n,t_data,Round,totaltimeslot);

%%%%%%%%%%%%%%%%%%%%%%%%%simulation network in this paper, it is better not
%%%%%%%%%%%%%%%%%%%%%%%%%to modifiy the code in this part%%%%%%%%%%%%%%%%%%



for iteration = 1:iteration_total
    for iNA_number = 1:5
        
        %% Path Preplanning for AUVs
        NA = iNA_number;
         
        [Ax,Ay,Az, len_xyz,timeplus]  =  generate_AUVpath(Nx,Ny,Nz, NA,path1);
      
        
        %% Static node localization
        
        [AUV_listen,AUV_detected,AUV_in]  =   Static_node_localization(Ax, Ay, Az,Ts, Nx, Ny, Nz,n,t_data,path1,N0,N1,N2,NA,N,Round,totaltimeslot,timeplus,path,node,len_xyz,dc);
        
        
        %% The network throughput of static nodes
        timestart = 10;
        timeend = 20;
        
        totaltime = (len_xyz.*timeplus-timestart-timeend);
        slotnumber_static = floor(totaltime./n./t_data);
        
        sendpacketnumber_static = slotnumber_static.*(N0+N1);
        
        
        %% Transmission Scheduling
        
        vutrafficload = 0.05;
        vmtrafficload = 0.5;
        betau = 0.5;
        vu0 = L_packet;
        alphau = 3;
        Tu = 5*n*t_data;
        betam = 0.3;
        vm0 = L_packet;
        alpham = 3;
        Tm = 10*n*t_data;
        
        
        
        for ivutrafficload = 1:length(vutrafficload)  %pkt/s
            
            vulambda = vutrafficload(ivutrafficload);
            
            [vuAUVrun_detected,vuarrivaltime,vuAUVrun_in,vuAUVrun_knownnodes,vuAUVrun_knownnodesin,vuAUVrun_unknownnodes,vuAUVrun_unknownnodesin,vuAUVtimebefore,vuAUVtimeplus,vutimeslot_number,vutimeslot_numberplus,vudiscardnodecount,vuknown_nodescount,vuknown_nodesincount,vuqueuearrivaltime,vuqueuecount,vuqueueindex,vuqueueout,vuqueuesendtime,vusendauvrun_time,vusendauvrun_timecase1,vusendauvrun_timecase2,vusendauvrun_timecase3,vusendauvrun_timecase4,vusenddelay,vutime_interval,vutotaldelay,vuunknown_nodescount, vuunknown_nodesincount,vuAUVrun_distancenode,vucollisioncount,vutotaltimearrival,vuvoi_choose,vuvoi_chooseindex,vuvoi_choosenode,vuvoitotaldelay_choose,vuvoicollisioncount,vuvoi,vuauvcollisioncount,vuvoisenddelay_choose,vuvoiauvcollisioncount_choose,vuvoi_collision]  =  VBPSIalgorithm(Nx,Ny,Nz,node,path,timeslot,len_xyz,dis_table,Ax,Ay,Az,AUV_in,AUV_detected,vulambda,NA,N,timeplus,n,t_data,N0,N1,N2,vs,vutrafficload,ivutrafficload,iteration,betau,vu0,alphau,Tu,L_packet,timestart,timeend);
            for ivmtrafficload = 1:length(vmtrafficload)
                vmlambda = vmtrafficload(ivmtrafficload); %pkt/s
                [vmAUVrun_detected,vmarrivaltime,vmAUVrun_in,vmAUVrun_knownnodes,vmAUVrun_knownnodesin,vmAUVrun_unknownnodes,vmAUVrun_unknownnodesin,vmAUVtimebefore,vmAUVtimeplus,vmtimeslot_number,vmtimeslot_numberplus,vmdiscardnodecount,vmknown_nodescount,vmknown_nodesincount,vmqueuearrivaltime,vmqueuecount,vmqueueindex,vmqueueout,vmqueuesendtime,vmsendauvrun_time,vmsendauvrun_timecase1,vmsendauvrun_timecase2,vmsendauvrun_timecase3,vmsendauvrun_timecase4,vmsenddelay,vmtime_interval,vmtotaldelay,vmtotaldelay_choose, vmunknown_nodescount,vmAUVrun_distancenode,vmcollisioncount,vmtotaltimearrival,vmvoi_choose,vmvoi_chooseindex,vmvoi_choosenode,vmvoitotaldelay_choose,vmvoicollisioncount,vmvoi,vmauvcollisioncount,vmvoisenddelay_choose,vmvoiauvcollisioncount_choose,vmvoi_collision]  =  VBPSIalgorithm(Nx,Ny,Nz,node,path,timeslot,len_xyz,dis_table,Ax,Ay,Az,AUV_in,AUV_detected,vmlambda,NA,N,timeplus,n,t_data,N0,N1,N2,vs,vmtrafficload,ivmtrafficload,iteration,betam,vm0,alpham,Tm,L_packet,timestart,timeend);
                
                
                %% Calculate VoI, network throughput, congestion ratio, and  network delay
                
                vuvoidiscardnodecount = vudiscardnodecount;
                
                for iNA_plot = 1:NA
                    vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0; % for vu
                    vuvoitotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vuvoitotaldelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vuvoisenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vuvoisenddelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    
                    vuvoidiscardnodecount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vuvoicollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vuvoiauvcollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vuvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vuvoi_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vuvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vuvoi_collision_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vuvoithroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vuvoitotalthroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    
                    
                    
                    
                    for ivuarrival = 1:vutotaltimearrival(iNA_plot,ivutrafficload)
                        if vuvoitotaldelay_choose(iNA_plot,ivuarrival)<1000
                            vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                            vuvoitotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoitotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoitotaldelay_choose(iNA_plot,ivuarrival);
                            vuvoisenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoisenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoisenddelay_choose(iNA_plot,ivuarrival);
                            
                            vuvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoi_collision(iNA_plot,ivuarrival);
                            vuvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoi_choose(iNA_plot,ivuarrival);
                            
                        end
                    end
                    if vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)~= 0
                        vuvoitotaldelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoitotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        vuvoisenddelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoisenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        
                        vuvoi_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        vuvoi_collision_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        
                        vuvoidiscardnodecount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vuvoidiscardnodecount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload);
                        vuvoicollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vuvoicollisioncount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload);
                        vuvoiauvcollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vuvoiauvcollisioncount_choose(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload);
                        
                        vuvoithroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = (vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)-sum(vuvoicollisioncount(iNA_plot,:))-sum(vuvoiauvcollisioncount_choose(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);
                        vuvoitotalthroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = (sendpacketnumber_static(iNA_plot)-sum(vuvoicollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);
                        
                        
                    end
                end
                
                
                vmvoidiscardnodecount = vmdiscardnodecount;
                
                for iNA_plot = 1:NA
                    vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0; % for vu
                    vmvoitotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vmvoitotaldelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vmvoisenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vmvoisenddelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    
                    vmvoidiscardnodecount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vmvoicollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vmvoiauvcollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vmvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vmvoi_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vmvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vmvoi_collision_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vmvoithroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    vmvoitotalthroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                    
                    
                    
                    for ivmarrival = 1:vmtotaltimearrival(iNA_plot,ivmtrafficload)
                        if vmvoitotaldelay_choose(iNA_plot,ivmarrival)<1000
                            vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                            vmvoitotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoitotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoitotaldelay_choose(iNA_plot,ivmarrival);
                            vmvoisenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoisenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoitotaldelay_choose(iNA_plot,ivmarrival);
                            
                            vmvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoi_choose(iNA_plot,ivmarrival);
                            vmvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoi_collision(iNA_plot,ivmarrival);
                            
                        end
                    end
                    if vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)~= 0
                        vmvoitotaldelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoitotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        vmvoisenddelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoisenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        
                        
                        vmvoidiscardnodecount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vmvoidiscardnodecount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload);
                        vmvoicollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vmvoicollisioncount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload);
                        vmvoiauvcollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vmvoiauvcollisioncount_choose(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload);
                        vmvoithroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = (vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)-sum(vmvoicollisioncount(iNA_plot,:))-sum(vmvoiauvcollisioncount_choose(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);
                        vmvoitotalthroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = (sendpacketnumber_static(iNA_plot)-sum(vmvoicollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);
                        
                        vmvoi_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        vmvoi_collision_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        
                        
                    end
                end
                
            end
            
        end
    end
    
%     %     % save data
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoi_choose'),'vuvoi_choose');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoi_collision'),'vuvoi_collision');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoi_plot'),'vuvoi_plot');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoi_collision_plot'),'vuvoi_collision_plot');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoi_plot_sum'),'vuvoi_plot_sum');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoi_collision_plot_sum'),'vuvoi_collision_plot_sum');
%     
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoitotaldelay_choose'),'vuvoitotaldelay_choose');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoitotaldelay_plot'),'vuvoitotaldelay_plot');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoisenddelay_choose'),'vuvoisenddelay_choose');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoisenddelay_plot'),'vuvoisenddelay_plot');
%     
%     
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoidiscardnodecount_plot'),'vuvoidiscardnodecount_plot');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoicollisioncount_plot'),'vuvoicollisioncount_plot');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoiauvcollisioncount_plot'),'vuvoiauvcollisioncount_plot');
%     
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoithroughput_plot'),'vuvoithroughput_plot');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoitotalthroughput_plot'),'vuvoitotalthroughput_plot');
%     
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoiarrivaltime'),'vuvoiarrivaltime');
%     
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoitotaldelay_plot_sum'),'vuvoitotaldelay_plot_sum');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoisenddelay_plot_sum'),'vuvoisenddelay_plot_sum');
%     
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoi_choose'),'vmvoi_choose');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoi_collision'),'vmvoi_collision');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoi_plot'),'vmvoi_plot');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoi_collision_plot'),'vmvoi_collision_plot');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoi_plot_sum'),'vmvoi_plot_sum');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoi_collision_plot_sum'),'vmvoi_collision_plot_sum');
%     
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoitotaldelay_choose'),'vmvoitotaldelay_choose');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoitotaldelay_plot'),'vmvoitotaldelay_plot');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoisenddelay_choose'),'vmvoisenddelay_choose');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoisenddelay_plot'),'vmvoisenddelay_plot');
%     
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoidiscardnodecount_plot'),'vmvoidiscardnodecount_plot');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoicollisioncount_plot'),'vmvoicollisioncount_plot');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoiauvcollisioncount_plot'),'vmvoiauvcollisioncount_plot');
%     
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoithroughput_plot'),'vmvoithroughput_plot');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoitotalthroughput_plot'),'vmvoitotalthroughput_plot');
%     
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoiarrivaltime'),'vmvoiarrivaltime');
%     
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoitotaldelay_plot_sum'),'vmvoitotaldelay_plot_sum');
%         save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoisenddelay_plot_sum'),'vmvoisenddelay_plot_sum');
     iteration
end





%% Average
for iNA_number = 1:5
    for iNA_plot = 1:NA
        for ivutrafficload = 1:length(vutrafficload)
            for ivmtrafficload = 1:length(vmtrafficload)
                vuvoitotaldelay_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoitotaldelay_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoitotaldelay_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoitotaldelay_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoisenddelay_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoisenddelay_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoisenddelay_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoisenddelay_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                
                
                vuvoidiscardnodecount_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoidiscardnodecount_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoidiscardnodecount_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoidiscardnodecount_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoicollisioncount_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoicollisioncount_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoicollisioncount_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoicollisioncount_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoiauvcollisioncount_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoiauvcollisioncount_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoiauvcollisioncount_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoiauvcollisioncount_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoithroughput_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoithroughput_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoithroughput_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoithroughput_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoi_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoi_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoi_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoi_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoi_collision_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoi_collision_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoi_collision_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoi_collision_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoitotalthroughput_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoitotalthroughput_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoitotalthroughput_vbps_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoitotalthroughput_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            end
        end
    end
end




% % %%save data
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoi_vbps_average'),'vuvoi_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoi_collision_vbps_average'),'vuvoi_collision_vbps_average');
% 
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoitotaldelay_vbps_average'),'vuvoitotaldelay_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoisenddelay_vbps_average'),'vuvoisenddelay_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoidiscardnodecount_vbps_average'),'vuvoidiscardnodecount_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoicollisioncount_vbps_average'),'vuvoicollisioncount_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoiauvcollisioncount_vbps_average'),'vuvoiauvcollisioncount_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoithroughput_vbps_average'),'vuvoithroughput_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vuvoitotalthroughput_vbps_average'),'vuvoitotalthroughput_vbps_average');
% 
% 
% 
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoi_vbps_average'),'vmvoi_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoi_collision_vbps_average'),'vmvoi_collision_vbps_average');
% 
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoitotaldelay_vbps_average'),'vmvoitotaldelay_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoisenddelay_vbps_average'),'vmvoisenddelay_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoidiscardnodecount_vbps_average'),'vmvoidiscardnodecount_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoicollisioncount_vbps_average'),'vmvoicollisioncount_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoiauvcollisioncount_vbps_average'),'vmvoiauvcollisioncount_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoithroughput_vbps_average'),'vmvoithroughput_vbps_average');
% save(fullfile(path2,'20211220vbpsauvnodenumbers_vmvoitotalthroughput_vbps_average'),'vmvoitotalthroughput_vbps_average');
% 
%%
clr(1,:) = [0 113 188]/255;
clr(2,:) = [216 82 24]/255;
clr(3,:) = [236 176 31]/255;
clr(4,:) = [125 46 141]/255;
clr(5,:) = [118 171 47]/255;
clr(6,:) = [76 189 237]/255;
clr(7,:) = [161 19 46]/255;
clr(8,:) = [0 0 0]/255;
clr(9,:) = [25 35 45]/255;
clr(10,:) = [100 100 100]/255;
clr(11,:) = [150 150 150]/255;

%collision probability
figure1 = figure;
for iNA_number = 1:5%1:length(vmtrafficload)
    
    A5(iNA_number) = sum(vmvoiauvcollisioncount_vbps_average(iNA_number,:));
end
plot(NArange,A5,'*-b','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));

backColor  =  [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',1:1:5);
xlabel('\fontsize{16}Numbers of AUVs ');
ylabel('\fontsize{16}Collision Probability among AUVs');



%%VoI
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
%%%Note: To avoid displaying extraneous data in the figure, it may be desirable to alter the blue color utilized in constructing the legend and x-y coordinates to a white color.








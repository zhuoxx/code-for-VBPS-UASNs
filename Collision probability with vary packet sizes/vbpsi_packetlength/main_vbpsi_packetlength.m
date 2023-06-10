%%%%%%%%%%%%%%%%%%%%%% Main file for VBPS-I with vary packet sizes

clc;clear;close all;

%%-------------------------Parameters----------------------------------------
dc  =  3500;%m
f  =  26000;%Hz
R  =  13900;%bps
vs  =  1500;%m/s
v_A  =  10;%m/s


%% Generate locations for predeployed static nodes
N0 = 3;
N1 = N0*5;
N2 = N0*5*3;
N = N0+N1+N2;


path1 = 'C:\Users\zxx\Desktop\upload\Collision probability with vary packet sizes\laccm_packetlength\networkdata'; %networkdata file path
path2 = 'C:\Users\zxx\Desktop\upload\Collision probability with vary packet sizes\laccm_packetlength\plot_data'; %plotdata file path

tic
[Nx,Ny,Nz]  =  generate_staticnode_location(N0,N1,N2,dc,path1); %Generate locations for predeployed static nodes
toc

 
%% Path Preplanning for AUVs
NA = 2;
tic
[Ax,Ay,Az, len_xyz,timeplus]  =  generate_AUVpath(Nx,Ny,Nz, NA,path1);
toc
%%
iteration_total  = 2000;
N0_length  =  1;
N0range  =  3;
L_packetrange  =  500:500:5000;
for iteration  =  1:iteration_total
    for  iL_packet  =  1:length(L_packetrange)
        L_packet  =  L_packetrange(iL_packet);%bit
        t_data  =  L_packet/R;
        n  =  ceil((dc/vs+t_data)/t_data);
        
        
        %% Generate the network topology and sending time for each static nodes based on TDMA protocol
        Round   =   50;  %Communication round
        max_pathcount   =   120;
        totaltimeslot   =   Round*max_pathcount;
        [delay_table,dis_table,TM,IM,timeslot,sequence,node,path,Ts]   =   generate_TDMA(N,Nx,Ny,Nz,dc,vs,N0,N1,N2,path1,n,t_data,Round,totaltimeslot,iL_packet);
        
        %% Static node localization
        tic
        [AUV_listen,AUV_detected,AUV_in]   =    Static_node_localization(Ax, Ay, Az,Ts, Nx, Ny, Nz,n,t_data,path1,N0,N1,N2,NA,N,Round,totaltimeslot,timeplus,path,node,len_xyz,dc,iL_packet);
        toc
        
        %% The network throughput of static nodes
        timestart  =  10;
        timeend  =  20;
        
        totaltime  =  (len_xyz.*timeplus-timestart-timeend);
        slotnumber_static  =  floor(totaltime./n./t_data);
        
        sendpacketnumber_static  =  slotnumber_static.*(N0+N1);
        
        
        %% Transmission Scheduling
        
        vutrafficload  =  0.05;
        vmtrafficload  =  0.5;
        
        
        
        
        betau  =  0.5;
        vu0  =  L_packet;
        alphau  =  3;
        Tu  =  5*n*t_data;
        
        
        betam  =  0.3;
        vm0  =  L_packet;
        alpham  =  3;
        Tm  =  10*n*t_data;
        
        
        
        for ivutrafficload  =  1:length(vutrafficload)  %pkt/s
          
            vulambda  =  vutrafficload(ivutrafficload);
            [vuAUVrun_detected,vuarrivaltime,vuAUVrun_in,vuAUVrun_knownnodes,vuAUVrun_knownnodesin,vuAUVrun_unknownnodes,vuAUVrun_unknownnodesin,vuAUVtimebefore,vuAUVtimeplus,vutimeslot_number,vutimeslot_numberplus,vudiscardnodecount,vuknown_nodescount,vuknown_nodesincount,vuqueuearrivaltime,vuqueuecount,vuqueueindex,vuqueueout,vuqueuesendtime,vusendauvrun_time,vusendauvrun_timecase1,vusendauvrun_timecase2,vusendauvrun_timecase3,vusendauvrun_timecase4,vusenddelay,vutime_interval,vutotaldelay,vuunknown_nodescount, vuunknown_nodesincount,vuAUVrun_distancenode,vucollisioncount,vutotaltimearrival,vuvoi_choose,vuvoi_chooseindex,vuvoi_choosenode,vuvoitotaldelay_choose,vuvoicollisioncount,vuvoi,vuauvcollisioncount,vuvoisenddelay_choose,vuvoiauvcollisioncount_choose,vuvoi_collision]   =   VBPSIalgorithm(Nx,Ny,Nz,node,path,timeslot,len_xyz,dis_table,Ax,Ay,Az,AUV_in,AUV_detected,vulambda,NA,N,timeplus,n,t_data,N0,N1,N2,vs,vutrafficload,ivutrafficload,iteration,betau,vu0,alphau,Tu,L_packet,timestart,timeend);
            for ivmtrafficload  =  1:length(vmtrafficload)
                vmlambda  =  vmtrafficload(ivmtrafficload); %pkt/s
                [vmAUVrun_detected,vmarrivaltime,vmAUVrun_in,vmAUVrun_knownnodes,vmAUVrun_knownnodesin,vmAUVrun_unknownnodes,vmAUVrun_unknownnodesin,vmAUVtimebefore,vmAUVtimeplus,vmtimeslot_number,vmtimeslot_numberplus,vmdiscardnodecount,vmknown_nodescount,vmknown_nodesincount,vmqueuearrivaltime,vmqueuecount,vmqueueindex,vmqueueout,vmqueuesendtime,vmsendauvrun_time,vmsendauvrun_timecase1,vmsendauvrun_timecase2,vmsendauvrun_timecase3,vmsendauvrun_timecase4,vmsenddelay,vmtime_interval,vmtotaldelay,vmtotaldelay_choose, vmunknown_nodescount,vmAUVrun_distancenode,vmcollisioncount,vmtotaltimearrival,vmvoi_choose,vmvoi_chooseindex,vmvoi_choosenode,vmvoitotaldelay_choose,vmvoicollisioncount,vmvoi,vmauvcollisioncount,vmvoisenddelay_choose,vmvoiauvcollisioncount_choose,vmvoi_collision]   =   VBPSIalgorithm(Nx,Ny,Nz,node,path,timeslot,len_xyz,dis_table,Ax,Ay,Az,AUV_in,AUV_detected,vmlambda,NA,N,timeplus,n,t_data,N0,N1,N2,vs,vmtrafficload,ivmtrafficload,iteration,betam,vm0,alpham,Tm,L_packet,timestart,timeend);
                
                
                %% Calculate VoI, network throughput, congestion ratio, and  network delay
                
                vuvoidiscardnodecount  =  vudiscardnodecount;
                
                for iNA_plot  =  1:NA
                    vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0; % for vu
                    vuvoitotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vuvoitotaldelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vuvoisenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vuvoisenddelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    
                    vuvoidiscardnodecount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vuvoicollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vuvoiauvcollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vuvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vuvoi_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vuvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vuvoi_collision_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vuvoithroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;%%zxx1207
                    vuvoitotalthroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    
                    
                    
                    
                    for ivuarrival  =  1:vutotaltimearrival(iNA_plot,ivutrafficload)
                        if vuvoitotaldelay_choose(iNA_plot,ivuarrival)<1000
                            vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                            vuvoitotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vuvoitotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoitotaldelay_choose(iNA_plot,ivuarrival);
                            vuvoisenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vuvoisenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoisenddelay_choose(iNA_plot,ivuarrival);
                            
                            vuvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vuvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoi_collision(iNA_plot,ivuarrival);
                            vuvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vuvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoi_choose(iNA_plot,ivuarrival);
                            
                            
                        end
                    end
                    if vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)~=  0
                        vuvoitotaldelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vuvoitotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        vuvoisenddelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vuvoisenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        
                        vuvoi_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vuvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        vuvoi_collision_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vuvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        
                        vuvoidiscardnodecount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  sum(vuvoidiscardnodecount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload); %周围没有节点，也没法传输 rate
                        vuvoicollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  sum(vuvoicollisioncount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload);%由于不是globalmap导致冲突 rate
                        vuvoiauvcollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  sum(vuvoiauvcollisioncount_choose(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload);%由于不是globalmap导致冲突 rate
                        
                        vuvoithroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  (vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)-sum(vuvoicollisioncount(iNA_plot,:))-sum(vuvoiauvcollisioncount_choose(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);%吞吐量
                        vuvoitotalthroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  (sendpacketnumber_static(iNA_plot)-sum(vuvoicollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);%吞吐量
                        
                    end
                end
                
                
                vmvoidiscardnodecount  =  vmdiscardnodecount;
                
                for iNA_plot  =  1:NA
                    vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0; % for vu
                    vmvoitotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vmvoitotaldelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vmvoisenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vmvoisenddelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    
                    vmvoidiscardnodecount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vmvoicollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vmvoiauvcollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vmvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vmvoi_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vmvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vmvoi_collision_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    vmvoithroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;%%zxx1207
                    vmvoitotalthroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  0;
                    
                    
                    for ivmarrival  =  1:vmtotaltimearrival(iNA_plot,ivmtrafficload)
                        if vmvoitotaldelay_choose(iNA_plot,ivmarrival)<1000
                            vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                            vmvoitotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vmvoitotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoitotaldelay_choose(iNA_plot,ivmarrival);
                            vmvoisenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vmvoisenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoitotaldelay_choose(iNA_plot,ivmarrival);
                            
                            vmvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vmvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoi_choose(iNA_plot,ivmarrival);
                            vmvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vmvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoi_collision(iNA_plot,ivmarrival);
                            
                        end
                    end
                    if vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)~=  0
                        vmvoitotaldelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vmvoitotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        vmvoisenddelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  vmvoisenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
     
                        vmvoidiscardnodecount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)  =  sum(vmvoidiscardnodecount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload); %周围没有节点，也没法传输 rate zxx1209
                        vmvoicollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vmvoicollisioncount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload);%由于不是globalmap导致冲突 rate zxx1209
                        vmvoiauvcollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vmvoiauvcollisioncount_choose(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload);%由于不是globalmap导致冲突 rate zxx1209
                        vmvoithroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = (vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)-sum(vmvoicollisioncount(iNA_plot,:))-sum(vmvoiauvcollisioncount_choose(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);%吞吐量
                        vmvoitotalthroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = (sendpacketnumber_static(iNA_plot)-sum(vmvoicollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);%吞吐量
                        
                        vmvoi_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        vmvoi_collision_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                        
                        
                    end
                end
                
            end
        end
    end
    
    % % %%save data
    % %
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoi_choose'),'vuvoi_choose');
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoi_collision'),'vuvoi_collision');
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoi_plot'),'vuvoi_plot');
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoi_collision_plot'),'vuvoi_collision_plot');
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoi_plot_sum'),'vuvoi_plot_sum');
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoi_collision_plot_sum'),'vuvoi_collision_plot_sum');
    %
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoitotaldelay_choose'),'vuvoitotaldelay_choose');
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoitotaldelay_plot'),'vuvoitotaldelay_plot');
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoisenddelay_choose'),'vuvoisenddelay_choose');
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoisenddelay_plot'),'vuvoisenddelay_plot');
     %
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoidiscardnodecount_plot'),'vuvoidiscardnodecount_plot');
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoicollisioncount_plot'),'vuvoicollisioncount_plot');
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoiauvcollisioncount_plot'),'vuvoiauvcollisioncount_plot');
    %
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoithroughput_plot'),'vuvoithroughput_plot');
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoitotalthroughput_plot'),'vuvoitotalthroughput_plot');
    %
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoiarrivaltime'),'vuvoiarrivaltime');
    %
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoitotaldelay_plot_sum'),'vuvoitotaldelay_plot_sum');
    % save(fullfile(path2,'20211220vbpspacketlength_vuvoisenddelay_plot_sum'),'vuvoisenddelay_plot_sum');
     %
    
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoi_choose'),'vmvoi_choose');
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoi_collision'),'vmvoi_collision');
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoi_plot'),'vmvoi_plot');
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoi_collision_plot'),'vmvoi_collision_plot');
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoi_plot_sum'),'vmvoi_plot_sum');
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoi_collision_plot_sum'),'vmvoi_collision_plot_sum');
    %
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoitotaldelay_choose'),'vmvoitotaldelay_choose');
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoitotaldelay_plot'),'vmvoitotaldelay_plot');
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoisenddelay_choose'),'vmvoisenddelay_choose');
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoisenddelay_plot'),'vmvoisenddelay_plot');
     %
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoidiscardnodecount_plot'),'vmvoidiscardnodecount_plot');
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoicollisioncount_plot'),'vmvoicollisioncount_plot');
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoiauvcollisioncount_plot'),'vmvoiauvcollisioncount_plot');
    %
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoithroughput_plot'),'vmvoithroughput_plot');
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoitotalthroughput_plot'),'vmvoitotalthroughput_plot');
    %
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoiarrivaltime'),'vmvoiarrivaltime');
    %
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoitotaldelay_plot_sum'),'vmvoitotaldelay_plot_sum');
    % save(fullfile(path2,'20211220vbpspacketlength_vmvoisenddelay_plot_sum'),'vmvoisenddelay_plot_sum');
     iteration
end





%delay和 voi 取平均
for iL_packet = 1:length(L_packetrange)
    for iNA_plot = 1:NA
        for ivutrafficload = 1:length(vutrafficload)
            for ivmtrafficload = 1:length(vmtrafficload)
                vuvoitotaldelay_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoitotaldelay_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoitotaldelay_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoitotaldelay_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoisenddelay_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoisenddelay_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoisenddelay_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoisenddelay_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                
                
                vuvoidiscardnodecount_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoidiscardnodecount_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoidiscardnodecount_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoidiscardnodecount_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoicollisioncount_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoicollisioncount_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoicollisioncount_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoicollisioncount_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoiauvcollisioncount_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoiauvcollisioncount_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoiauvcollisioncount_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoiauvcollisioncount_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoithroughput_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoithroughput_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoithroughput_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoithroughput_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoi_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoi_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoi_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoi_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoi_collision_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoi_collision_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoi_collision_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoi_collision_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuvoitotalthroughput_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoitotalthroughput_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoitotalthroughput_vbps_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoitotalthroughput_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                
            end
        end
    end
end
%
% %%save data
% save(fullfile(path2,'20211220vbpspacketlength_vuvoi_vbps_average'),'vuvoi_vbps_average');
% save(fullfile(path2,'20211220vbpspacketlength_vuvoi_collision_vbps_average'),'vuvoi_collision_vbps_average');
%
% save(fullfile(path2,'20211220vbpspacketlength_vuvoitotaldelay_vbps_average'),'vuvoitotaldelay_vbps_average');
% save(fullfile(path2,'20211220vbpspacketlength_vuvoisenddelay_vbps_average'),'vuvoisenddelay_vbps_average');
 % save(fullfile(path2,'20211220vbpspacketlength_vuvoidiscardnodecount_vbps_average'),'vuvoidiscardnodecount_vbps_average');
% save(fullfile(path2,'20211220vbpspacketlength_vuvoicollisioncount_vbps_average'),'vuvoicollisioncount_vbps_average');
% save(fullfile(path2,'20211220vbpspacketlength_vuvoiauvcollisioncount_vbps_average'),'vuvoiauvcollisioncount_vbps_average');
% save(fullfile(path2,'20211220vbpspacketlength_vuvoithroughput_vbps_average'),'vuvoithroughput_vbps_average');
% save(fullfile(path2,'20211220vbpspacketlength_vuvoitotalthroughput_vbps_average'),'vuvoitotalthroughput_vbps_average');
%
%
%
% save(fullfile(path2,'20211220vbpspacketlength_vmvoi_vbps_average'),'vmvoi_vbps_average');
% save(fullfile(path2,'20211220vbpspacketlength_vmvoi_collision_vbps_average'),'vmvoi_collision_vbps_average');
%
% save(fullfile(path2,'20211220vbpspacketlength_vmvoitotaldelay_vbps_average'),'vmvoitotaldelay_vbps_average');
% save(fullfile(path2,'20211220vbpspacketlength_vmvoisenddelay_vbps_average'),'vmvoisenddelay_vbps_average');
 % save(fullfile(path2,'20211220vbpspacketlength_vmvoidiscardnodecount_vbps_average'),'vmvoidiscardnodecount_vbps_average');
% save(fullfile(path2,'20211220vbpspacketlength_vmvoicollisioncount_vbps_average'),'vmvoicollisioncount_vbps_average');
% save(fullfile(path2,'20211220vbpspacketlength_vmvoiauvcollisioncount_vbps_average'),'vmvoiauvcollisioncount_vbps_average');
% save(fullfile(path2,'20211220vbpspacketlength_vmvoithroughput_vbps_average'),'vmvoithroughput_vbps_average');
% save(fullfile(path2,'20211220vbpspacketlength_vmvoitotalthroughput_vbps_average'),'vmvoitotalthroughput_vbps_average');




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

for iL_packet = 1:length(L_packetrange)%1:length(vmtrafficload)
    A1(iL_packet) = sum(vuvoicollisioncount_vbps_average(iL_packet,:));
    
end


for iL_packet = 1:length(L_packetrange)%1:length(vmtrafficload)
    
    A5(iL_packet) = sum(vmvoicollisioncount_vbps_average(iL_packet,:));
end





figure1 = figure;

plot(500:500:5000,A1+A5,'s-b','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));


backColor  =  [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',500:500:5000);
xlabel('\fontsize{16} Packet Sizes (Bits) ');
ylabel('\fontsize{16} Collision Probability for Static nodes');




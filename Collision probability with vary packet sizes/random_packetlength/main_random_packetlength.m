%%%%%%%%%%%%%%%%%%%%%% Main file for random with vary packet sizes

clc;clear;close all;

%%-------------------------Parameters----------------------------------------
dc = 3500;%m
f = 26000;%Hz
R = 13900;%bps
vs = 1500;%m/s
v_A = 10;%m/s


%% Generate locations for predeployed static nodes
N0=3;
N1=N0*5;
N2=N0*5*3;
N=N0+N1+N2;

path1 = 'C:\Users\zxx\Desktop\upload\Collision probability with vary packet sizes\random_packetlength\networkdata'; %networkdata file path
path2 = 'C:\Users\zxx\Desktop\upload\Collision probability with vary packet sizes\random_packetlength\plot_data'; %plotdata file path


tic
[Nx,Ny,Nz] = generate_staticnode_location(N0,N1,N2,dc,path1); %Generate locations for predeployed static nodes
toc



%% Path Preplanning for AUVs
NA=2;
tic
[Ax,Ay,Az, len_xyz,timeplus] = generate_AUVpath(Nx,Ny,Nz, NA,path1);
toc
%%
iteration_total=2000;
N0_length=1;
N0range=3;
L_packetrange=500:500:5000;
for iteration=1:iteration_total
    for  iL_packet=1:length(L_packetrange)
        L_packet=L_packetrange(iL_packet);%bit
        t_data=L_packet/R;
        n=ceil((dc/vs+t_data)/t_data);
        
        %% Generate the network topology and sending time for each static nodes based on TDMA protocol
        Round   =   50;  %Communication round
        max_pathcount   =   120;
        totaltimeslot   =   Round*max_pathcount;
        [delay_table,dis_table,TM,IM,timeslot,sequence,node,path,Ts]   =   generate_TDMA(N,Nx,Ny,Nz,dc,vs,N0,N1,N2,path1,n,t_data,Round,totaltimeslot,iL_packet);
        
        
        %% Static node localization
        tic
        [AUV_listen,AUV_detected,AUV_in] =  Static_node_localization(Ax, Ay, Az,Ts, Nx, Ny, Nz,n,t_data,path1,N0,N1,N2,NA,N,Round,totaltimeslot,timeplus,path,node,len_xyz,dc,iL_packet);
        toc
        
        %% The network throughput of static nodes
        timestart=10;
        timeend=20;
        
        totaltime=(len_xyz.*timeplus-timestart-timeend);
        slotnumber_static=floor(totaltime./n./t_data);
        
        sendpacketnumber_static=slotnumber_static.*(N0+N1);
        
        
        %% Transmission Scheduling
        
        vutrafficload=0.05;
        vmtrafficload=0.5;
        
        
        
        
        betau=0.5;
        vu0=L_packet;
        alphau=3;
        Tu=5*n*t_data;
        
        
        betam=0.3;
        vm0=L_packet;
        alpham=3;
        Tm=10*n*t_data;
        
        
        
        
        for ivutrafficload=1:length(vutrafficload)  %pkt/s
            
            vulambda=vutrafficload(ivutrafficload);
            [vuAUVrun_detected,vuarrivaltime,vuAUVrun_in,vuAUVrun_knownnodes,vuAUVrun_knownnodesin,vuAUVrun_unknownnodes,vuAUVrun_unknownnodesin,vuAUVtimebefore,vuAUVtimeplus,vutimeslot_number,vutimeslot_numberplus,vudiscardnodecount,vuknown_nodescount,vuknown_nodesincount,vuqueuearrivaltime,vuqueuecount,vuqueueindex,vuqueueout,vuqueuesendtime,vusendauvrun_time_random,vusenddelay,vutime_interval,vutotaldelay, vuunknown_nodescount,vuunknown_nodesincount,vuAUVrun_distancenode,vucollisioncount,vutotaltimearrival,vuvoi,vusendauvrun_node_random,vuauvcollisioncount,vuvoi_collision,vupacketerrorate] = randomalgorithm(Nx,Ny,Nz,node,path,timeslot,len_xyz,dis_table,Ax,Ay,Az,AUV_in,AUV_detected,vulambda,NA,N,timeplus,n,t_data,N0,N1,N2,vs,vutrafficload,ivutrafficload,iteration,betau,vu0,alphau,Tu,L_packet,timestart,timeend);
            for ivmtrafficload=1:length(vmtrafficload)
                vmlambda=vmtrafficload(ivmtrafficload); %pkt/s
                [vmAUVrun_detected,vmarrivaltime,vmAUVrun_in,vmAUVrun_knownnodes,vmAUVrun_knownnodesin,vmAUVrun_unknownnodes,vmAUVrun_unknownnodesin,vmAUVtimebefore,vmAUVtimeplus,vmtimeslot_number,vmtimeslot_numberplus,vmdiscardnodecount,vmknown_nodescount,vmknown_nodesincount,vmqueuearrivaltime,vmqueuecount,vmqueueindex,vmqueueout,vmqueuesendtime,vmsendauvrun_time_random,vmsenddelay,vmtime_interval,vmtotaldelay, vmunknown_nodescount, vmunknown_nodesincount,vmAUVrun_distancenode,vmcollisioncount,vmtotaltimearrival,vmvoi,vmsendauvrun_node_random,vmauvcollisioncount,vmvoi_collision,vmpacketerrorate] = randomalgorithm(Nx,Ny,Nz,node,path,timeslot,len_xyz,dis_table,Ax,Ay,Az,AUV_in,AUV_detected,vmlambda,NA,N,timeplus,n,t_data,N0,N1,N2,vs,vmtrafficload,ivmtrafficload,iteration,betam,vm0,alpham,Tm,L_packet,timestart,timeend);
                
                %% Calculate VoI, network throughput, congestion ratio, and  network delay
                for iNA_plot=1:NA
                    vuarrivaltime_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0; % for vu
                    vutotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vutotaldelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vusenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vusenddelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vudiscardnodecount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vucollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vuauvcollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;%%zxx1207
                    
                    vuthroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vutotalthroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    for ivuarrival=1:vutotaltimearrival(iNA_plot,ivutrafficload)
                        if vutotaldelay(iNA_plot,ivuarrival)<1000
                            vuarrivaltime_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vuarrivaltime_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                            vutotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vutotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vutotaldelay(iNA_plot,ivuarrival);
                            vusenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vusenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vusenddelay(iNA_plot,ivuarrival);
                        end
                        
                    end
                    
                    vutotaldelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vutotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuarrivaltime_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                    vusenddelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vusenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuarrivaltime_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                    
                    vudiscardnodecount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=sum(vudiscardnodecount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload); %周围没有节点，也没法传输 rate
                    vucollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=sum(vucollisioncount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload);%由于不是globalmap导致冲突 rate
                    vuauvcollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=sum(vuauvcollisioncount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload);%由于不是globalmap导致冲突 rate
                    
                    vuthroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=(vuarrivaltime_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)-sum(vudiscardnodecount(iNA_plot,:))-sum(vucollisioncount(iNA_plot,:))-sum(vuauvcollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);%吞吐量
                    vutotalthroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=((len_xyz(iNA_plot)*timeplus-timestart-timeend)/n/t_data*(N0+N1)-sum(vucollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);%系统总吞吐量
                    
                end
                
                
                for iNA_plot=1:NA
                    vmarrivaltime_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0; % for vu
                    vmtotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vmtotaldelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vmsenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vmsenddelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    
                    vmdiscardnodecount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vmcollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vmauvcollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    
                    vmthroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vmtotalthroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    for ivmarrival=1:vmtotaltimearrival(iNA_plot,ivmtrafficload)
                        if vmtotaldelay(iNA_plot,ivmarrival)<1000
                            vmarrivaltime_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmarrivaltime_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                            vmtotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmtotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmtotaldelay(iNA_plot,ivmarrival);
                            vmsenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmsenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmsenddelay(iNA_plot,ivmarrival);
                        end
                    end
                    vmtotaldelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmtotaldelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmarrivaltime_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                    vmsenddelay_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vusenddelay_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmarrivaltime_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                    
                    vmdiscardnodecount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=sum(vmdiscardnodecount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload); %周围没有节点，也没法传输 rate
                    vmcollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=sum(vmcollisioncount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload);%由于不是globalmap导致冲突 rate
                    vmauvcollisioncount_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=sum(vmauvcollisioncount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload);%由于不是globalmap导致冲突 rate
                    vmthroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=(vmarrivaltime_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)-sum(vmcollisioncount(iNA_plot,:))-sum(vmauvcollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);%吞吐量
                    vmtotalthroughput_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=((len_xyz(iNA_plot)*timeplus-timestart-timeend)/n/t_data*(N0+N1)-sum(vmcollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);%系统总吞吐量
                end
                
                
                
                for iNA_plot=1:NA
                    vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0; % for vu
                    vuvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vuvoi_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vuvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vuvoi_collision_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    
                    for ivuarrival=2:vutotaltimearrival(iNA_plot,ivutrafficload)
                        if vutotaldelay(iNA_plot,ivuarrival)<1000
                            vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                            vuvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vuvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoi(iNA_plot,ivuarrival);
                            vuvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vuvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoi_collision(iNA_plot,ivuarrival);
                        end
                    end
                    vuvoi_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vuvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                    vuvoi_collision_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vuvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                    
                end
                
                
                
                
                for iNA_plot=1:NA
                    vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0; % for vu
                    vmvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vmvoi_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vmvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    vmvoi_collision_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                    
                    for ivmarrival=2:vmtotaltimearrival(iNA_plot,ivmtrafficload)
                        if vmtotaldelay(iNA_plot,ivmarrival)<1000
                            vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                            vmvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoi(iNA_plot,ivmarrival);
                            vmvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoi_collision(iNA_plot,ivmarrival);
                        end
                    end
                    vmvoi_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmvoi_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                    vmvoi_collision_plot(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmvoi_collision_plot_sum(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iL_packet,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                end
                
                
            end
            
        end
        
    end
    % save(fullfile(path2,'20211220randompacketlength_vuauvcollisioncount'),'vuauvcollisioncount');
    % save(fullfile(path2,'20211220randompacketlength_vusendauvrun_node_random'),'vusendauvrun_node_random');
    %
    % save(fullfile(path2,'20211220randompacketlength_vutotaldelay'),'vutotaldelay');
    % save(fullfile(path2,'20211220randompacketlength_vutotaldelay_plot'),'vutotaldelay_plot');
    % save(fullfile(path2,'20211220randompacketlength_vusenddelay'),'vusenddelay');
    % save(fullfile(path2,'20211220randompacketlength_vusenddelay_plot'),'vusenddelay_plot');
    %
    % save(fullfile(path2,'20211220randompacketlength_vudiscardnodecount_plot'),'vudiscardnodecount_plot');
    % save(fullfile(path2,'20211220randompacketlength_vucollisioncount_plot'),'vucollisioncount_plot');
    % save(fullfile(path2,'20211220randompacketlength_vuauvcollisioncount_plot'),'vuauvcollisioncount_plot');
    %
    % save(fullfile(path2,'20211220randompacketlength_vuthroughput_plot'),'vuthroughput_plot');
    % save(fullfile(path2,'20211220randompacketlength_vutotalthroughput_plot'),'vutotalthroughput_plot');
    %
    %
    % save(fullfile(path2,'20211220randompacketlength_vuarrivaltime_plot'),'vuarrivaltime_plot');
    %
    % save(fullfile(path2,'20211220randompacketlength_vutotaldelay_plot_sum'),'vutotaldelay_plot_sum');
    % save(fullfile(path2,'20211220randompacketlength_vusenddelay_plot_sum'),'vusenddelay_plot_sum');
    %
    % save(fullfile(path2,'20211220randompacketlength_vmauvcollisioncount'),'vmauvcollisioncount');
    % save(fullfile(path2,'20211220randompacketlength_vmsendauvrun_node_random'),'vmsendauvrun_node_random');
    %
    % save(fullfile(path2,'20211220randompacketlength_vmtotaldelay'),'vmtotaldelay');
    % save(fullfile(path2,'20211220randompacketlength_vmtotaldelay_plot'),'vmtotaldelay_plot');
    % save(fullfile(path2,'20211220randompacketlength_vmsenddelay'),'vmsenddelay');
    % save(fullfile(path2,'20211220randompacketlength_vmsenddelay_plot'),'vmsenddelay_plot');
    %
    % save(fullfile(path2,'20211220randompacketlength_vmdiscardnodecount_plot'),'vmdiscardnodecount_plot');
    % save(fullfile(path2,'20211220randompacketlength_vmcollisioncount_plot'),'vmcollisioncount_plot');
    % save(fullfile(path2,'20211220randompacketlength_vmauvcollisioncount_plot'),'vmauvcollisioncount_plot');
    %
    % save(fullfile(path2,'20211220randompacketlength_vmthroughput_plot'),'vmthroughput_plot');
    % save(fullfile(path2,'20211220randompacketlength_vmtotalthroughput_plot'),'vmtotalthroughput_plot');
    %
    % save(fullfile(path2,'20211220randompacketlength_vmarrivaltime'),'vmarrivaltime');
    %
    % save(fullfile(path2,'20211220randompacketlength_vmtotaldelay_plot_sum'),'vmtotaldelay_plot_sum');
    % save(fullfile(path2,'20211220randompacketlength_vmsenddelay_plot_sum'),'vmsenddelay_plot_sum');
    %
    % save(fullfile(path2,'20211220randompacketlength_vuvoi'),'vuvoi');
    % save(fullfile(path2,'20211220randompacketlength_vuvoi_collision'),'vuvoi_collision');
    % save(fullfile(path2,'20211220randompacketlength_vuvoi_plot'),'vuvoi_plot');
    % save(fullfile(path2,'20211220randompacketlength_vuvoi_collision_plot'),'vuvoi_collision_plot');
    % save(fullfile(path2,'20211220randompacketlength_vuvoi_plot_sum'),'vuvoi_plot_sum');
    % save(fullfile(path2,'20211220randompacketlength_vuvoi_collision_plot_sum'),'vuvoi_collision_plot_sum');
    %
    %
    % save(fullfile(path2,'20211220randompacketlength_vmvoi'),'vmvoi');
    % save(fullfile(path2,'20211220randompacketlength_vmvoi_collision'),'vmvoi_collision');
    % save(fullfile(path2,'20211220randompacketlength_vmvoi_plot'),'vmvoi_plot');
    % save(fullfile(path2,'20211220randompacketlength_vmvoi_collision_plot'),'vmvoi_collision_plot');
    % save(fullfile(path2,'20211220randompacketlength_vmvoi_plot_sum'),'vmvoi_plot_sum');
    % save(fullfile(path2,'20211220randompacketlength_vmvoi_collision_plot_sum'),'vmvoi_collision_plot_sum');
    
    iteration
end


%Average
for iL_packet=1:length(L_packetrange)
    for iNA_plot=1:NA
        for ivutrafficload=1:length(vutrafficload)
            for ivmtrafficload=1:length(vmtrafficload)
                vutotaldelay_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vutotaldelay_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmtotaldelay_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmtotaldelay_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vusenddelay_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vusenddelay_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmsenddelay_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmsenddelay_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                
                vudiscardnodecount_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vudiscardnodecount_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmdiscardnodecount_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmdiscardnodecount_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vucollisioncount_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vucollisioncount_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmcollisioncount_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmcollisioncount_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuauvcollisioncount_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vuauvcollisioncount_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmauvcollisioncount_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmauvcollisioncount_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vuthroughput_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vuthroughput_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmthroughput_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmthroughput_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vutotalthroughput_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vutotalthroughput_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmtotalthroughput_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmtotalthroughput_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                
                vuvoi_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vuvoi_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoi_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmvoi_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                
                vuvoi_collision_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vuvoi_collision_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
                vmvoi_collision_random_average(iL_packet,iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmvoi_collision_plot(iL_packet,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            end
        end
    end
end


% %%save data
% save(fullfile(path2,'20211220randompacketlength_vutotaldelay_random_average'),'vutotaldelay_random_average');
% save(fullfile(path2,'20211220randompacketlength_vusenddelay_random_average'),'vusenddelay_random_average');
%
% save(fullfile(path2,'20211220randompacketlength_vudiscardnodecount_random_average'),'vudiscardnodecount_random_average');
% save(fullfile(path2,'20211220randompacketlength_vuauvcollisioncount_random_average'),'vuauvcollisioncount_random_average');
% save(fullfile(path2,'20211220randompacketlength_vucollisioncount_random_average'),'vucollisioncount_random_average');
% save(fullfile(path2,'20211220randompacketlength_vuthroughput_random_average'),'vuthroughput_random_average');
% save(fullfile(path2,'20211220randompacketlength_vutotalthroughput_random_average'),'vutotalthroughput_random_average');
%
%
%
% save(fullfile(path2,'20211220randompacketlength_vmtotaldelay_random_average'),'vmtotaldelay_random_average');
% save(fullfile(path2,'20211220randompacketlength_vmsenddelay_random_average'),'vmsenddelay_random_average');
%
% save(fullfile(path2,'20211220randompacketlength_vmdiscardnodecount_random_average'),'vmdiscardnodecount_random_average');
% save(fullfile(path2,'20211220randompacketlength_vmcollisioncount_random_average'),'vmcollisioncount_random_average');
% save(fullfile(path2,'20211220randompacketlength_vmauvcollisioncount_random_average'),'vmauvcollisioncount_random_average');
% save(fullfile(path2,'20211220randompacketlength_vmthroughput_random_average'),'vmthroughput_random_average');
% save(fullfile(path2,'20211220randompacketlength_vmtotalthroughput_random_average'),'vmtotalthroughput_random_average');
%
%
%
%
% save(fullfile(path2,'20211220randompacketlength_vuvoi_random_average'),'vuvoi_random_average');
%
% save(fullfile(path2,'20211220randompacketlength_vuvoi_collision_random_average'),'vuvoi_collision_random_average');
%
% save(fullfile(path2,'20211220randompacketlength_vmvoi_random_average'),'vmvoi_random_average');
% save(fullfile(path2,'20211220randompacketlength_vmvoi_collision_random_average'),'vmvoi_collision_random_average');
%



%%
clr(1,:)=[0 113 188]/255;
clr(2,:)=[216 82 24]/255;
clr(3,:)=[236 176 31]/255;
clr(4,:)=[125 46 141]/255;
clr(5,:)=[118 171 47]/255;
clr(6,:)=[76 189 237]/255;
clr(7,:)=[161 19 46]/255;
clr(8,:)=[255 255 255]/255;
clr(9,:)=[25 35 45]/255;
clr(10,:)=[100 100 100]/255;
clr(11,:)=[150 150 150]/255;

for iL_packet=1:length(L_packetrange)%1:length(vmtrafficload)
    A1(iL_packet)=sum(vucollisioncount_random_average(iL_packet,:));
    
end


for iL_packet=1:length(L_packetrange)%1:length(vmtrafficload)
    
    A5(iL_packet)=sum(vmcollisioncount_random_average(iL_packet,:));
end





figure1=figure;

plot(500:500:5000,A1+A5,'s-b','LineWidth',1.8,'MarkerSize',8,'color',clr(1,:));


backColor = [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',500:500:5000);
xlabel('\fontsize{16} Packet Sizes (Bits) ');
ylabel('\fontsize{16} Collision Probability for Static nodes');


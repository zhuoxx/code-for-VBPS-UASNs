
%%%%%%%%%%%%%%%%%%%%%% Main file for Random with different numbers of AUVs

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

iteration_total =1;% 2000;
NArange = 1:5;


%%%%%%%%%%%%%%%%%%%%%%%%%simulation network in this paper, it is better not
%%%%%%%%%%%%%%%%%%%%%%%%%to modifiy this part%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate locations for predeployed static nodes
N0 = 3;        %Number of sink nodes
N1 = N0*5;     %Number of relay nodes
N2 = N0*5*3;   %Number of sensor nodes
N = N0+N1+N2;  %Total static node number

path1 = 'C:\Users\zxx\Desktop\upload\Network performance with different numbers of AUVs\random_AUVnumber\networkdata'; %networkdata file path
path2 = 'C:\Users\zxx\Desktop\upload\Network performance with different numbers of AUVs\random_AUVnumber\plot_data'; %plotdata file path


tic
[Nx,Ny,Nz]  =  generate_staticnode_location(N0,N1,N2,dc,path1); %Generate locations for predeployed static nodes
toc

%% Generate the network topology and sending time for each static nodes based on TDMA protocol
Round   =   50;  %Communication round
max_pathcount   =   120;
totaltimeslot  =  Round*max_pathcount;
[delay_table,dis_table,TM,IM,timeslot,sequence,node,path,Ts]  =  generate_TDMA(N,Nx,Ny,Nz,dc,vs,N0,N1,N2,path1,n,t_data,Round,totaltimeslot);
 


for iteration = 1:iteration_total
    for iNA_number = 1:5

        %% Path Preplanning for AUVs
        NA = iNA_number;
        tic
        [Ax,Ay,Az, len_xyz,timeplus]  =  generate_AUVpath(Nx,Ny,Nz, NA,path1);
        toc

        %% Static node localization
        tic
        [AUV_listen,AUV_detected,AUV_in]  =   Static_node_localization(Ax, Ay, Az,Ts, Nx, Ny, Nz,n,t_data,path1,N0,N1,N2,NA,N,Round,totaltimeslot,timeplus,path,node,len_xyz,dc);
        toc

        %% The network throughput of static nodes
        timestart = 10;
        timeend = 20;

        totaltime = (len_xyz.*timeplus-timestart-timeend);
        slotnumber_static = floor(totaltime./n./t_data);

        sendpacketnumber_static = slotnumber_static.*(N0+N1);%每个timeslot都有N0和N1在接收信号


        %% Transmission Scheduling

        vutrafficload  =  0.05;
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
           [vuAUVrun_detected,vuarrivaltime,vuAUVrun_in,vuAUVrun_knownnodes,vuAUVrun_knownnodesin,vuAUVrun_unknownnodes,vuAUVrun_unknownnodesin,vuAUVtimebefore,vuAUVtimeplus,vutimeslot_number,vutimeslot_numberplus,vudiscardnodecount,vuknown_nodescount,vuknown_nodesincount,vuqueuearrivaltime,vuqueuecount,vuqueueindex,vuqueueout,vuqueuesendtime,vusendauvrun_time_random,vusenddelay,vutime_interval,vutotaldelay, vuunknown_nodescount,vuunknown_nodesincount,vuAUVrun_distancenode,vucollisioncount,vutotaltimearrival,vuvoi,vusendauvrun_node_random,vuauvcollisioncount,vuvoi_collision,vupacketerrorate]  =  randomalgorithm(Nx,Ny,Nz,node,path,timeslot,len_xyz,dis_table,Ax,Ay,Az,AUV_in,AUV_detected,vulambda,NA,N,timeplus,n,t_data,N0,N1,N2,vs,vutrafficload,ivutrafficload,iteration,betau,vu0,alphau,Tu,L_packet,timestart,timeend);
        for ivmtrafficload = 1:length(vmtrafficload)    
            vmlambda = vmtrafficload(ivmtrafficload); %pkt/s
           [vmAUVrun_detected,vmarrivaltime,vmAUVrun_in,vmAUVrun_knownnodes,vmAUVrun_knownnodesin,vmAUVrun_unknownnodes,vmAUVrun_unknownnodesin,vmAUVtimebefore,vmAUVtimeplus,vmtimeslot_number,vmtimeslot_numberplus,vmdiscardnodecount,vmknown_nodescount,vmknown_nodesincount,vmqueuearrivaltime,vmqueuecount,vmqueueindex,vmqueueout,vmqueuesendtime,vmsendauvrun_time_random,vmsenddelay,vmtime_interval,vmtotaldelay, vmunknown_nodescount, vmunknown_nodesincount,vmAUVrun_distancenode,vmcollisioncount,vmtotaltimearrival,vmvoi,vmsendauvrun_node_random,vmauvcollisioncount,vmvoi_collision,vmpacketerrorate]  =  randomalgorithm(Nx,Ny,Nz,node,path,timeslot,len_xyz,dis_table,Ax,Ay,Az,AUV_in,AUV_detected,vmlambda,NA,N,timeplus,n,t_data,N0,N1,N2,vs,vmtrafficload,ivmtrafficload,iteration,betam,vm0,alpham,Tm,L_packet,timestart,timeend);

                  
                   
             %% Calculate VoI, network throughput, congestion ratio, and  network delay
           for iNA_plot = 1:NA
                vuarrivaltime_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0; % for vu
                vutotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vutotaldelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vusenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vusenddelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                  vudiscardnodecount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vucollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vuauvcollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;%%zxx1207

                vuthroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                 vutotalthroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                for ivuarrival = 1:vutotaltimearrival(iNA_plot,ivutrafficload)
                    if vutotaldelay(iNA_plot,ivuarrival)<1000
                        vuarrivaltime_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuarrivaltime_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                        vutotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vutotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vutotaldelay(iNA_plot,ivuarrival);
                        vusenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vusenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vusenddelay(iNA_plot,ivuarrival);
                          end
                    
                end
                
                vutotaldelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vutotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuarrivaltime_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);    
                vusenddelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vusenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuarrivaltime_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);    
               
                vudiscardnodecount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vudiscardnodecount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload); %周围没有节点，也没法传输 rate
                vucollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vucollisioncount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload);%由于不是globalmap导致冲突 rate
                 vuauvcollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vuauvcollisioncount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload);%由于不是globalmap导致冲突 rate
                 
                vuthroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = (vuarrivaltime_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)-sum(vudiscardnodecount(iNA_plot,:))-sum(vucollisioncount(iNA_plot,:))-sum(vuauvcollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);%吞吐量
                vutotalthroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = ((len_xyz(iNA_plot)*timeplus-timestart-timeend)/n/t_data*(N0+N1)-sum(vucollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);%系统总吞吐量

            end

            
            for iNA_plot = 1:NA
                vmarrivaltime_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0; % for vu
                vmtotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vmtotaldelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vmsenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vmsenddelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
               
                vmdiscardnodecount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vmcollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vmauvcollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;

                vmthroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vmtotalthroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                for ivmarrival = 1:vmtotaltimearrival(iNA_plot,ivmtrafficload)
                    if vmtotaldelay(iNA_plot,ivmarrival)<1000
                        vmarrivaltime_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmarrivaltime_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                        vmtotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmtotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmtotaldelay(iNA_plot,ivmarrival);
                          vmsenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmsenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmsenddelay(iNA_plot,ivmarrival);
                     end
                end
                vmtotaldelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmtotaldelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmarrivaltime_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
               vmsenddelay_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vusenddelay_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmarrivaltime_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);    
                
                vmdiscardnodecount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vmdiscardnodecount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload); %周围没有节点，也没法传输 rate
                vmcollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vmcollisioncount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload);%由于不是globalmap导致冲突 rate
                vmauvcollisioncount_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = sum(vmauvcollisioncount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload);%由于不是globalmap导致冲突 rate
                 vmthroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = (vmarrivaltime_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)-sum(vmcollisioncount(iNA_plot,:))-sum(vmauvcollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);%吞吐量
                vmtotalthroughput_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = ((len_xyz(iNA_plot)*timeplus-timestart-timeend)/n/t_data*(N0+N1)-sum(vmcollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);%系统总吞吐量
            end
            
            %% voi和delay的计算
            
            for iNA_plot = 1:NA
                vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0; % for vu
                vuvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vuvoi_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                 vuvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vuvoi_collision_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                
                for ivuarrival = 2:vutotaltimearrival(iNA_plot,ivutrafficload)
                    if vutotaldelay(iNA_plot,ivuarrival)<1000
                        vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                        vuvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoi(iNA_plot,ivuarrival);
                         vuvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoi_collision(iNA_plot,ivuarrival);
                    end
                end
                vuvoi_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
               vuvoi_collision_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vuvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                
            end
            
            
            
            
            for iNA_plot = 1:NA
                vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0; % for vu
                vmvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vmvoi_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                 vmvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                vmvoi_collision_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = 0;
                  
                for ivmarrival = 2:vmtotaltimearrival(iNA_plot,ivmtrafficload)
                    if vmtotaldelay(iNA_plot,ivmarrival)<1000
                        vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                        vmvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoi(iNA_plot,ivmarrival);
                        vmvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoi_collision(iNA_plot,ivmarrival);
                    end
                end
                 vmvoi_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoi_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                 vmvoi_collision_plot(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload) = vmvoi_collision_plot_sum(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iNA_number,iNA_plot,iteration,ivutrafficload,ivmtrafficload);
            end
            
            
        end
        
    end

end
% save(fullfile(path2,'20211220randomauvnodenumbers_vuauvcollisioncount'),'vuauvcollisioncount');
% save(fullfile(path2,'20211220randomauvnodenumbers_vusendauvrun_node_random'),'vusendauvrun_node_random');
% 
% save(fullfile(path2,'20211220randomauvnodenumbers_vutotaldelay'),'vutotaldelay');
% save(fullfile(path2,'20211220randomauvnodenumbers_vutotaldelay_plot'),'vutotaldelay_plot');
% save(fullfile(path2,'20211220randomauvnodenumbers_vusenddelay'),'vusenddelay');
% save(fullfile(path2,'20211220randomauvnodenumbers_vusenddelay_plot'),'vusenddelay_plot');
 % 
% save(fullfile(path2,'20211220randomauvnodenumbers_vudiscardnodecount_plot'),'vudiscardnodecount_plot');
% save(fullfile(path2,'20211220randomauvnodenumbers_vucollisioncount_plot'),'vucollisioncount_plot');
% save(fullfile(path2,'20211220randomauvnodenumbers_vuauvcollisioncount_plot'),'vuauvcollisioncount_plot');
% 
% save(fullfile(path2,'20211220randomauvnodenumbers_vuthroughput_plot'),'vuthroughput_plot');
% save(fullfile(path2,'20211220randomauvnodenumbers_vutotalthroughput_plot'),'vutotalthroughput_plot');
% 
% 
% save(fullfile(path2,'20211220randomauvnodenumbers_vuarrivaltime_plot'),'vuarrivaltime_plot');
% 
% save(fullfile(path2,'20211220randomauvnodenumbers_vutotaldelay_plot_sum'),'vutotaldelay_plot_sum');
% save(fullfile(path2,'20211220randomauvnodenumbers_vusenddelay_plot_sum'),'vusenddelay_plot_sum');
 % 
% save(fullfile(path2,'20211220randomauvnodenumbers_vmauvcollisioncount'),'vmauvcollisioncount');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmsendauvrun_node_random'),'vmsendauvrun_node_random');
% 
% save(fullfile(path2,'20211220randomauvnodenumbers_vmtotaldelay'),'vmtotaldelay');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmtotaldelay_plot'),'vmtotaldelay_plot');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmsenddelay'),'vmsenddelay');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmsenddelay_plot'),'vmsenddelay_plot');
 % 
% save(fullfile(path2,'20211220randomauvnodenumbers_vmdiscardnodecount_plot'),'vmdiscardnodecount_plot');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmcollisioncount_plot'),'vmcollisioncount_plot');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmauvcollisioncount_plot'),'vmauvcollisioncount_plot');
% 
% save(fullfile(path2,'20211220randomauvnodenumbers_vmthroughput_plot'),'vmthroughput_plot');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmtotalthroughput_plot'),'vmtotalthroughput_plot');
% 
% save(fullfile(path2,'20211220randomauvnodenumbers_vmarrivaltime'),'vmarrivaltime');
% 
% save(fullfile(path2,'20211220randomauvnodenumbers_vmtotaldelay_plot_sum'),'vmtotaldelay_plot_sum');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmsenddelay_plot_sum'),'vmsenddelay_plot_sum');
 % 
% save(fullfile(path2,'20211220randomauvnodenumbers_vuvoi'),'vuvoi');
% save(fullfile(path2,'20211220randomauvnodenumbers_vuvoi_collision'),'vuvoi_collision');
% save(fullfile(path2,'20211220randomauvnodenumbers_vuvoi_plot'),'vuvoi_plot');
% save(fullfile(path2,'20211220randomauvnodenumbers_vuvoi_collision_plot'),'vuvoi_collision_plot');
% save(fullfile(path2,'20211220randomauvnodenumbers_vuvoi_plot_sum'),'vuvoi_plot_sum');
% save(fullfile(path2,'20211220randomauvnodenumbers_vuvoi_collision_plot_sum'),'vuvoi_collision_plot_sum');
% 
% 
% save(fullfile(path2,'20211220randomauvnodenumbers_vmvoi'),'vmvoi');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmvoi_collision'),'vmvoi_collision');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmvoi_plot'),'vmvoi_plot');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmvoi_collision_plot'),'vmvoi_collision_plot');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmvoi_plot_sum'),'vmvoi_plot_sum');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmvoi_collision_plot_sum'),'vmvoi_collision_plot_sum');
iteration
end


%average
for iNA_number = 1:5 %%%%zxx1209
for iNA_plot = 1:NA
    for ivutrafficload = 1:length(vutrafficload)
        for ivmtrafficload = 1:length(vmtrafficload)
            vutotaldelay_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vutotaldelay_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmtotaldelay_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmtotaldelay_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
             vusenddelay_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vusenddelay_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmsenddelay_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmsenddelay_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
           
            vudiscardnodecount_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vudiscardnodecount_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmdiscardnodecount_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmdiscardnodecount_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vucollisioncount_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vucollisioncount_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmcollisioncount_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmcollisioncount_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vuauvcollisioncount_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuauvcollisioncount_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmauvcollisioncount_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmauvcollisioncount_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vuthroughput_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuthroughput_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmthroughput_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmthroughput_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vutotalthroughput_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vutotalthroughput_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmtotalthroughput_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmtotalthroughput_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
           
            vuvoi_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoi_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmvoi_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoi_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            
              vuvoi_collision_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vuvoi_collision_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmvoi_collision_random_average(iNA_number,iNA_plot,ivutrafficload,ivmtrafficload) = sum(vmvoi_collision_plot(iNA_number,iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
        end
    end
end
end


% %%delay
% save(fullfile(path2,'20211220randomauvnodenumbers_vutotaldelay_random_average'),'vutotaldelay_random_average');
% save(fullfile(path2,'20211220randomauvnodenumbers_vusenddelay_random_average'),'vusenddelay_random_average');
 % 
% save(fullfile(path2,'20211220randomauvnodenumbers_vudiscardnodecount_random_average'),'vudiscardnodecount_random_average');
% save(fullfile(path2,'20211220randomauvnodenumbers_vuauvcollisioncount_random_average'),'vuauvcollisioncount_random_average');
% save(fullfile(path2,'20211220randomauvnodenumbers_vucollisioncount_random_average'),'vucollisioncount_random_average');
% save(fullfile(path2,'20211220randomauvnodenumbers_vuthroughput_random_average'),'vuthroughput_random_average');
% save(fullfile(path2,'20211220randomauvnodenumbers_vutotalthroughput_random_average'),'vutotalthroughput_random_average');
% 
% 
% 
% save(fullfile(path2,'20211220randomauvnodenumbers_vmtotaldelay_random_average'),'vmtotaldelay_random_average');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmsenddelay_random_average'),'vmsenddelay_random_average');
 % 
% save(fullfile(path2,'20211220randomauvnodenumbers_vmdiscardnodecount_random_average'),'vmdiscardnodecount_random_average');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmcollisioncount_random_average'),'vmcollisioncount_random_average');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmauvcollisioncount_random_average'),'vmauvcollisioncount_random_average');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmthroughput_random_average'),'vmthroughput_random_average');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmtotalthroughput_random_average'),'vmtotalthroughput_random_average');
% 
% 
% %%voi和delay
% 
% 
% save(fullfile(path2,'20211220randomauvnodenumbers_vuvoi_random_average'),'vuvoi_random_average');
% 
 % save(fullfile(path2,'20211220randomauvnodenumbers_vuvoi_collision_random_average'),'vuvoi_collision_random_average');
% 
% save(fullfile(path2,'20211220randomauvnodenumbers_vmvoi_random_average'),'vmvoi_random_average');
% save(fullfile(path2,'20211220randomauvnodenumbers_vmvoi_collision_random_average'),'vmvoi_collision_random_average');
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

    A5(iNA_number) = sum(vmauvcollisioncount_random_average(iNA_number,:));
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
    vuvoi_vbps(iNA_number,:) = vuvoi_collision_random_average(iNA_number,:)+vmvoi_collision_random_average(iNA_number,:);
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




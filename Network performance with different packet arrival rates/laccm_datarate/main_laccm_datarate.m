%%%%%%%%%%%%%%%%%%%%%% Main file for LACCM with different datarate

clc;clear;close all;

%%-------------------------Parameters----------------------------------------
dc = 3500;%m
f = 26000;%Hz
R = 13900;%bps
L_packet = 1000;%bit
t_data = L_packet/R;
vs = 1500;%m/s
n = ceil((dc/vs+t_data)/t_data);

v_A = 10;%m/s

%%%%%%%%%%%%%%%%%%%%%%%%%simulation network in this paper, it is better not
%%%%%%%%%%%%%%%%%%%%%%%%%to modifiy this part%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate static nodes
N0 =   2;        %Number of sink nodes
N1 =   10;       %Number of relay nodes
N2 =   30;       %Number of sensor nodes
N =   N0+N1+N2;  %Total number

path1 = 'C:\Users\zxx\Desktop\upload\Network performance with different packet arrival rates\laccm_datarate\networkdata'; %networkdata file path
path2 = 'C:\Users\zxx\Desktop\upload\Network performance with different packet arrival rates\laccm_datarate\plot_data'; %plotdata file path


tic
[Nx,Ny,Nz] =   generate_staticnode_location(N0,N1,N2,dc,path1); %Generate locations for predeployed static nodes
toc

%% Generate the network topology and sending time for each static nodes based on TDMA protocol
Round = 50;  %Communication round
max_pathcount = 120;
totaltimeslot = Round*max_pathcount;
[delay_table,dis_table,TM,IM,timeslot,sequence,node,path,Ts] = generate_TDMA(N,Nx,Ny,Nz,dc,vs,N0,N1,N2,path1,n,t_data,Round,totaltimeslot);



%% Path Preplanning for AUVs
NA =  2;
tic
[Ax,Ay,Az, len_xyz,timeplus] =   generate_AUVpath(Nx,Ny,Nz, NA,path1);
toc

%% Static node localization
tic
[AUV_listen,AUV_detected,AUV_in] = Static_node_localization(Ax, Ay, Az,Ts, Nx, Ny, Nz,n,t_data,path1,N0,N1,N2,NA,N,Round,totaltimeslot,timeplus,path,node,len_xyz,dc);
toc 
%%%%%%%%%%%%%%%%%%%%%%%%%simulation network in this paper, it is better not
%%%%%%%%%%%%%%%%%%%%%%%%%to modifiy the code in this part%%%%%%%%%%%%%%%%%%
%% The network throughput of static nodes
timestart=100;
timeend=500;

totaltime=(len_xyz.*timeplus-timestart-timeend);
slotnumber_static=floor(totaltime./n./t_data);

sendpacketnumber_static=slotnumber_static.*(N0+N1); %The sink nodes and relay nodes receive data packet at each time slot
%% Topology construction and transmission scheduling


%%--------------------------Parameters----------------------
%The packet arrival rate
tic
vutrafficload=0.01:0.01:0.1;% for abnormal data
vmtrafficload=0.1:0.1:1;%for normal data

iteration_total=2000;

betau=0.5;
vu0=L_packet;
alphau=3;
Tu=5*n*t_data;


betam=0.3;
vm0=L_packet;
alpham=3;
Tm=10*n*t_data;

%%------------Toplogy construction and transmission scheduling-------------
for iteration=1:iteration_total

    for ivutrafficload=1:length(vutrafficload)  %pkt/s
        for ivmtrafficload=1:length(vmtrafficload)
            %delaycomputing function for vu
            vulambda=vutrafficload(ivutrafficload);
            [vuAUVrun_detected,vuarrivaltime,vuAUVrun_in,vuAUVrun_knownnodes,vuAUVrun_knownnodesin,vuAUVrun_unknownnodes,vuAUVrun_unknownnodesin,vuAUVtimebefore,vuAUVtimeplus,vutimeslot_number,vutimeslot_numberplus,vudiscardnodecount,vuknown_nodescount,vuknown_nodesincount,vuqueuearrivaltime,vuqueuecount,vuqueueindex,vuqueueout,vuqueuesendtime,vusendauvrun_time_random,vusenddelay,vutime_interval,vutotaldelay, vuunknown_nodescount,vuunknown_nodesincount,vuAUVrun_distancenode,vucollisioncount,vutotaltimearrival,vuvoi,vusendauvrun_node_random,vutotaldiscardnodecount,vuauvcollisioncount,vuvoi_collision] = laccmalgorithm(Nx,Ny,Nz,node,path,timeslot,len_xyz,dis_table,Ax,Ay,Az,AUV_in,AUV_detected,vulambda,NA,N,timeplus,n,t_data,N0,N1,N2,vs,vutrafficload,ivutrafficload,iteration,betau,vu0,alphau,Tu,L_packet,timestart,timeend);

            vmlambda=vmtrafficload(ivmtrafficload); %pkt/s
            [vmAUVrun_detected,vmarrivaltime,vmAUVrun_in,vmAUVrun_knownnodes,vmAUVrun_knownnodesin,vmAUVrun_unknownnodes,vmAUVrun_unknownnodesin,vmAUVtimebefore,vmAUVtimeplus,vmtimeslot_number,vmtimeslot_numberplus,vmdiscardnodecount,vmknown_nodescount,vmknown_nodesincount,vmqueuearrivaltime,vmqueuecount,vmqueueindex,vmqueueout,vmqueuesendtime,vmsendauvrun_time_random,vmsenddelay,vmtime_interval,vmtotaldelay, vmunknown_nodescount, vmunknown_nodesincount,vmAUVrun_distancenode,vmcollisioncount,vmtotaltimearrival,vmvoi,vmsendauvrun_node_random,vmtotaldiscardnodecount,vmauvcollisioncount,vmvoi_collision] = laccmalgorithm(Nx,Ny,Nz,node,path,timeslot,len_xyz,dis_table,Ax,Ay,Az,AUV_in,AUV_detected,vmlambda,NA,N,timeplus,n,t_data,N0,N1,N2,vs,vmtrafficload,ivmtrafficload,iteration,betam,vm0,alpham,Tm,L_packet,timestart,timeend);

            %% Calculate VoI, network throughput, congestion ratio, and  network delay
            for iNA_plot=1:NA
                vuarrivaltime_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0; % for vu
                vutotaldelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vutotaldelay_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vusenddelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vusenddelay_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vudiscardnodecount_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vucollisioncount_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vuauvcollisioncount_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;

                vuthroughput_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vutotalthroughput_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;




                for ivuarrival=1:vutotaltimearrival(iNA_plot,ivutrafficload)
                    if vutotaldelay(iNA_plot,ivuarrival)<1000
                        vuarrivaltime_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vuarrivaltime_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                        vutotaldelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vutotaldelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vutotaldelay(iNA_plot,ivuarrival);
                        vusenddelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vusenddelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vusenddelay(iNA_plot,ivuarrival);

                    end

                end

                vutotaldelay_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vutotaldelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuarrivaltime_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                vusenddelay_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vusenddelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuarrivaltime_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload);


                vudiscardnodecount_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=sum(vudiscardnodecount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload);
                vucollisioncount_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=sum(vucollisioncount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload);
                vuauvcollisioncount_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=sum(vuauvcollisioncount(iNA_plot,:))/vutotaltimearrival(iNA_plot,ivutrafficload);

                vuthroughput_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=(vuarrivaltime_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)-sum(vucollisioncount(iNA_plot,:))-sum(vuauvcollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);
                vutotalthroughput_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=((len_xyz(iNA_plot)*timeplus-timestart-timeend)/n/t_data*(N0+N1)-sum(vucollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);

            end





            for iNA_plot=1:NA
                vmarrivaltime_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0; % for vu
                vmtotaldelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vmtotaldelay_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vmsenddelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vmsenddelay_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;


                vmdiscardnodecount_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vmcollisioncount_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vmauvcollisioncount_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;

                vmthroughput_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vmtotalthroughput_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                for ivmarrival=1:vmtotaltimearrival(iNA_plot,ivmtrafficload)
                    if vmtotaldelay(iNA_plot,ivmarrival)<1000
                        vmarrivaltime_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmarrivaltime_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                        vmtotaldelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmtotaldelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmtotaldelay(iNA_plot,ivmarrival);
                        vmsenddelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmsenddelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmsenddelay(iNA_plot,ivmarrival);
                    end
                end
                vmtotaldelay_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmtotaldelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmarrivaltime_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                vmsenddelay_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmsenddelay_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmarrivaltime_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload);

                vmdiscardnodecount_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=sum(vmdiscardnodecount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload);
                vmcollisioncount_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=sum(vmcollisioncount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload);
                vmauvcollisioncount_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=sum(vmauvcollisioncount(iNA_plot,:))/vmtotaltimearrival(iNA_plot,ivmtrafficload);%

                vmthroughput_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=(vmarrivaltime_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)-sum(vmcollisioncount(iNA_plot,:))-sum(vmauvcollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);
                vmtotalthroughput_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=((len_xyz(iNA_plot)*timeplus-timestart-timeend)/n/t_data*(N0+N1)-sum(vmcollisioncount(iNA_plot,:)))*t_data/(len_xyz(iNA_plot)*timeplus-timestart-timeend);

            end



            for iNA_plot=1:NA
                vuvoiarrivaltime(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0; % for vu
                vuvoi_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vuvoi_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vuvoi_collision_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vuvoi_collision_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;


                for ivuarrival=2:vutotaltimearrival(iNA_plot,ivutrafficload)
                    if vutotaldelay(iNA_plot,ivuarrival)<1000
                        vuvoiarrivaltime(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vuvoiarrivaltime(iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                        vuvoi_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vuvoi_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoi(iNA_plot,ivuarrival);
                        vuvoi_collision_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vuvoi_collision_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vuvoi_collision(iNA_plot,ivuarrival);
                    end
                end
                vuvoi_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vuvoi_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                vuvoi_collision_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vuvoi_collision_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vuvoiarrivaltime(iNA_plot,iteration,ivutrafficload,ivmtrafficload);

            end




            for iNA_plot=1:NA
                vmvoiarrivaltime(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0; % for vu
                vmvoi_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vmvoi_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vmvoi_collision_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;
                vmvoi_collision_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=0;

                for ivmarrival=2:vmtotaltimearrival(iNA_plot,ivmtrafficload)
                    if vmtotaldelay(iNA_plot,ivmarrival)<1000
                        vmvoiarrivaltime(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmvoiarrivaltime(iNA_plot,iteration,ivutrafficload,ivmtrafficload)+1;
                        vmvoi_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmvoi_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoi(iNA_plot,ivmarrival);
                        vmvoi_collision_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmvoi_collision_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)+vmvoi_collision(iNA_plot,ivmarrival);
                    end
                end
                vmvoi_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmvoi_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iNA_plot,iteration,ivutrafficload,ivmtrafficload);
                vmvoi_collision_plot(iNA_plot,iteration,ivutrafficload,ivmtrafficload)=vmvoi_collision_plot_sum(iNA_plot,iteration,ivutrafficload,ivmtrafficload)/vmvoiarrivaltime(iNA_plot,iteration,ivutrafficload,ivmtrafficload);
            end


        end

    end
    %% save data
    % save(fullfile(path2,'20211220laccmatarate_vuauvcollisioncount'),'vuauvcollisioncount');
    % save(fullfile(path2,'20211220laccmatarate_vusendauvrun_node_random'),'vusendauvrun_node_random');
    %
    % save(fullfile(path2,'20211220laccmatarate_vutotaldelay'),'vutotaldelay');
    % save(fullfile(path2,'20211220laccmatarate_vutotaldelay_plot'),'vutotaldelay_plot');
    % save(fullfile(path2,'20211220laccmatarate_vusenddelay'),'vusenddelay');
    % save(fullfile(path2,'20211220laccmatarate_vusenddelay_plot'),'vusenddelay_plot');

    %
    % save(fullfile(path2,'20211220laccmatarate_vudiscardnodecount_plot'),'vudiscardnodecount_plot');
    % save(fullfile(path2,'20211220laccmatarate_vucollisioncount_plot'),'vucollisioncount_plot');
    % save(fullfile(path2,'20211220laccmatarate_vuauvcollisioncount_plot'),'vuauvcollisioncount_plot');
    %
    % save(fullfile(path2,'20211220laccmatarate_vuthroughput_plot'),'vuthroughput_plot');
    % save(fullfile(path2,'20211220laccmatarate_vutotalthroughput_plot'),'vutotalthroughput_plot');
    %
    % save(fullfile(path2,'20211220laccmatarate_vuarrivaltime_plot'),'vuarrivaltime_plot');
    %
    % save(fullfile(path2,'20211220laccmatarate_vutotaldelay_plot_sum'),'vutotaldelay_plot_sum');
    % save(fullfile(path2,'20211220laccmatarate_vusenddelay_plot_sum'),'vusenddelay_plot_sum');

    %
    % save(fullfile(path2,'20211220laccmatarate_vmauvcollisioncount'),'vmauvcollisioncount');
    % save(fullfile(path2,'20211220laccmatarate_vmsendauvrun_node_random'),'vmsendauvrun_node_random');
    %
    % save(fullfile(path2,'20211220laccmatarate_vmtotaldelay'),'vmtotaldelay');
    % save(fullfile(path2,'20211220laccmatarate_vmtotaldelay_plot'),'vmtotaldelay_plot');
    % save(fullfile(path2,'20211220laccmatarate_vmsenddelay'),'vmsenddelay');
    % save(fullfile(path2,'20211220laccmatarate_vmsenddelay_plot'),'vmsenddelay_plot');


    %
    % save(fullfile(path2,'20211220laccmatarate_vmdiscardnodecount_plot'),'vmdiscardnodecount_plot');
    % save(fullfile(path2,'20211220laccmatarate_vmcollisioncount_plot'),'vmcollisioncount_plot');
    % save(fullfile(path2,'20211220laccmatarate_vmauvcollisioncount_plot'),'vmauvcollisioncount_plot');
    %
    % save(fullfile(path2,'20211220laccmatarate_vmthroughput_plot'),'vmthroughput_plot');
    % save(fullfile(path2,'20211220laccmatarate_vmtotalthroughput_plot'),'vmtotalthroughput_plot');
    %
    % save(fullfile(path2,'20211220laccmatarate_vmarrivaltime_plot'),'vmarrivaltime_plot');
    %
    % save(fullfile(path2,'20211220laccmatarate_vmtotaldelay_plot_sum'),'vmtotaldelay_plot_sum');
    % save(fullfile(path2,'20211220laccmatarate_vmsenddelay_plot_sum'),'vmsenddelay_plot_sum');

    % save(fullfile(path2,'20211220laccmatarate_vuvoi'),'vuvoi');
    % save(fullfile(path2,'20211220laccmatarate_vuvoi_collision'),'vuvoi_collision');
    % save(fullfile(path2,'20211220laccmatarate_vuvoi_plot'),'vuvoi_plot');
    % save(fullfile(path2,'20211220laccmatarate_vuvoi_collision_plot'),'vuvoi_collision_plot');
    % save(fullfile(path2,'20211220laccmatarate_vuvoi_plot_sum'),'vuvoi_plot_sum');
    % save(fullfile(path2,'20211220laccmatarate_vuvoi_collision_plot_sum'),'vuvoi_collision_plot_sum');
    %
    % save(fullfile(path2,'20211220laccmatarate_vmvoi'),'vmvoi');
    % save(fullfile(path2,'20211220laccmatarate_vmvoi_collision'),'vmvoi_collision');
    % save(fullfile(path2,'20211220laccmatarate_vmvoi_plot'),'vmvoi_plot');
    % save(fullfile(path2,'20211220laccmatarate_vmvoi_collision_plot'),'vmvoi_collision_plot');
    % save(fullfile(path2,'20211220laccmatarate_vmvoi_plot_sum'),'vmvoi_plot_sum');
    % save(fullfile(path2,'20211220laccmatarate_vmvoi_collision_plot_sum'),'vmvoi_collision_plot_sum');
    iteration

end


for iNA_plot=1:NA
    for ivutrafficload=1:length(vutrafficload)
        for ivmtrafficload=1:length(vmtrafficload)
            vutotaldelay_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vutotaldelay_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmtotaldelay_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmtotaldelay_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vusenddelay_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vusenddelay_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmsenddelay_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmsenddelay_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;



            vudiscardnodecount_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vudiscardnodecount_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmdiscardnodecount_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmdiscardnodecount_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vucollisioncount_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vucollisioncount_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vuauvcollisioncount_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vuauvcollisioncount_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmcollisioncount_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmcollisioncount_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmauvcollisioncount_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmauvcollisioncount_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vuthroughput_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vuthroughput_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmthroughput_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmthroughput_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vutotalthroughput_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vutotalthroughput_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmtotalthroughput_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmtotalthroughput_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;


            vuvoi_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vuvoi_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmvoi_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmvoi_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vuvoi_collision_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vuvoi_collision_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
            vmvoi_collision_laccm_average(iNA_plot,ivutrafficload,ivmtrafficload)=sum(vmvoi_collision_plot(iNA_plot,:,ivutrafficload,ivmtrafficload))/iteration_total;
        end
    end
end




%%delay
% save(fullfile(path2,'20211220laccmatarate_vutotaldelay_laccm_average'),'vutotaldelay_laccm_average');
% save(fullfile(path2,'20211220laccmatarate_vusenddelay_laccm_average'),'vusenddelay_laccm_average');
%
% save(fullfile(path2,'20211220laccmatarate_vudiscardnodecount_laccm_average'),'vudiscardnodecount_laccm_average');
% save(fullfile(path2,'20211220laccmatarate_vucollisioncount_laccm_average'),'vucollisioncount_laccm_average');
% save(fullfile(path2,'20211220laccmatarate_vuauvcollisioncount_laccm_average'),'vuauvcollisioncount_laccm_average');
%
% save(fullfile(path2,'20211220laccmatarate_vuthroughput_laccm_average'),'vuthroughput_laccm_average');
% save(fullfile(path2,'20211220laccmatarate_vutotalthroughput_laccm_average'),'vutotalthroughput_laccm_average');
%
% save(fullfile(path2,'20211220laccmatarate_vmtotaldelay_laccm_average'),'vmtotaldelay_laccm_average');
% save(fullfile(path2,'20211220laccmatarate_vmsenddelay_laccm_average'),'vmsenddelay_laccm_average');
%
% save(fullfile(path2,'20211220laccmatarate_vmdiscardnodecount_laccm_average'),'vmdiscardnodecount_laccm_average');
% save(fullfile(path2,'20211220laccmatarate_vmcollisioncount_laccm_average'),'vmcollisioncount_laccm_average');
% save(fullfile(path2,'20211220laccmatarate_vmauvcollisioncount_laccm_average'),'vmauvcollisioncount_laccm_average');
%
% save(fullfile(path2,'20211220laccmatarate_vmthroughput_laccm_average'),'vmthroughput_laccm_average');
% save(fullfile(path2,'20211220laccmatarate_vmtotalthroughput_laccm_average'),'vmtotalthroughput_laccm_average');
%
% save(fullfile(path2,'20211220laccmatarate_vuvoi_laccm_average'),'vuvoi_laccm_average');
% save(fullfile(path2,'20211220laccmatarate_vuvoi_collision_laccm_average'),'vuvoi_collision_laccm_average');

% save(fullfile(path2,'20211220laccmatarate_vmvoi_laccm_average'),'vmvoi_laccm_average');
% save(fullfile(path2,'20211220laccmatarate_vmvoi_collision_laccm_average'),'vmvoi_collision_laccm_average');

%%
clr(1,:)=[0 113 188]/255;
clr(2,:)=[216 82 24]/255;
clr(3,:)=[236 176 31]/255;
clr(4,:)=[125 46 141]/255;
clr(5,:)=[118 171 47]/255;
clr(6,:)=[76 189 237]/255;
clr(7,:)=[161 19 46]/255;
clr(8,:)=[0 0 0]/255;
clr(9,:)=[25 35 45]/255;
clr(10,:)=[100 100 100]/255;
clr(11,:)=[150 150 150]/255;


%% Network throughput
figure1=figure;
for iNA=1:NA
    for iteration=1:iteration_total
        for ivutrafficload=1:length(vutrafficload)
            vulambda=vutrafficload(ivutrafficload);
            for ivmtrafficload=1:length(vmtrafficload)
                vutotaltimearrival(iNA,iteration,ivutrafficload,ivmtrafficload)=floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vulambda)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);

            end
        end
    end
end
for iNA=1:NA
    for ivutrafficload=1:length(vutrafficload)
        for ivmtrafficload=1:length(vmtrafficload)
            vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)=sum(vutotaltimearrival(iNA,:,ivutrafficload,ivmtrafficload))/iteration_total;

            vuthroughput_laccm_average(iNA,ivutrafficload,ivmtrafficload)=(vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vudiscardnodecount_laccm_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vucollisioncount_laccm_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vuauvcollisioncount_laccm_average(iNA,ivutrafficload,ivmtrafficload)*vutotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload))*t_data/(len_xyz(iNA)*timeplus-timestart-timeend);

        end
    end
end

for iNA=1:NA
    for iteration=1:iteration_total
        for ivutrafficload=1:length(vutrafficload)

            for ivmtrafficload=1:length(vmtrafficload)
                vmlambda=vmtrafficload(ivmtrafficload);
                vmtotaltimearrival(iNA,iteration,ivutrafficload,ivmtrafficload)=floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vmlambda)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);

            end
        end
    end
end
for iNA=1:NA
    for ivutrafficload=1:length(vutrafficload)
        for ivmtrafficload=1:length(vmtrafficload)
            vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)=sum(vmtotaltimearrival(iNA,:,ivutrafficload,ivmtrafficload))/iteration_total;

            vmthroughput_laccm_average(iNA,ivutrafficload,ivmtrafficload)=(vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmdiscardnodecount_laccm_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmcollisioncount_laccm_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload)-vmauvcollisioncount_laccm_average(iNA,ivutrafficload,ivmtrafficload)*vmtotaltimearrival_average(iNA,ivutrafficload,ivmtrafficload))*t_data/(len_xyz(iNA)*timeplus-timestart-timeend);

        end
    end
end



for ivutrafficload=length(vutrafficload)%1:length(vmtrafficload)
    for i=1:length(vmtrafficload)
        A(i)=sum(vmthroughput_laccm_average(:,ivutrafficload,i)/t_data*L_packet);
    end
    plot(vmtrafficload(2:end),A(2:end).','s-b','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
    hold on
end
backColor = [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',0.1:0.1:1);

xlabel('\fontsize{16}Packet Arrival Rate for Normal Data (Packets/s) ');
ylabel('\fontsize{16}Network Throughput (Bits/s)');
legend1=legend('\fontsize{16}\fontname{Times New Roman}VBPS-I','Location','NorthEast');

%% VoI
figure2=figure;

for ivutrafficload=length(vutrafficload)%1:length(vmtrafficload)
    for i=1:length(vmtrafficload)
        A(i)=sum(vmvoi_collision_laccm_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A(2:end).','s-b','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
    hold on
end

backColor = [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',0.1:0.1:1);
xlabel('\fontsize{16}Data traffic loads for normal data (pkts/s)');
ylabel('\fontsize{16}Cumulative VoI ');
legend1=legend('\fontsize{16}\fontname{Times New Roman}VBPS-I','Location','NorthEast');

%% Network delay
figure3=figure;

for ivutrafficload=length(vutrafficload)%1:length(vmtrafficload)
    for i=1:length(vmtrafficload)
        A8(i)=sum(vmtotaldelay_laccm_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A8(2:end).','s-r','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
    hold on
end

backColor = [245 249 253]/255;
set(gca, 'color','none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xaxislocation','bottom');
set(gca,'xtick',0.1:0.1:1);
axis([-inf ,inf,16,40])
xlabel('\fontsize{16}Data traffic loads for normal data(pkts/s)');
ylabel('\fontsize{16} Average end-to-end delay (s)');

legend1=legend('\fontsize{16}\fontname{Times New Roman}VBPS-I','Location','NorthEast');

%% Collision Probability
figure4=figure;
for ivutrafficload=length(vutrafficload)%1:length(vmtrafficload)
    for i=1:length(vmtrafficload)
        A5(i)=sum(vmauvcollisioncount_laccm_average(:,ivutrafficload,i))+sum(vmcollisioncount_laccm_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A5(2:end).','s-b','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
    hold on
end
backColor = [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',0.1:0.1:1);
axis([-inf ,inf,0,0.18]);
xlabel('\fontsize{16}Packet Arrival Rate for Normal Data (Packets/s) ');
ylabel('\fontsize{16}Collision Probability  ');
legend1=legend('\fontsize{16}\fontname{Times New Roman}VBPS-I','Location','NorthEast');

%% Congestion ratio
figure5=figure;

for ivmtrafficload=length(vmtrafficload)%1:length(vmtrafficload)
    plot(vmtrafficload(2:end),sum(vudiscardnodecount_laccm_average(:,2:end,ivmtrafficload)),'s-r','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
    hold on
end


for ivutrafficload=length(vutrafficload)%1:length(vmtrafficload)
    for i=1:length(vmtrafficload)
        A(i)=sum(vmdiscardnodecount_laccm_average(:,ivutrafficload,i));
    end
    plot(vmtrafficload(2:end),A(2:end).','s--b','LineWidth',1.8,'MarkerSize',8,'color',clr(4,:));
    hold on
end

backColor = [245 249 253]/255;
set(gca, 'color', 'none');
grid on; set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'xtick',0.1:0.1:1);
axis([-inf ,inf,0.01,0.08]);
xlabel('\fontsize{16}Packet Arrival Rate (Packets/s)  ');
ylabel('\fontsize{16}Congestion Ratio');
legend1=legend('\fontsize{15}\fontname{Times New Roman}VBPS-I\_u','\fontsize{15}\fontname{Times New Roman}VBPS-I\_m','Location','NorthEast');

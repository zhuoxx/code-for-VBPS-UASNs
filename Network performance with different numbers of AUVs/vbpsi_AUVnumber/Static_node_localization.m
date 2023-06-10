function [AUV_listen,AUV_detected,AUV_in] =  Static_node_localization(Ax, Ay, Az,Ts, Nx, Ny, Nz,n,t_data,path1,N0,N1,N2,NA,N,Round,totaltimeslot,timeplus,path,node,len_xyz,dc)
%%Static node localization as illustrated in the paper


% By receiving information from a node, the next hop can be determined. 
% For each Autonomous Underwater Vehicle (AUV) located at a specific position, 
% the time at which the node sent the packet and the time at which the AUV 
% received the packet is recorded over several rounds. 
% The node's location can then be calculated. 
% Initially, the AUV needs to calculate the nodes within its range and 
% determine whether those nodes sent data during the time frame. 
% The static nodes have Ax, Ay, Az coordinates, Ts sending time, Nx, Ny, Nz information. 
% Each node sends data packets according to a time slot. 
% The AUV can determine the number of nodes in the cluster based on the time 
% difference between two consecutive data packets sent from a node. 
% By listening to the time difference between two nodes, the time slot length 
% can be calculated. 



%%---------------------Initialization--------------------------------------
timeslot_length =  n*t_data;
for iN =  N0+1:N0+N1+N2
    AUV_na(iN).round =  [];
    AUV_na(iN).tstime =  [];
    for iNA =  1:NA
        AUV(iNA).listen(iN) =  0;
        AUV(iNA).round =  [];
        AUV(iNA).tstime =  [];
        AUV_na(iN).count =  0;
        AUV(iNA).locationawareness(iN) =  0;
        AUV(iNA).timecostawareness(iN) =  0;
        AUV(iNA).hopawareness(iN) =  0;
    end
end
for itimeslot =  1:totaltimeslot
    AUV(itimeslot).in =  [];
end
AUV_in =  zeros(NA,totaltimeslot,N);
AUV_listen =  zeros(NA,totaltimeslot,N);
flag_itimeslot =  [];


%%---------Static node localization and topology construction--------------
for iNA =  1:NA
    for iN =  N0+1:N0+N1+N2
        for iround =  1:Round
            vuAUVtimebefore =  floor(Ts(iN,iround)/timeplus);
            vuAUVtimeplus =  mod(Ts(iN,iround),timeplus);
            itimeslot =  (iround-1)*path(iN).count+(node(iN).seq-1);
            flag_itimeslot =  [flag_itimeslot,itimeslot];
            if vuAUVtimebefore>0
                if vuAUVtimebefore+1<=  len_xyz(iNA)
                    Ax_time =  Ax(iNA,vuAUVtimebefore)+(Ax(iNA,vuAUVtimebefore+1)-Ax(iNA,vuAUVtimebefore))/timeplus*vuAUVtimeplus; %The location of AUV when static nodes send data packets节点发送信息的时候AUV的位置
                    Ay_time =  Ay(iNA,vuAUVtimebefore)+(Ay(iNA,vuAUVtimebefore+1)-Ay(iNA,vuAUVtimebefore))/timeplus*vuAUVtimeplus;
                    Az_time =  Az(iNA,vuAUVtimebefore)+(Az(iNA,vuAUVtimebefore+1)-Az(iNA,vuAUVtimebefore))/timeplus*vuAUVtimeplus;
                    
                    if (Ax_time-Nx(iN))^2+(Ay_time-Ny(iN))^2+(Az_time-Nz(iN))^2<=dc^2
                        AUV(iNA).listen(iN) =  AUV(iNA).listen(iN)+1;
                        AUV_listen(iNA,itimeslot,iN) =  1;%The listened static nodes in each time slot 
                        AUV(itimeslot).in =  union(AUV(itimeslot).in,iN);%The static nodes in the communication range of AUV
                        % if AUV(iNA).listen(iN)<=3
                        %  AUV(iNA).lasttime(iN)=Ts(iN,iround);
                        % end
                        AUV_na(iN).round =  [AUV_na(iN).round,iround]; 
                        AUV_na(iN).tstime =  [AUV_na(iN).tstime,Ts(iN,iround)];
                        
                        AUV(iNA).round =  [AUV(iNA).round,iround]; % The communication round of static nodes 
                        AUV(iNA).tstime =  [AUV(iNA).tstime,Ts(iN,iround)];%The transmission time of static nodes
                        
                    end
                end
            end
        end
        
    end
    for itimeslot =  1:totaltimeslot
        AUV_in(iNA,itimeslot,1:length(AUV(itimeslot).in)) =   AUV(itimeslot).in;%The static nodes in the communication range of AUV
    end
end


%%--------------Multi-AUV cooperation--------------------------------------
for iN =  N0+1:N0+N1+N2
    for itimeslot =  1:totaltimeslot
        for iNA1 =  1:NA
            if AUV_listen(iNA1,itimeslot,iN)==1
                for iNA2 =  1:NA
                    AUV_listen(iNA2,itimeslot,iN) =  1;
                end
            end
        end
    end
end
%%-----------The static nodes that located---------------------------------
AUV_detected =  zeros(NA,totaltimeslot,N);
for iNA =  1:NA
    for iN =  N0:N0+N1+N2
        for itimeslot =  1:totaltimeslot
            if sum(AUV_listen(iNA,1:itimeslot,iN))>=3
                AUV_detected(iNA,itimeslot,iN) =  1;
            end
        end
    end
end


%%-------------------Save and load the locations--------------------------
%Since the path of AUVs and transmission time and locations of static nodes 
% are predefined, we directly load the obtained localization of static
% nodes in the simulation.

% If you  have changed the path of AUV or the transmission time and 
% locations of static ndoes, you can save the new obtained localization of 
% static nodes.
% save(fullfile(path1,strcat('AUV_listen',num2str(NA))),'AUV_listen');
% save(fullfile(path1,strcat('AUV_detected',num2str(NA))),'AUV_detected');
% save(fullfile(path1,strcat('AUV_in',num2str(NA))),'AUV_in');

 AUV_listen = cell2mat(struct2cell(load(fullfile(path1,strcat('AUV_listen',num2str(NA))))));
AUV_detected = cell2mat(struct2cell(load(fullfile(path1,strcat('AUV_detected',num2str(NA))))));
AUV_in = cell2mat(struct2cell(load(fullfile(path1,strcat('AUV_in',num2str(NA))))));



end
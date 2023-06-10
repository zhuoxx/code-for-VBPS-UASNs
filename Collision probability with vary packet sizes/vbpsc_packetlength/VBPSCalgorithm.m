function [vuAUVrun_detected,vuarrivaltime,vuAUVrun_in,vuAUVrun_knownnodes,vuAUVrun_knownnodesin,vuAUVrun_unknownnodes,vuAUVrun_unknownnodesin,vuAUVtimebefore,vuAUVtimeplus,vutimeslot_number,vutimeslot_numberplus,vudiscardnodecount,vuknown_nodescount,vuknown_nodesincount,vuqueuearrivaltime,vuqueuecount,vuqueueindex,vuqueueout,vuqueuesendtime,vusendauvrun_time,vusendauvrun_timecase1,vusendauvrun_timecase2,vusendauvrun_timecase3,vusendauvrun_timecase4,vusenddelay,vutime_interval,vutotaldelay, vuunknown_nodescount,vuunknown_nodesincount,vuAUVrun_distancenode,vucollisioncount,vutotaltimearrival,vuvoi_choose,vuvoi_chooseindex,vuvoi_choosenode,vuvoitotaldelay_choose,vuvoicollisioncount,vuvoi,vuauvcollisioncount,vuvoisenddelay_choose,vuvoiauvcollisioncount_choose,vuvoi_collision]  =  VBPSCalgorithm(Nx,Ny,Nz,node,path,timeslot,len_xyz,dis_table,Ax,Ay,Az,AUV_in,AUV_detected,vulambda,NA,N,timeplus,n,t_data,N0,N1,N2,vs,vutrafficload,ivutrafficload,iteration,betau,vu0,alphau,Tu,L_packet,timestart,timeend)

%%Topology construction and transmission scheduling
%The data packets of AUV arrive with Poisson distribution. 
% When each data packets arrives, AUV will calculate the locations of 
% listened static nodes and topology. Then AUV calculate the sending time 
% of data packets and the next-hop node. Different from the VBPS-I
% algorithm, the VBPS-C algorithm considered the collision among AUVs

%%----------------------------data arrival time---------------------------- 
for iNA = 1:NA
    vutotaltimearrival(iNA,:)  =  floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vutrafficload)); %[51,51];%floor((len_xyz*timeplus-timestart-timeend)/time_interval);
end

%%----------------------------Initilization---------------------------- 
vutime_interval(:,:,iteration,ivutrafficload) = zeros(NA,max(max(vutotaltimearrival)));
vuAUVrun_detected = zeros(NA,max(max(vutotaltimearrival)),N);
vuAUVrun_in = zeros(NA,max(max(vutotaltimearrival)),N);
vuknown_nodescount = zeros(NA,max(max(vutotaltimearrival)));
vuunknown_nodescount = zeros(NA,max(max(vutotaltimearrival)));
vuAUVrun_knownnodes = zeros(NA,max(max(vutotaltimearrival)),N);
vuAUVrun_unknownnodes = zeros(NA,max(max(vutotaltimearrival)),N);
vuknown_nodesincount = zeros(NA,max(max(vutotaltimearrival)));
vuunknown_nodesincount = zeros(NA,max(max(vutotaltimearrival)));
vuAUVrun_knownnodesin = zeros(NA,max(max(vutotaltimearrival)),N);
vuAUVrun_unknownnodesin = zeros(NA,max(max(vutotaltimearrival)),N);
vusendauvrun_timecase1 = zeros(NA,max(max(vutotaltimearrival)),N);
vusendauvrun_timecase2 = zeros(NA,max(max(vutotaltimearrival)),N);
vusendauvrun_timecase3 = zeros(NA,max(max(vutotaltimearrival)),N);
vusendauvrun_timecase4 = zeros(NA,max(max(vutotaltimearrival)),N);
vusendauvrun_time = zeros(NA,max(max(vutotaltimearrival)),N);
vudiscardnodecount = zeros(NA,max(max(vutotaltimearrival)));
vuAUVrun_distancenode = zeros(NA,max(max(vutotaltimearrival)),N);
vuAx_packetarrivaltime = zeros(NA,max(max(vutotaltimearrival)));
vuAy_packetarrivaltime = zeros(NA,max(max(vutotaltimearrival)));
vuAz_packetarrivaltime = zeros(NA,max(max(vutotaltimearrival)));
vuqueuecount = ones(NA,max(max(vutotaltimearrival)));
vuarrivaltime = zeros(NA,max(max(vutotaltimearrival)))+timestart;
vuqueueout = ones(NA,max(max(vutotaltimearrival)));
vuqueuearrivaltime = zeros(NA,max(max(vutotaltimearrival)));
vuqueueindex = ones(NA,max(max(vutotaltimearrival)));
vuqueuesendtime = zeros(NA,max(max(vutotaltimearrival)),N);
vuAUVtimebefore = zeros(NA,max(max(vutotaltimearrival)));
vuAUVtimeplus = zeros(NA,max(max(vutotaltimearrival)));
vutimeslot_number = zeros(NA,max(max(vutotaltimearrival)));
vutimeslot_numberplus = zeros(NA,max(max(vutotaltimearrival)));

for iNA = 1:NA
    vutime_interval(iNA,1:vutotaltimearrival(iNA,ivutrafficload),iteration,ivutrafficload) = exprnd(1/vulambda,1,vutotaltimearrival(iNA,ivutrafficload));
    for iarrival = 2:vutotaltimearrival(iNA,ivutrafficload) %Each data packets arrival time  
        vuarrivaltime(iNA,iarrival) = vuarrivaltime(iNA,iarrival-1)+vutime_interval(iNA,iarrival,iteration,ivutrafficload);
        vuqueuecount(iNA,iarrival) = vuqueuecount(iNA,iarrival-1)+1;
        vuqueuearrivaltime(iNA,vuqueuecount(iNA,iarrival)) = vuarrivaltime(iNA,iarrival); 
        
        vuAUVtimebefore(iNA,iarrival) = floor(vuarrivaltime(iNA,iarrival)/timeplus);
        vuAUVtimeplus(iNA,iarrival) = mod(vuarrivaltime(iNA,iarrival),timeplus);
        vutimeslot_number(iNA,iarrival) = floor(vuarrivaltime(iNA,iarrival)/(n*t_data));
        vutimeslot_numberplus(iNA,iarrival) = mod(vuarrivaltime(iNA,iarrival),(n*t_data));
        
        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)-1
            vuAx_packetarrivaltime(iNA,iarrival) = Ax(iNA,vuAUVtimebefore(iNA,iarrival))+(Ax(iNA,vuAUVtimebefore(iNA,iarrival)+1)-Ax(iNA,vuAUVtimebefore(iNA,iarrival)))/timeplus*vuAUVtimeplus(iNA,iarrival);  
            vuAy_packetarrivaltime(iNA,iarrival) = Ay(iNA,vuAUVtimebefore(iNA,iarrival))+(Ay(iNA,vuAUVtimebefore(iNA,iarrival)+1)-Ay(iNA,vuAUVtimebefore(iNA,iarrival)))/timeplus*vuAUVtimeplus(iNA,iarrival);
            vuAz_packetarrivaltime(iNA,iarrival) = Az(iNA,vuAUVtimebefore(iNA,iarrival))+(Az(iNA,vuAUVtimebefore(iNA,iarrival)+1)-Az(iNA,vuAUVtimebefore(iNA,iarrival)))/timeplus*vuAUVtimeplus(iNA,iarrival);
            
            vuAUVrun_in(iNA,iarrival,:) = AUV_in(iNA,vutimeslot_number(iNA,iarrival),:);
            vuAUVrun_detected(iNA,iarrival,:) = AUV_detected(iNA,vutimeslot_number(iNA,iarrival),:);
            
            %Localization for static nodes and topology construction
            for i = 1:length(vuAUVrun_in(iNA,iarrival,:))
                if vuAUVrun_in(iNA,iarrival,i)~=  0 %If the static noeds is in the communication range of AUV
                    if vuAUVrun_detected(iNA,iarrival,vuAUVrun_in(iNA,iarrival,i))==  1%If the locations of static nodes are obtained
                        vuknown_nodesincount(iNA,iarrival) = vuknown_nodesincount(iNA,iarrival)+1;
                        vuAUVrun_knownnodesin(iNA,iarrival,vuknown_nodesincount(iNA,iarrival)) = vuAUVrun_in(iNA,iarrival,i);
                    else
                        vuunknown_nodesincount(iNA,iarrival) = vuunknown_nodesincount(iNA,iarrival)+1;
                        vuAUVrun_unknownnodesin(iNA,iarrival,vuunknown_nodesincount(iNA,iarrival)) = vuAUVrun_in(iNA,iarrival,i);
                    end
                else
                    break;
                end
            end
            
            for iN = N0+1:N0+N1+N2   %If the other static nodes are known
                if vuAUVrun_detected(iNA,iarrival,iN)==  1 && ismember(iN, vuAUVrun_in(iNA,iarrival,i))==  0
                    vuknown_nodescount(iNA,iarrival) = vuknown_nodescount(iNA,iarrival)+1;
                    vuAUVrun_knownnodes(iNA,iarrival,vuknown_nodescount(iNA,iarrival)) = iN;
                elseif vuAUVrun_detected(iNA,iarrival,iN)==  0 && ismember(iN, vuAUVrun_in(iNA,iarrival,i))==  0
                    vuunknown_nodescount(iNA,iarrival) = vuunknown_nodescount(iNA,iarrival)+1;
                    vuAUVrun_unknownnodes(iNA,iarrival,vuunknown_nodescount(iNA,iarrival)) = vuAUVrun_in(iNA,iarrival,i);
                end
            end
            
            %% Packets Transmission 
            if vuknown_nodesincount(iNA,iarrival)>=  1
                for i = 1:vuknown_nodesincount(iNA,iarrival)
                    
                    vuAUVrun_distancenode(iNA,iarrival,i) =  sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx(vuAUVrun_knownnodesin(iNA,iarrival,i)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny(vuAUVrun_knownnodesin(iNA,iarrival,i)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz(vuAUVrun_knownnodesin(iNA,iarrival,i)))^2);
                    %Packet Collision Constraint
                    % Case1: Receiving-Sending Collisions at  the next-hop node
                    vureceiveauvrun_timecase1 = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,i)/vs; 
                    if abs(vureceiveauvrun_timecase1+t_data-(vutimeslot_number(iNA,iarrival)+1)*n*t_data)<t_data && ismember(vuAUVrun_knownnodesin(iNA,iarrival,i), timeslot(vutimeslot_number(iNA,iarrival)+1).nodeset)==  1%The data packets is arrived in the current time slot but arrive in the next time slot 
                        vusendauvrun_timecase1(iNA,iarrival,i) = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival)+(n*t_data+t_data-vutimeslot_numberplus(iNA,iarrival)-vuAUVrun_distancenode(iNA,iarrival,i)/vs);
                    else
                        vusendauvrun_timecase1(iNA,iarrival,i) = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                    end

                    % Case 2: Receiving-Receiving Collisions at the  next-hop node
                    if vuAUVrun_knownnodesin(iNA,iarrival,i)<=  N0+N1   
                        for jn = 1:length(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).sonnode)  
                            vusendnoderun_temp = path(vuAUVrun_knownnodesin(iNA,iarrival,i)).sonnode(jn);
                            if ismember(vusendnoderun_temp, timeslot(vutimeslot_number(iNA,iarrival)).nodeset)==  1 
                                vureceivesnrun_timecase2 = vutimeslot_number(iNA,iarrival)*n*t_data+dis_table(vusendnoderun_temp,vuAUVrun_knownnodesin(iNA,iarrival,i))/vs; 
                                vureceiveauvrun_timecase2 = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,i)/vs; 
                                vusendauvrun_timecase2(iNA,iarrival,i) = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival)+max(vureceivesnrun_timecase2+t_data-vureceiveauvrun_timecase2,0); 
                            else
                                vusendauvrun_timecase2(iNA,iarrival,i) = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                            end
                        end
                    else
                        vusendauvrun_timecase2(iNA,iarrival,i) = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                    end
                    
                    % Case 3: Receiving-Sending Collisions at  the  interference  nodes
                    vusendauvrun_timecase3_max = 0;
                    for jn = 1:vuknown_nodesincount(iNA,iarrival)   
                        if jn~=  i
                            if ismember(vuAUVrun_knownnodesin(iNA,iarrival,jn), timeslot(vutimeslot_number(iNA,iarrival)+1).nodeset)==  1 
                                vuAUVrun_distancenode(iNA,iarrival,jn) = sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2);
                                vureceiveauvrun_timecase3 = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,jn)/vs; 
                                if vureceiveauvrun_timecase3>=  (vutimeslot_number(iNA,iarrival)+1)*n*t_data
                                    vusendauvrun_timecase3_temp = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival)+(n*t_data+t_data-vutimeslot_numberplus(iNA,iarrival)-vuAUVrun_distancenode(iNA,iarrival,i)/vs);
                                else
                                    vusendauvrun_timecase3_temp = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                                end
                            else
                                vusendauvrun_timecase3_temp = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                            end
                            if vusendauvrun_timecase3_temp>vusendauvrun_timecase3_max
                                vusendauvrun_timecase3_max = vusendauvrun_timecase3_temp;
                            end
                        else
                            vusendauvrun_timecase3_max = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                        end
                    end
                    vusendauvrun_timecase3(iNA,iarrival,i) = vusendauvrun_timecase3_max;
                    
                    % Case 4: Receiving-Receiving Collisions at the interference  nodes
                    vusendauvrun_timecase4_max = 0;
                    for jn = 1:vuknown_nodesincount(iNA,iarrival)   
                        if jn~=  i
                            if vuAUVrun_knownnodesin(iNA,iarrival,jn)<=  N0+N1  
                                for kn = 1:length(path(vuAUVrun_knownnodesin(iNA,iarrival,jn)).sonnode)  
                                    vusendnoderun_temp = path(vuAUVrun_knownnodesin(iNA,iarrival,jn)).sonnode(kn);
                                    if ismember(vusendnoderun_temp, timeslot(vutimeslot_number(iNA,iarrival)).nodeset)==  1 
                                        vuAUVrun_distancenode(iNA,iarrival,jn) = sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2);
                                        vureceivesnrun_timecase4 = vutimeslot_number(iNA,iarrival)*n*t_data+dis_table(vusendnoderun_temp,vuAUVrun_knownnodesin(iNA,iarrival,jn))/vs; 
                                        vureceiveauvrun_timecase4 = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,jn)/vs; 
                                        vusendauvrun_timecase4_temp = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival)+max(vureceivesnrun_timecase4+t_data-vureceiveauvrun_timecase4,0); 
                                        
                                        if vusendauvrun_timecase4_temp>vusendauvrun_timecase4_max
                                            vusendauvrun_timecase4_max = vusendauvrun_timecase4_temp;
                                        end
                                    else
                                        vusendauvrun_timecase4_max = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                                    end
                                end
                            else
                                vusendauvrun_timecase4_max = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                            end
                        else
                            vusendauvrun_timecase4_max = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                        end
                    end
                    vusendauvrun_timecase4(iNA,iarrival,i) = vusendauvrun_timecase4_max;
                    vusendauvrun_time(iNA,iarrival,i) = max([vusendauvrun_timecase1(iNA,iarrival,i),vusendauvrun_timecase2(iNA,iarrival,i),vusendauvrun_timecase3(iNA,iarrival,i),vusendauvrun_timecase4(iNA,iarrival,i)]);
                    
                    vuqueueout(iNA,iarrival) =  vuqueueout(iNA,iarrival-1)+1;
                    vuqueueindex(iNA,iarrival) = vuqueuecount(iNA,iarrival)-(vuqueuecount(iNA,iarrival-1)-vuqueueout(iNA,iarrival-1));
                    vuqueuesendtime(iNA,vuqueueindex(iNA,iarrival),i) = vusendauvrun_time(iNA,iarrival,i);
                    
                end
                
                
                 %% If there are no static nodes around the AUV, the transmission will be delayed 
            elseif vuknown_nodesincount(iNA,iarrival)==  0 && iarrival< vutotaltimearrival(iNA,ivutrafficload)
                
                vunexttimearrivaltimeslotnumber = floor(( vuarrivaltime(iNA,iarrival)+vutime_interval(iNA,iarrival+1,iteration,ivutrafficload))/(n*t_data));
                
                if vunexttimearrivaltimeslotnumber>vutimeslot_number(iNA,iarrival)
                    vuflag_knownnodesincount = 0;
                    for vutimeslot_number_temp = vutimeslot_number(iNA,iarrival)+1:1:vunexttimearrivaltimeslotnumber 
                        
                        
                       %Recalculate the time  
                        vutimeslot_number(iNA,iarrival) = vutimeslot_number_temp+1;
                        vutimeslot_numberplus(iNA,iarrival) = 0;
                        vuAUVtimebefore(iNA,iarrival) = floor(vutimeslot_number(iNA,iarrival)*n*t_data/timeplus);
                        vuAUVtimeplus(iNA,iarrival) = mod(vutimeslot_number(iNA,iarrival)*n*t_data,timeplus);
                        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)-1
                            vuAx_packetarrivaltime(iNA,iarrival) = Ax(iNA,vuAUVtimebefore(iNA,iarrival))+(Ax(iNA,vuAUVtimebefore(iNA,iarrival)+1)-Ax(iNA,vuAUVtimebefore(iNA,iarrival)))/timeplus*vuAUVtimeplus(iNA,iarrival);  
                            vuAy_packetarrivaltime(iNA,iarrival) = Ay(iNA,vuAUVtimebefore(iNA,iarrival))+(Ay(iNA,vuAUVtimebefore(iNA,iarrival)+1)-Ay(iNA,vuAUVtimebefore(iNA,iarrival)))/timeplus*vuAUVtimeplus(iNA,iarrival);
                            vuAz_packetarrivaltime(iNA,iarrival) = Az(iNA,vuAUVtimebefore(iNA,iarrival))+(Az(iNA,vuAUVtimebefore(iNA,iarrival)+1)-Az(iNA,vuAUVtimebefore(iNA,iarrival)))/timeplus*vuAUVtimeplus(iNA,iarrival);
                            
                            vuAUVrun_in(iNA,iarrival,:) = AUV_in(iNA,vutimeslot_number(iNA,iarrival),:);
                            vuAUVrun_detected(iNA,iarrival,:) = AUV_detected(iNA,vutimeslot_number(iNA,iarrival),:);
                            
                            %Recalculate the localization and topology
                            for i = 1:length(vuAUVrun_in(iNA,iarrival,:))
                                if vuAUVrun_in(iNA,iarrival,i)~=  0  
                                    if vuAUVrun_detected(iNA,iarrival,vuAUVrun_in(iNA,iarrival,i))==  1 
                                        vuknown_nodesincount(iNA,iarrival) = vuknown_nodesincount(iNA,iarrival)+1;
                                        vuAUVrun_knownnodesin(iNA,iarrival,vuknown_nodesincount(iNA,iarrival)) = vuAUVrun_in(iNA,iarrival,i);
                                    else
                                        vuunknown_nodesincount(iNA,iarrival) = vuunknown_nodesincount(iNA,iarrival)+1;
                                        vuAUVrun_unknownnodesin(iNA,iarrival,vuunknown_nodesincount(iNA,iarrival)) = vuAUVrun_in(iNA,iarrival,i);
                                    end
                                else
                                    break;
                                end
                            end
                            
                            for iN = N0+1:N0+N1+N2    
                                if vuAUVrun_detected(iNA,iarrival,iN)==  1 && ismember(iN, vuAUVrun_in(iNA,iarrival,i))==  0
                                    vuknown_nodescount(iNA,iarrival) = vuknown_nodescount(iNA,iarrival)+1;
                                    vuAUVrun_knownnodes(iNA,iarrival,vuknown_nodescount(iNA,iarrival)) = iN;
                                elseif vuAUVrun_detected(iNA,iarrival,iN)==  0 && ismember(iN, vuAUVrun_in(iNA,iarrival,i))==  0
                                    vuunknown_nodescount(iNA,iarrival) = vuunknown_nodescount(iNA,iarrival)+1;
                                    vuAUVrun_unknownnodes(iNA,iarrival,vuunknown_nodescount(iNA,iarrival)) = vuAUVrun_in(iNA,iarrival,i);
                                end
                            end
                            
                            if vuknown_nodesincount(iNA,iarrival)>=  1
                                
                                for i = 1:vuknown_nodesincount(iNA,iarrival)
                                    vuAUVrun_distancenode(iNA,iarrival,i) =  sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx(vuAUVrun_knownnodesin(iNA,iarrival,i)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny(vuAUVrun_knownnodesin(iNA,iarrival,i)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz(vuAUVrun_knownnodesin(iNA,iarrival,i)))^2);
                                    
                                    % Case1: Receiving-Sending Collisions at  the next-hop node
                                    vusendauvrun_timecase1(iNA,iarrival,i) = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                                    
                                    % Case 2: Receiving-Receiving Collisions at the  next-hop node
                                    if vuAUVrun_knownnodesin(iNA,iarrival,i)<=  N0+N1  
                                        for jn = 1:length(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).sonnode) 
                                            vusendnoderun_temp = path(vuAUVrun_knownnodesin(iNA,iarrival,i)).sonnode(jn);
                                            if ismember(vusendnoderun_temp, timeslot(vutimeslot_number(iNA,iarrival)).nodeset)==  1
                                                vureceivesnrun_timecase2 = vutimeslot_number(iNA,iarrival)*n*t_data+dis_table(vusendnoderun_temp,vuAUVrun_knownnodesin(iNA,iarrival,i))/vs;
                                                vureceiveauvrun_timecase2 = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,i)/vs;
                                                vusendauvrun_timecase2(iNA,iarrival,i) = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival)+max(vureceivesnrun_timecase2+t_data-vureceiveauvrun_timecase2,0);
                                            end
                                        end
                                    else
                                        vusendauvrun_timecase2(iNA,iarrival,i) = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                                    end
                                    
                                    % Case 3: Receiving-Sending Collisions at  the interference  nodes 
                                    vusendauvrun_timecase3(iNA,iarrival,i) = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                                    
                                    % Case 4: Receiving-Receiving Collisions at the interference  nodes
                                    vusendauvrun_timecase4_max = 0;
                                    for jn = 1:vuknown_nodesincount(iNA,iarrival)   
                                        if jn~=  i
                                            if vuAUVrun_knownnodesin(iNA,iarrival,jn)<=  N0+N1 
                                                for kn = 1:length(path(vuAUVrun_knownnodesin(iNA,iarrival,jn)).sonnode) 
                                                    vusendnoderun_temp = path(vuAUVrun_knownnodesin(iNA,iarrival,jn)).sonnode(kn);
                                                    if ismember(vusendnoderun_temp, timeslot(vutimeslot_number(iNA,iarrival)).nodeset)==  1
                                                        vuAUVrun_distancenode(iNA,iarrival,jn) = sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2);
                                                        vureceivesnrun_timecase4 = vutimeslot_number(iNA,iarrival)*n*t_data+dis_table(vusendnoderun_temp,vuAUVrun_knownnodesin(iNA,iarrival,jn))/vs;
                                                        vureceiveauvrun_timecase4 = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,jn)/vs;
                                                        vusendauvrun_timecase4_temp = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival)+max(vureceivesnrun_timecase4+t_data-vureceiveauvrun_timecase4,0);
                                                        
                                                        if vusendauvrun_timecase4_temp>vusendauvrun_timecase4_max
                                                            vusendauvrun_timecase4_max = vusendauvrun_timecase4_temp;
                                                        end
                                                    else
                                                        vusendauvrun_timecase4_max = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                                                    end
                                                end
                                            else
                                                vusendauvrun_timecase4_max = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                                            end
                                        else
                                            vusendauvrun_timecase4_max = vutimeslot_number(iNA,iarrival)*n*t_data+vutimeslot_numberplus(iNA,iarrival);
                                        end
                                    end
                                    vusendauvrun_timecase4(iNA,iarrival,i) = vusendauvrun_timecase4_max;
                                    vusendauvrun_time(iNA,iarrival,i) = max([vusendauvrun_timecase1(iNA,iarrival,i),vusendauvrun_timecase2(iNA,iarrival,i),vusendauvrun_timecase3(iNA,iarrival,i),vusendauvrun_timecase4(iNA,iarrival,i)]);
                                    
                                    vuqueueout(iNA,iarrival) =  vuqueueout(iNA,iarrival-1)+1;
                                    vuqueueindex(iNA,iarrival) = vuqueuecount(iNA,iarrival)-(vuqueuecount(iNA,iarrival-1)-vuqueueout(iNA,iarrival-1));
                                    vuqueuesendtime(iNA,vuqueueindex(iNA,iarrival),i) = vusendauvrun_time(iNA,iarrival,i);
                                    
                                end
                                vuflag_knownnodesincount = 1;
                                break;
                                
                            end
                        end
                    end
                    if vuflag_knownnodesincount==  0
                        vudiscardnodecount(iNA,iarrival) = 1;
                        vuqueueout(iNA,iarrival) =  vuqueueout(iNA,iarrival-1);
                        vuqueueindex(iNA,iarrival) = vuqueuecount(iNA,iarrival)-(vuqueuecount(iNA,iarrival-1)-vuqueueout(iNA,iarrival-1));
                        
                    end
                else
                    vudiscardnodecount(iNA,iarrival) = 1;
                    vuqueueout(iNA,iarrival) =  vuqueueout(iNA,iarrival-1);
                    vuqueueindex(iNA,iarrival) = vuqueuecount(iNA,iarrival)-(vuqueuecount(iNA,iarrival-1)-vuqueueout(iNA,iarrival-1));
                    
                end
                
            end
            
            
            
        end
    end
end

% Case 5: Receiving-Receiving Collisions due to multiple AUVs
vuauvcollisioncount1 = zeros(NA,max(max(vutotaltimearrival)),N);
vucollisioncount = zeros(NA,max(max(vutotaltimearrival)));

for iNA = 1:NA
    for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload) 
        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)
            
            
            for jNA = iNA+1:NA
                for jarrival = 1:vutotaltimearrival(jNA,ivutrafficload)
                    if vuAUVtimebefore(jNA,jarrival)>0 && vuAUVtimebefore(jNA,jarrival)<len_xyz(jNA)
                        AUVrun_distanceauv(iNA,iarrival,jNA,jarrival) = sqrt((vuAx_packetarrivaltime(iNA,iarrival)-vuAx_packetarrivaltime(jNA,jarrival))^2+(vuAy_packetarrivaltime(iNA,iarrival)-vuAy_packetarrivaltime(jNA,jarrival))^2+(vuAz_packetarrivaltime(iNA,iarrival)-vuAz_packetarrivaltime(jNA,jarrival))^2);
                        if AUVrun_distanceauv(iNA,iarrival,jNA,jarrival)<3500
                            if vuAUVtimebefore(iNA,iarrival)>=  vuAUVtimebefore(jNA,jarrival)
                                for i = 1:vuknown_nodesincount(iNA,iarrival)
                                    vusendauvrun_time(iNA,iarrival,i) = vusendauvrun_time(iNA,iarrival,i)+AUVrun_distanceauv(iNA,iarrival,jNA,jarrival)/vs;
                                    for j = 1:vuknown_nodesincount(jNA,jarrival)
                                        if vuAUVrun_knownnodesin(iNA,iarrival,i)==  vuAUVrun_knownnodesin(jNA,jarrival,j)
                                            if vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)+t_data>vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j) && vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)<vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)
                                                vusendauvrun_time(jNA,jarrival,j) = vusendauvrun_time(jNA,jarrival,j)+(vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)+t_data-vusendauvrun_time(jNA,jarrival,j)-vuAUVrun_distancenode(jNA,jarrival,j));
                                                vuqueuesendtime(jNA,vuqueueindex(jNA,jarrival),j) = vusendauvrun_time(jNA,jarrival,j);
                                                vuauvcollisioncount1(iNA,iarrival,i) = vuauvcollisioncount1(iNA,iarrival,i)+1;
                                                vuauvcollisioncount1(jNA,jarrival,j) = vuauvcollisioncount1(jNA,jarrival,j)+1;
                                                %                                         vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                                                %                                         vucollisioncount(jNA,jarrival) = vucollisioncount(jNA,jarrival)+1;
                                                
                                            elseif vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)<vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)+t_data && vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)>vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)
                                                vusendauvrun_time(jNA,jarrival,j) = vusendauvrun_time(jNA,jarrival,j)+(vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)+t_data-(vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)));
                                                vuqueuesendtime(jNA,vuqueueindex(jNA,jarrival),j) = vusendauvrun_time(jNA,jarrival,j);
                                                vuauvcollisioncount1(iNA,iarrival,i) = vuauvcollisioncount1(iNA,iarrival,i)+1;
                                                vuauvcollisioncount1(jNA,jarrival,j) = vuauvcollisioncount1(jNA,jarrival,j)+1;
                                                %                                         vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                                                %                                         vucollisioncount(jNA,jarrival) = vucollisioncount(jNA,jarrival)+1;
                                            end
                                        end
                                    end
                                end
                            elseif  vuAUVtimebefore(iNA,iarrival)< vuAUVtimebefore(jNA,jarrival) 
                                for j = 1:vuknown_nodesincount(jNA,jarrival)
                                    vusendauvrun_time(jNA,jarrival,j) = vusendauvrun_time(jNA,jarrival,j)+AUVrun_distanceauv(iNA,iarrival,jNA,jarrival)/vs;
                                    for i = 1:vuknown_nodesincount(iNA,iarrival)
                                        if vuAUVrun_knownnodesin(iNA,iarrival,i)==  vuAUVrun_knownnodesin(jNA,jarrival,j)
                                            if vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)+t_data>vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j) && vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)<vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)
                                               vusendauvrun_time(iNA,iarrival,i) = vusendauvrun_time(iNA,iarrival,i)+(vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)+t_data)-(vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i));
                                                vuqueuesendtime(iNA,vuqueueindex(iNA,iarrival),i) = vusendauvrun_time(iNA,iarrival,i);
                                                vuauvcollisioncount1(iNA,iarrival,i) = vuauvcollisioncount1(iNA,iarrival,i)+1;
                                                vuauvcollisioncount1(jNA,jarrival,j) = vuauvcollisioncount1(jNA,jarrival,j)+1; 
                                                %                                         vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                                                %                                         vucollisioncount(jNA,jarrival) = vucollisioncount(jNA,jarrival)+1;
                                                
                                            elseif vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)<vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)+t_data && vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)>vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)
                                                vusendauvrun_time(iNA,iarrival,i) = vusendauvrun_time(iNA,iarrival,i)+vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)+t_data-(vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i));
                                                vuqueuesendtime(iNA,vuqueueindex(iNA,iarrival),i) = vusendauvrun_time(iNA,iarrival,i);
                                                vuauvcollisioncount1(iNA,iarrival,i) = vuauvcollisioncount1(iNA,iarrival,i)+1;
                                                vuauvcollisioncount1(jNA,jarrival,j) = vuauvcollisioncount1(jNA,jarrival,j)+1;
                                                %                                         vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                                                %                                         vucollisioncount(jNA,jarrival) = vucollisioncount(jNA,jarrival)+1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                         
                     end
                    
                end
            end
        end
    end
end

%Calculate the collision among multiple AUVs
vuauvcollisioncount = zeros(NA,max(max(vutotaltimearrival)),N);
for iNA = 1:NA
    for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload) %Arrival time
        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)
            
            
            for jNA = iNA+1:NA
                for jarrival = 1:vutotaltimearrival(jNA,ivutrafficload)
                    if vuAUVtimebefore(jNA,jarrival)>0 && vuAUVtimebefore(jNA,jarrival)<len_xyz(iNA)
                     
                        
                        for i = 1:vuknown_nodesincount(iNA,iarrival)
                            for j = 1:vuknown_nodesincount(jNA,iarrival)
                                if vuAUVrun_knownnodesin(iNA,iarrival,i)==  vuAUVrun_knownnodesin(jNA,iarrival,j)
                                    if  vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)/vs+t_data>vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)/vs && vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)/vs<vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)/vs
                                        vuauvcollisioncount(iNA,iarrival,i) = vuauvcollisioncount(iNA,iarrival,i)+1;
                                        vuauvcollisioncount(jNA,jarrival,j) = vuauvcollisioncount(jNA,jarrival,j)+1;
                                     elseif vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)/vs<vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)/vs+t_data && vusendauvrun_time(iNA,iarrival,i)+vuAUVrun_distancenode(iNA,iarrival,i)/vs>vusendauvrun_time(jNA,jarrival,j)+vuAUVrun_distancenode(jNA,jarrival,j)/vs
                                         vuauvcollisioncount(iNA,iarrival,i) = vuauvcollisioncount(iNA,iarrival,i)+1;
                                        vuauvcollisioncount(jNA,jarrival,j) = vuauvcollisioncount(jNA,jarrival,j)+1;
                                     end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
%%  Calculate the network throughput, network delay, VoI, congestion ratio, and collision probability
vuqueuesendtimebefore = zeros(NA,max(max(vutotaltimearrival)),N);
vuqueuesendtimeplus = zeros(NA,max(max(vutotaltimearrival)),N);
vuAx_packetsendtime = zeros(NA,max(max(vutotaltimearrival)),N);
vuAy_packetsendtime = zeros(NA,max(max(vutotaltimearrival)),N);
vuAz_packetsendtime = zeros(NA,max(max(vutotaltimearrival)),N);
vuAUVrun_sendtimedistancenode = zeros(NA,max(max(vutotaltimearrival)),N);
 
for iNA = 1:NA
    
        for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload)  
            if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)
                for i = 1:vuknown_nodesincount(iNA,iarrival)
                    if vuqueuesendtime(iNA,iarrival,i)>0
                        vuqueuesendtimebefore(iNA,iarrival,i) = floor(vuqueuesendtime(iNA,iarrival,i)/timeplus);
                        vuqueuesendtimeplus(iNA,iarrival,i) = mod(vuqueuesendtime(iNA,iarrival),timeplus);    
                        
                        if vuqueuesendtime(iNA,iarrival,i)>0 && vuqueuesendtime(iNA,iarrival,i)<len_xyz(iNA)*timeplus-timestart-timeend
                            vuAx_packetsendtime(iNA,iarrival,i) = Ax(iNA,vuqueuesendtimebefore(iNA,iarrival,i))+(Ax(iNA,vuqueuesendtimebefore(iNA,iarrival,i)+1)-Ax(iNA,vuqueuesendtimebefore(iNA,iarrival,i)))/timeplus*vuqueuesendtimeplus(iNA,iarrival,i); %节点发�?�信息的时�?�AUV的位�?
                            vuAy_packetsendtime(iNA,iarrival,i) = Ay(iNA,vuqueuesendtimebefore(iNA,iarrival,i))+(Ay(iNA,vuqueuesendtimebefore(iNA,iarrival,i)+1)-Ay(iNA,vuqueuesendtimebefore(iNA,iarrival,i)))/timeplus*vuqueuesendtimeplus(iNA,iarrival,i);
                            vuAz_packetsendtime(iNA,iarrival,i) = Az(iNA,vuqueuesendtimebefore(iNA,iarrival,i))+(Az(iNA,vuqueuesendtimebefore(iNA,iarrival,i)+1)-Az(iNA,vuqueuesendtimebefore(iNA,iarrival,i)))/timeplus*vuqueuesendtimeplus(iNA,iarrival,i);
                            vuAUVrun_sendtimedistancenode(iNA,iarrival,i) =  sqrt((vuAx_packetsendtime(iNA,iarrival,i)-Nx(vuAUVrun_knownnodesin(iNA,iarrival,i)))^2+(vuAy_packetsendtime(iNA,iarrival,i)-Ny(vuAUVrun_knownnodesin(iNA,iarrival,i)))^2+(vuAz_packetsendtime(iNA,iarrival,i)-Nz(vuAUVrun_knownnodesin(iNA,iarrival,i)))^2);
                            
                        end
                    end
                end
            end
        end
    
end

%Calculate delay
vutotaldelay = zeros(NA,max(max(vutotaltimearrival)),N);
packeterrorate = zeros(NA,max(max(vutotaltimearrival)),N);
vusenddelay = zeros(NA,max(max(vutotaltimearrival)),N);
for iNA = 1:NA
    for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload)  
        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)
            
            for i = 1:vuknown_nodesincount(iNA,iarrival)
                if vuqueuesendtime(iNA,iarrival,i)>0 && vuqueuesendtime(iNA,iarrival,i)<len_xyz(iNA)*timeplus-timestart-timeend-1
                    vusenddelay(iNA,iarrival,i) = vuqueuesendtime(iNA,iarrival,i)-vuqueuearrivaltime(iNA,iarrival);
                    packeterrorate(iNA,iarrival,i) = packeterroratefunction(vuAUVrun_sendtimedistancenode(iNA,iarrival,i),L_packet);  
                    if vuAUVrun_knownnodesin(iNA,iarrival,i)<=  N0  %The sink node
                        vutotaldelay(iNA,iarrival,i) =  vusenddelay(iNA,iarrival,i);
                        packeterrorate(iNA,iarrival,i) = packeterrorate(iNA,iarrival,i);
                    elseif vuAUVrun_knownnodesin(iNA,iarrival,i)>N0 && vuAUVrun_knownnodesin(iNA,iarrival,i)<=  N0+N1  %The relay node
                        vutimeslot_temp = floor(vuqueuesendtime(iNA,iarrival,i)/(n*t_data));
                        vuround_temp = floor(vutimeslot_temp/path(vuAUVrun_knownnodesin(iNA,iarrival,i)).count);
                        vuroundplus_temp = mod(vutimeslot_temp,path(vuAUVrun_knownnodesin(iNA,iarrival,i)).count);
                        if node(vuAUVrun_knownnodesin(iNA,iarrival,i)).seq < vuroundplus_temp
                            vutotaldelay(iNA,iarrival,i) = vuqueuesendtime(iNA,iarrival,i)+(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).count-vuroundplus_temp)*n*t_data+(node(vuAUVrun_knownnodesin(iNA,iarrival,i)).seq-1)*n*t_data+dis_table(vuAUVrun_knownnodesin(iNA,iarrival,i),path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop)/vs-vuqueuearrivaltime(iNA,iarrival);
                        elseif node(vuAUVrun_knownnodesin(iNA,iarrival,i)).seq >=   vuroundplus_temp
                            vutotaldelay(iNA,iarrival,i) = vuqueuesendtime(iNA,iarrival,i)+(node(vuAUVrun_knownnodesin(iNA,iarrival,i)).seq-vuroundplus_temp)*n*t_data+dis_table(vuAUVrun_knownnodesin(iNA,iarrival,i),path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop)/vs-vuqueuearrivaltime(iNA,iarrival);
                        end
                        packeterrorate(iNA,iarrival,i) = packeterrorate(iNA,iarrival,i)*packeterroratefunction(dis_table(vuAUVrun_knownnodesin(iNA,iarrival,i),path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop),L_packet);  
                    elseif vuAUVrun_knownnodesin(iNA,iarrival,i)>N0+N1    %The sensor node
                        vutimeslot_temp = floor(vuqueuesendtime(iNA,iarrival,i)/(n*t_data));
                        vuround_temp = floor(vutimeslot_temp/path(vuAUVrun_knownnodesin(iNA,iarrival,i)).count);
                        vuroundplus_temp = mod(vutimeslot_temp,path(vuAUVrun_knownnodesin(iNA,iarrival,i)).count);
                        packeterrorate(iNA,iarrival,i) = packeterrorate(iNA,iarrival,i)*packeterroratefunction(dis_table(vuAUVrun_knownnodesin(iNA,iarrival,i),path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop),L_packet)*packeterroratefunction(dis_table(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop,path(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).hop),L_packet);  
                        if node(vuAUVrun_knownnodesin(iNA,iarrival,i)).seq <=   vuroundplus_temp
                            vuhop1time = vuqueuesendtime(iNA,iarrival,i)+(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).count-vuroundplus_temp)*n*t_data+(node(vuAUVrun_knownnodesin(iNA,iarrival,i)).seq-1)*n*t_data+dis_table(vuAUVrun_knownnodesin(iNA,iarrival,i),path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop)/vs;
                            vuhop2timeslot = floor(vuhop1time/(n*t_data));
                            vuhop2round = floor(vuhop2timeslot/path(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).count);
                            vuhop2roundplus = mod(vuhop2timeslot,path(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).count);
                            if node(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).seq <=   vuhop2roundplus
                                vuhop2time = (path(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).count-vuhop2roundplus)*n*t_data+(node(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).seq-1)*n*t_data+dis_table(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop,path(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).hop)/vs;
                                vutotaldelay(iNA,iarrival,i) = vuhop1time+vuhop2time-vuqueuearrivaltime(iNA,iarrival);
                            elseif node(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).seq > vuhop2roundplus
                                vuhop2time = (node(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).seq-vuhop2roundplus)*n*t_data+dis_table(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop,path(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).hop)/vs;
                                vutotaldelay(iNA,iarrival,i) = vuhop1time+vuhop2time-vuqueuearrivaltime(iNA,iarrival);
                            end
                        elseif node(vuAUVrun_knownnodesin(iNA,iarrival,i)).seq > vuroundplus_temp
                            vuhop1time = vuqueuesendtime(iNA,iarrival,i)+(node(vuAUVrun_knownnodesin(iNA,iarrival,i)).seq-vuroundplus_temp)*n*t_data+dis_table(vuAUVrun_knownnodesin(iNA,iarrival,i),path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop)/vs;
                            vuhop2timeslot = floor(vuhop1time/(n*t_data));
                            vuhop2round = floor(vuhop2timeslot/path(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).count);
                            vuhop2roundplus = mod(vuhop2timeslot,path(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).count);
                            if node(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).seq <=   vuhop2roundplus
                                vuhop2time = (path(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).count-vuhop2roundplus)*n*t_data+(node(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).seq-1)*n*t_data+dis_table(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop,path(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).hop)/vs;
                                vutotaldelay(iNA,iarrival,i) = vuhop1time+vuhop2time-vuqueuearrivaltime(iNA,iarrival);
                            elseif node(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).seq > vuhop2roundplus
                                vuhop2time = (node(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).seq-vuhop2roundplus)*n*t_data+dis_table(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop,path(path(vuAUVrun_knownnodesin(iNA,iarrival,i)).hop).hop)/vs;
                                vutotaldelay(iNA,iarrival,i) = vuhop1time+vuhop2time-vuqueuearrivaltime(iNA,iarrival);
                            end
                        end
                    end
                end
            end
        end
    end
end


 
%% Choose the sending time and next-hop nodes to maximize value of information 
vut = zeros(NA,max(max(vutotaltimearrival)),N);
vuvoi = zeros(NA,max(max(vutotaltimearrival)),N);
for iNA = 1:NA
    for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload)  
        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)
             for i = 1:vuknown_nodesincount(iNA,iarrival)
                 if vuqueuesendtime(iNA,iarrival,i)>0 && vuqueuesendtime(iNA,iarrival,i)<len_xyz(iNA)*timeplus-timestart-timeend   %%%zxx1209
                    vusenddelay(iNA,iarrival,i) = vuqueuesendtime(iNA,iarrival,i)-vuqueuearrivaltime(iNA,iarrival);
                    vut(iNA,iarrival,i)  =  valueofinformation(vusenddelay(iNA,iarrival,i),betau,vu0,alphau,Tu);
                    vuvoi(iNA,iarrival,i)  =  vut(iNA,iarrival,i)*(1-packeterrorate(iNA,iarrival,i));
                end
            end
        end
    end
end

% Choose the maximum VoI
vuvoi_choose = zeros(NA,max(max(vutotaltimearrival)));
vuvoi_choosenode = zeros(NA,max(max(vutotaltimearrival)));
vuvoi_chooseindex = zeros(NA,max(max(vutotaltimearrival)));
vuvoitotaldelay_choose = 10000000000*ones(NA,max(max(vutotaltimearrival)));
vuvoisenddelay_choose = 10000000000*ones(NA,max(max(vutotaltimearrival)));
vuvoiauvcollisioncount_choose = zeros(NA,max(max(vutotaltimearrival))); 
for iNA = 1:NA
    for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload)  
       if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)
        for i = 1:vuknown_nodesincount(iNA,iarrival)
             if vuqueuesendtime(iNA,iarrival,i)>0 && vuqueuesendtime(iNA,iarrival,i)<len_xyz(iNA)*timeplus-timestart-timeend    
              
            if vuvoi(iNA,iarrival,i)>vuvoi_choose(iNA,iarrival)
                vuvoi_choose(iNA,iarrival) = vuvoi(iNA,iarrival,i);
                vuvoi_chooseindex(iNA,iarrival) = i;
                vuvoi_choosenode(iNA,iarrival) = vuAUVrun_knownnodesin(iNA,iarrival,i);
                 vuvoitotaldelay_choose(iNA,iarrival) = vutotaldelay(iNA,iarrival,i);
                vuvoisenddelay_choose(iNA,iarrival) = vusenddelay(iNA,iarrival,i);  
                vuvoiauvcollisioncount_choose(iNA,iarrival) = vuauvcollisioncount(iNA,iarrival,i);  
            end
             end
        end
       end
    end
end

 

 

            


%%  The collision for the unknown static nodes  voi

vuvoicollisioncount = zeros(NA,max(max(vutotaltimearrival)));
for iNA = 1:NA
    for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload)  
 
        % Case1: Receiving-Sending Collisions at  the next-hop node
        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)-1 &&vuvoi_chooseindex(iNA,iarrival)~=  0
            
            for jn = 1:vuunknown_nodesincount(iNA,iarrival)
                vuAUVrun_undistancenode(iNA,iarrival,jn) =  sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx(vuAUVrun_unknownnodesin(iNA,iarrival,jn)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny(vuAUVrun_unknownnodesin(iNA,iarrival,jn)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz(vuAUVrun_unknownnodesin(iNA,iarrival,jn)))^2);
                
                % Case1: Receiving-Sending Collisions at  the next-hop node
 
                vuunreceiveauvrun_timecase1 = vusendauvrun_time(iNA,iarrival,vuvoi_chooseindex(iNA,iarrival))+vuAUVrun_undistancenode(iNA,iarrival,jn)/vs; 
                if abs(vuunreceiveauvrun_timecase1+t_data-(vutimeslot_number(iNA,iarrival)+1)*n*t_data)<t_data && ismember(vuAUVrun_unknownnodesin(iNA,iarrival,jn), timeslot(vutimeslot_number(iNA,iarrival)+1).nodeset)==  1%当前slot数据包到达，但是发�?�之后，节点可能在下�?个slot接收�?
                    vuvoicollisioncount(iNA,iarrival) = 1;
                    break;
                end
                
                  % Case 2: Receiving-Receiving Collisions at the  next-hop
                % node  need not consider

                % Case 3: Receiving-Sending Collisions at  the interference
                % nodes  need not consider
                
                % Case 4: Receiving-Receiving Collisions at the interference  nodes
 
                if vuAUVrun_unknownnodesin(iNA,iarrival,jn)<=  N0+N1  
                    for sn = 1:length(path(vuAUVrun_unknownnodesin(iNA,iarrival,jn)).sonnode) 
                        vusendnoderun_temp = path(vuAUVrun_unknownnodesin(iNA,iarrival,jn)).sonnode(sn);
                        if ismember(vusendnoderun_temp, timeslot(vutimeslot_number(iNA,iarrival)).nodeset)==  1 
                            vuunreceivesnrun_timecase2 = vutimeslot_number(iNA,iarrival)*n*t_data+dis_table(vusendnoderun_temp,vuAUVrun_unknownnodesin(iNA,iarrival,jn))/vs; 
                            vuunreceiveauvrun_timecase2 = vusendauvrun_time(iNA,iarrival,vuvoi_chooseindex(iNA,iarrival))+vuAUVrun_distancenode(iNA,iarrival,jn)/vs; 
                            if abs(vuunreceivesnrun_timecase2-vuunreceiveauvrun_timecase2)<t_data
                                vuvoicollisioncount(iNA,iarrival) = 1;
                                break;
                            end
                            
                        end
                    end
                end
            end
        end
    end
end

vuvoi_collision = vuvoi_choose;
for iNA = 1:NA
    for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload) 
        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)-1 &&vuvoi_choosenode(iNA,iarrival)~=  0
            
            if vuvoicollisioncount(iNA,iarrival)~=  0 || vuvoiauvcollisioncount_choose(iNA,iarrival)~=  0
                vuvoi_collision(iNA,iarrival)  =  0;
            end
        end
    end
end

 
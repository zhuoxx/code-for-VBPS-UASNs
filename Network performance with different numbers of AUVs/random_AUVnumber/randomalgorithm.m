%function [vuAUVrun_detected,vuarrivaltime,vuAUVrun_in,vuAUVrun_knownnodes,vuAUVrun_knownnodesin,vuAUVrun_unknownnodes,vuAUVrun_unknownnodesin,vuAUVtimebefore(iNA,iarrival),vuAUVtimeplus(iNA,iarrival),vudiscardnodecount,vuknown_nodescount,vuknown_nodesincount,vuqueuearrivaltime,vuqueuecount,vuqueueindex,vuqueueout,vuqueuesendtime,vusendauvrun_time,vusendauvrun_timecase1,vusendauvrun_timecase2,vusendauvrun_timecase3,vusendauvrun_timecase4,vusenddelay,vutime_interval,vutotaldelay,vutotaldelay_choose, vutotaldelay_chooseindex, vutotaldelay_choosenode, vutotaltimearrival, vuunknown_nodescount, vuunknown_nodesincount ]  =  delaycomputing(Nx,Ny,Nz,node,path,sequence,timeslot,TM,IM,Ts,len_xyz,Local_N,dis_table,delay_table,Ax,Ay,Az,AUV_in, AUV_listen,AUV_detected,vulambda)
function [vuAUVrun_detected,vuarrivaltime,vuAUVrun_in,vuAUVrun_knownnodes,vuAUVrun_knownnodesin,vuAUVrun_unknownnodes,vuAUVrun_unknownnodesin,vuAUVtimebefore,vuAUVtimeplus,vutimeslot_number,vutimeslot_numberplus,vudiscardnodecount,vuknown_nodescount,vuknown_nodesincount,vuqueuearrivaltime,vuqueuecount,vuqueueindex,vuqueueout,vuqueuesendtime,vusendauvrun_time_random,vusenddelay,vutime_interval,vutotaldelay, vuunknown_nodescount,vuunknown_nodesincount,vuAUVrun_distancenode,vucollisioncount,vutotaltimearrival,vuvoi,vusendauvrun_node_random,vuauvcollisioncount,vuvoi_collision,vupacketerrorate]  =  randomalgorithm(Nx,Ny,Nz,node,path,timeslot,len_xyz,dis_table,Ax,Ay,Az,AUV_in,AUV_detected,vulambda,NA,N,timeplus,n,t_data,N0,N1,N2,vs,vutrafficload,ivutrafficload,iteration,betau,vu0,alphau,Tu,L_packet,timestart,timeend)
%%Topology construction and transmission scheduling
%The data packets of AUV arrive with Poisson distribution. 
% When each data packets arrives, AUV will send packet at the random delayed time. 

%%----------------------------data arrival time---------------------------- 
for iNA = 1:NA
    vutotaltimearrival(iNA,:) = floor((len_xyz(iNA).*timeplus-timestart-timeend)./(1./vutrafficload));  
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

vusendauvrun_time_random = zeros(NA,max(max(vutotaltimearrival)));
vusendauvrun_node_random = zeros(NA,max(max(vutotaltimearrival)));
vusendauvrun_timeslot_number = zeros(NA,max(max(vutotaltimearrival)));
vusendauvrun_timeslot_numberplus = zeros(NA,max(max(vutotaltimearrival)));
vureceiveauvrun_time = zeros(NA,max(max(vutotaltimearrival)));
vureceiveauvrun_timecase1timeslot_number = zeros(NA,max(max(vutotaltimearrival)));
vureceiveauvrun_timecase1timeslot_numberplus = zeros(NA,max(max(vutotaltimearrival)));

vuAUVrun_distancenode = zeros(NA,max(max(vutotaltimearrival)),N);

vudiscardnodecount = zeros(NA,max(max(vutotaltimearrival)));
vutotaldiscardnodecount = zeros(NA,max(max(vutotaltimearrival))); 

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

vucollisioncount = zeros(NA,max(max(vutotaltimearrival)));
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
                if vuAUVrun_in(iNA,iarrival,i)~= 0 %Whether the static noeds is in the communication range of AUV
                    if vuAUVrun_detected(iNA,iarrival,vuAUVrun_in(iNA,iarrival,i))==1%If the locations of static nodes are obtained
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
            
            for iN = N0+1:N0+N1+N2    %Whether the other static nodes are known
                if vuAUVrun_detected(iNA,iarrival,iN)==1 && ismember(iN, vuAUVrun_in(iNA,iarrival,i))==0
                    vuknown_nodescount(iNA,iarrival) = vuknown_nodescount(iNA,iarrival)+1;
                    vuAUVrun_knownnodes(iNA,iarrival,vuknown_nodescount(iNA,iarrival)) = iN;
                elseif vuAUVrun_detected(iNA,iarrival,iN)==0 && ismember(iN, vuAUVrun_in(iNA,iarrival,i))==0
                    vuunknown_nodescount(iNA,iarrival) = vuunknown_nodescount(iNA,iarrival)+1;
                    vuAUVrun_unknownnodes(iNA,iarrival,vuunknown_nodescount(iNA,iarrival)) = vuAUVrun_in(iNA,iarrival,i);
                end
            end
            
            
           %% Packets Transmission 
            if vuknown_nodesincount(iNA,iarrival)>= 1
                %select a random node to send at random time
                vusendauvrun_time_random(iNA,iarrival) = vuarrivaltime(iNA,iarrival) + (5*n*t_data) * rand(1,1);
                vusendauvrun_node_random(iNA,iarrival) = vuAUVrun_knownnodesin(iNA,iarrival,randi([1,vuknown_nodesincount(iNA,iarrival)],1,1));
                
                vusendauvrun_timeslot_number(iNA,iarrival) = floor(vusendauvrun_time_random(iNA,iarrival)/(n*t_data));
                vusendauvrun_timeslot_numberplus(iNA,iarrival) = mod(vusendauvrun_time_random(iNA,iarrival),(n*t_data));
                
                vuAUVrun_distancenode(iNA,iarrival,vusendauvrun_node_random(iNA,iarrival)) =  sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx( vusendauvrun_node_random(iNA,iarrival)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny( vusendauvrun_node_random(iNA,iarrival)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz( vusendauvrun_node_random(iNA,iarrival)))^2);
                
                vureceiveauvrun_time(iNA,iarrival) = vusendauvrun_time_random(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,vusendauvrun_node_random(iNA,iarrival))/vs; 
                
                vureceiveauvrun_timecase1timeslot_number(iNA,iarrival) = floor(vureceiveauvrun_time(iNA,iarrival)/(n*t_data));
                vureceiveauvrun_timecase1timeslot_numberplus(iNA,iarrival) = mod(vureceiveauvrun_time(iNA,iarrival),(n*t_data));
                
                %Calculate the packet collision constraint
                
                % Case1: Receiving-Sending Collisions at  the next-hop node
                
                
                if  ismember(vusendauvrun_node_random(iNA,iarrival), timeslot(vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)).nodeset)==1%The data packet is sent by AUV at the current slot, but may be receieved by static nodes at the next time slot 
                    if abs(vureceiveauvrun_time(iNA,iarrival)-vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)*n*t_data)< t_data  
                        vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                    end
                end
                
                
                
                
                
                % Case 2: Receiving-Receiving Collisions at the  next-hop node
                if vusendauvrun_node_random(iNA,iarrival)<= N0+N1   
                    for jn = 1:length(path(vusendauvrun_node_random(iNA,iarrival)).sonnode)  
                        vusendnoderun_temp = path(vusendauvrun_node_random(iNA,iarrival)).sonnode(jn);
                        if ismember(vusendnoderun_temp, timeslot(vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)).nodeset)==1 
                            vureceivesnrun_timecase2 = vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)*n*t_data+dis_table(vusendnoderun_temp,vusendauvrun_node_random(iNA,iarrival))/vs; 
                            if abs(vureceiveauvrun_time(iNA,iarrival)-vureceivesnrun_timecase2)<t_data
                                vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                            end
                        end
                    end
                end
                
                % Case 3: Receiving-Sending Collisions at  the  interference  nodes
                
                for jn = 1:vuknown_nodesincount(iNA,iarrival)    
                    if vuAUVrun_knownnodesin(iNA,iarrival,jn)~= vusendauvrun_node_random(iNA,iarrival)
                        if ismember(vuAUVrun_knownnodesin(iNA,iarrival,jn), timeslot(vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)).nodeset)==1 
                            vuAUVrun_distancenode(iNA,iarrival,jn) = sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2);
                            vureceiveauvrun_timecase3 = vusendauvrun_time_random(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,jn)/vs; 
                            if abs(vureceiveauvrun_timecase3-floor(vureceiveauvrun_timecase3/n/t_data)*n*t_data)<t_data
                                vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                            end
                        end
                    end
                end
                
                
                % Case 4: Receiving-Receiving Collisions at the interference  nodes
                
                for jn = 1:vuknown_nodesincount(iNA,iarrival)    
                    if vuAUVrun_knownnodesin(iNA,iarrival,jn)~= vusendauvrun_node_random(iNA,iarrival)
                        if vuAUVrun_knownnodesin(iNA,iarrival,jn)<= N0+N1   
                            for kn = 1:length(path(vuAUVrun_knownnodesin(iNA,iarrival,jn)).sonnode)  
                                vusendnoderun_temp = path(vuAUVrun_knownnodesin(iNA,iarrival,jn)).sonnode(kn);
                                if ismember(vusendnoderun_temp, timeslot(vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)).nodeset)==1 
                                    vuAUVrun_distancenode(iNA,iarrival,jn) = sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2);
                                    vureceivesnrun_timecase4 = vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)*n*t_data+dis_table(vusendnoderun_temp,vuAUVrun_knownnodesin(iNA,iarrival,jn))/vs; 
                                    vureceiveauvrun_timecase4 = vusendauvrun_time_random(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,jn)/vs; 
                                    
                                    if abs(vureceiveauvrun_timecase4-vureceivesnrun_timecase4)<t_data
                                        vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                                    end
                                end
                            end
                        end
                    end
                end
                
                vuqueueout(iNA,iarrival) =  vuqueueout(iNA,iarrival-1)+1;
                vuqueueindex(iNA,iarrival) = vuqueuecount(iNA,iarrival)-(vuqueuecount(iNA,iarrival-1)-vuqueueout(iNA,iarrival-1));
                vuqueuesendtime(iNA,vuqueueindex(iNA,iarrival)) = vusendauvrun_time_random(iNA,iarrival);
                
                 
            elseif vuknown_nodesincount(iNA,iarrival)==0 && iarrival< vutotaltimearrival(iNA,ivutrafficload)
                
                vunexttimearrivaltimeslotnumber = floor(( vuarrivaltime(iNA,iarrival)+vutime_interval(iNA,iarrival+1,iteration,ivutrafficload))/(n*t_data));
                vuflag_knownnodesincount = 0;
                if vunexttimearrivaltimeslotnumber>vutimeslot_number(iNA,iarrival)
                    vuflag_knownnodesincount = 0;
                    for vutimeslot_number_temp = vutimeslot_number(iNA,iarrival)+1:1:vunexttimearrivaltimeslotnumber%%%%%%%%%%%%%%%
                        
                        
                        
                         
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
                            
                          
                            for i = 1:length(vuAUVrun_in(iNA,iarrival,:))
                                if vuAUVrun_in(iNA,iarrival,i)~= 0  
                                    if vuAUVrun_detected(iNA,iarrival,vuAUVrun_in(iNA,iarrival,i))==1 
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
                                if vuAUVrun_detected(iNA,iarrival,iN)==1 && ismember(iN, vuAUVrun_in(iNA,iarrival,i))==0
                                    vuknown_nodescount(iNA,iarrival) = vuknown_nodescount(iNA,iarrival)+1;
                                    vuAUVrun_knownnodes(iNA,iarrival,vuknown_nodescount(iNA,iarrival)) = iN;
                                elseif vuAUVrun_detected(iNA,iarrival,iN)==0 && ismember(iN, vuAUVrun_in(iNA,iarrival,i))==0
                                    vuunknown_nodescount(iNA,iarrival) = vuunknown_nodescount(iNA,iarrival)+1;
                                    vuAUVrun_unknownnodes(iNA,iarrival,vuunknown_nodescount(iNA,iarrival)) = vuAUVrun_in(iNA,iarrival,i);
                                end
                            end
                            
                            if vuknown_nodesincount(iNA,iarrival)>= 1
                                vusendauvrun_time_random(iNA,iarrival) = vutimeslot_number(iNA,iarrival)*n*t_data + (5*n*t_data) * rand(1,1);
                                vusendauvrun_node_random(iNA,iarrival) = vuAUVrun_knownnodesin(iNA,iarrival,randi([1,vuknown_nodesincount(iNA,iarrival)],1,1));
                                
                                vusendauvrun_timeslot_number(iNA,iarrival) = floor(vusendauvrun_time_random(iNA,iarrival)/(n*t_data));
                                vusendauvrun_timeslot_numberplus(iNA,iarrival) = mod(vusendauvrun_time_random(iNA,iarrival),(n*t_data));
                                
                                vuAUVrun_distancenode(iNA,iarrival,vusendauvrun_node_random(iNA,iarrival)) =  sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx( vusendauvrun_node_random(iNA,iarrival)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny( vusendauvrun_node_random(iNA,iarrival)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz( vusendauvrun_node_random(iNA,iarrival)))^2);
                                
                                vureceiveauvrun_time(iNA,iarrival) = vusendauvrun_time_random(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,vusendauvrun_node_random(iNA,iarrival))/vs; 
                                
                                vureceiveauvrun_timecase1timeslot_number(iNA,iarrival) = floor(vureceiveauvrun_time(iNA,iarrival)/(n*t_data));
                                vureceiveauvrun_timecase1timeslot_numberplus(iNA,iarrival) = mod(vureceiveauvrun_time(iNA,iarrival),(n*t_data));
                                
                                
                                
                                % Case1: Receiving-Sending Collisions at  the next-hop node
                                
                                
                                if  ismember(vusendauvrun_node_random(iNA,iarrival), timeslot(vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)).nodeset)==1 
                                    if abs(vureceiveauvrun_time(iNA,iarrival)-vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)*n*t_data)< t_data  
                                        vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                                    end
                                end
                                
                                
                                
                                
                                
                                % Case 2: Receiving-Receiving Collisions at the  next-hop node
                                if vusendauvrun_node_random(iNA,iarrival)<= N0+N1   
                                    for jn = 1:length(path(vusendauvrun_node_random(iNA,iarrival)).sonnode) 
                                        vusendnoderun_temp = path(vusendauvrun_node_random(iNA,iarrival)).sonnode(jn);
                                        if vusendnoderun_temp~= 0
                                        if ismember(vusendnoderun_temp, timeslot(vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)).nodeset)==1 
                                            vureceivesnrun_timecase2 = vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)*n*t_data+dis_table(vusendnoderun_temp,vusendauvrun_node_random(iNA,iarrival))/vs; 
                                            if abs(vureceiveauvrun_time(iNA,iarrival)-vureceivesnrun_timecase2)<t_data
                                                vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                                            end
                                        end
                                        end
                                    end
                                end
                                
                                % Case 3: Receiving-Sending Collisions at  the  interference  nodes
                                
                                for jn = 1:vuknown_nodesincount(iNA,iarrival)   
                                    if vuAUVrun_knownnodesin(iNA,iarrival,jn)~= vusendauvrun_node_random(iNA,iarrival)
                                        if ismember(vuAUVrun_knownnodesin(iNA,iarrival,jn), timeslot(vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)).nodeset)==1 
                                            vuAUVrun_distancenode(iNA,iarrival,jn) = sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2);
                                            vureceiveauvrun_timecase3 = vusendauvrun_time_random(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,jn)/vs; 
                                            if abs(vureceiveauvrun_timecase3-floor(vureceiveauvrun_timecase3/n/t_data)*n*t_data)<t_data
                                                vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                                            end
                                        end
                                    end
                                end
                                
                                
                                % Case 4: Receiving-Receiving Collisions at the interference  nodes
                                
                                for jn = 1:vuknown_nodesincount(iNA,iarrival)   
                                    if vuAUVrun_knownnodesin(iNA,iarrival,jn)~= vusendauvrun_node_random(iNA,iarrival)
                                        if vuAUVrun_knownnodesin(iNA,iarrival,jn)<= N0+N1   
                                            for kn = 1:length(path(vuAUVrun_knownnodesin(iNA,iarrival,jn)).sonnode)  
                                                vusendnoderun_temp = path(vuAUVrun_knownnodesin(iNA,iarrival,jn)).sonnode(kn);
                                                if vusendnoderun_temp~= 0
                                                if ismember(vusendnoderun_temp, timeslot(vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)).nodeset)==1 
                                                    vuAUVrun_distancenode(iNA,iarrival,jn) = sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz(vuAUVrun_knownnodesin(iNA,iarrival,jn)))^2);
                                                    vureceivesnrun_timecase4 = vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)*n*t_data+dis_table(vusendnoderun_temp,vuAUVrun_knownnodesin(iNA,iarrival,jn))/vs;
                                                    vureceiveauvrun_timecase4 = vusendauvrun_time_random(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,jn)/vs; 
                                                    
                                                    if abs(vureceiveauvrun_timecase4-vureceivesnrun_timecase4)<t_data
                                                        vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                                                    end
                                                end
                                                end
                                            end
                                        end
                                    end
                                end
                                vuqueueout(iNA,iarrival) =  vuqueueout(iNA,iarrival-1)+1;
                                vuqueueindex(iNA,iarrival) = vuqueuecount(iNA,iarrival)-(vuqueuecount(iNA,iarrival-1)-vuqueueout(iNA,iarrival-1));
                                vuqueuesendtime(iNA,vuqueueindex(iNA,iarrival)) = vusendauvrun_time_random(iNA,iarrival);
                                
                            end
                            vuflag_knownnodesincount = 1;
                            break; 
                            
                        end
                    end
                end
                if vuflag_knownnodesincount==0
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

 
%   Collision Calculation  among multiple AUVs
vuauvcollisioncount = zeros(NA,max(max(vutotaltimearrival)));
for iNA = 1:NA
    for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload)  
        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)
            
            
            for jNA = iNA+1:NA
                for jarrival = 1:vutotaltimearrival(jNA,ivutrafficload)
                                   
                        if vutimeslot_number(iNA,iarrival)==vutimeslot_number(jNA,jarrival) &&vusendauvrun_node_random(iNA,iarrival)~= 0&& vusendauvrun_node_random(jNA,jarrival)~= 0
                            if abs(vusendauvrun_time_random(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,vusendauvrun_node_random(iNA,iarrival))/vs-vusendauvrun_time_random(jNA,jarrival)-vuAUVrun_distancenode(jNA,jarrival,vusendauvrun_node_random(jNA,jarrival))/vs)<t_data  
                             vuauvcollisioncount(iNA,iarrival) = vuauvcollisioncount(iNA,iarrival)+1;
                             vuauvcollisioncount(jNA,jarrival) = vuauvcollisioncount(jNA,jarrival)+1;
                            end
                        end
                end
            end
        end
    end
end



%%  Calculate the network throughput, network delay, VoI, congestion ratio, and collision probability


for iNA = 1:NA
    
        for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload) 
            if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)
                 if vusendauvrun_node_random(iNA,iarrival)~= 0
                    if vuqueuesendtime(iNA,iarrival)>0
                        vuqueuesendtimebefore(iNA,iarrival) = floor(vuqueuesendtime(iNA,iarrival)/timeplus);
                        vuqueuesendtimeplus(iNA,iarrival) = mod(vuqueuesendtime(iNA,iarrival),timeplus);
                        
                        if vuqueuesendtime(iNA,iarrival)>0 && vuqueuesendtime(iNA,iarrival)<len_xyz(iNA)*timeplus-timestart-timeend
                            vuAx_packetsendtime(iNA,iarrival) = Ax(iNA,vuqueuesendtimebefore(iNA,iarrival))+(Ax(iNA,vuqueuesendtimebefore(iNA,iarrival)+1)-Ax(iNA,vuqueuesendtimebefore(iNA,iarrival)))/timeplus*vuqueuesendtimeplus(iNA,iarrival); 
                            vuAy_packetsendtime(iNA,iarrival) = Ay(iNA,vuqueuesendtimebefore(iNA,iarrival))+(Ay(iNA,vuqueuesendtimebefore(iNA,iarrival)+1)-Ay(iNA,vuqueuesendtimebefore(iNA,iarrival)))/timeplus*vuqueuesendtimeplus(iNA,iarrival);
                            vuAz_packetsendtime(iNA,iarrival) = Az(iNA,vuqueuesendtimebefore(iNA,iarrival))+(Az(iNA,vuqueuesendtimebefore(iNA,iarrival)+1)-Az(iNA,vuqueuesendtimebefore(iNA,iarrival)))/timeplus*vuqueuesendtimeplus(iNA,iarrival);
                            vuAUVrun_sendtimedistancenode(iNA,iarrival) =  sqrt((vuAx_packetsendtime(iNA,iarrival)-Nx(vusendauvrun_node_random(iNA,iarrival)))^2+(vuAy_packetsendtime(iNA,iarrival)-Ny(vusendauvrun_node_random(iNA,iarrival)))^2+(vuAz_packetsendtime(iNA,iarrival)-Nz(vusendauvrun_node_random(iNA,iarrival)))^2);
                            
                        end
                    end
                 end
            end
        end
     
end


vutotaldelay = zeros(NA,max(max(vutotaltimearrival)));
vupacketerrorate = zeros(NA,max(max(vutotaltimearrival)));
vusenddelay = zeros(NA,max(max(vutotaltimearrival)));
for iNA = 1:NA
    for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload) %每次数据到达的时候时间
        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)
            if vusendauvrun_node_random(iNA,iarrival)~= 0
                if vuqueuesendtime(iNA,iarrival)>0 && vuqueuesendtime(iNA,iarrival)<len_xyz(iNA)*timeplus-timestart-timeend
                    vusenddelay(iNA,iarrival) = vuqueuesendtime(iNA,iarrival)-vuqueuearrivaltime(iNA,iarrival);
                    vupacketerrorate(iNA,iarrival) = packeterroratefunction(vuAUVrun_sendtimedistancenode(iNA,iarrival),L_packet); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if vusendauvrun_node_random(iNA,iarrival)<= N0  %第一层
                        vutotaldelay(iNA,iarrival) =  vusenddelay(iNA,iarrival);
                        vupacketerrorate(iNA,iarrival) = vupacketerrorate(iNA,iarrival);
                    elseif vusendauvrun_node_random(iNA,iarrival)>N0 && vusendauvrun_node_random(iNA,iarrival)<= N0+N1  %在第二层
                        vutimeslot_temp = floor(vuqueuesendtime(iNA,iarrival)/(n*t_data));
                        vuround_temp = floor(vutimeslot_temp/path(vusendauvrun_node_random(iNA,iarrival)).count);
                        vuroundplus_temp = mod(vutimeslot_temp,path(vusendauvrun_node_random(iNA,iarrival)).count);
                        if node(vusendauvrun_node_random(iNA,iarrival)).seq < vuroundplus_temp
                            vutotaldelay(iNA,iarrival) = vuqueuesendtime(iNA,iarrival)+(path(vusendauvrun_node_random(iNA,iarrival)).count-vuroundplus_temp)*n*t_data+(node(vusendauvrun_node_random(iNA,iarrival)).seq-1)*n*t_data+dis_table(vusendauvrun_node_random(iNA,iarrival),path(vusendauvrun_node_random(iNA,iarrival)).hop)/vs-vuqueuearrivaltime(iNA,iarrival);
                        elseif node(vusendauvrun_node_random(iNA,iarrival)).seq >=  vuroundplus_temp
                            vutotaldelay(iNA,iarrival) = vuqueuesendtime(iNA,iarrival)+(node(vusendauvrun_node_random(iNA,iarrival)).seq-vuroundplus_temp)*n*t_data+dis_table(vusendauvrun_node_random(iNA,iarrival),path(vusendauvrun_node_random(iNA,iarrival)).hop)/vs-vuqueuearrivaltime(iNA,iarrival);
                        end
                        vupacketerrorate(iNA,iarrival) = vupacketerrorate(iNA,iarrival)*packeterroratefunction(dis_table(vusendauvrun_node_random(iNA,iarrival),path(vusendauvrun_node_random(iNA,iarrival)).hop),L_packet); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    elseif vusendauvrun_node_random(iNA,iarrival)>N0+N1   %在第三层
                        vutimeslot_temp = floor(vuqueuesendtime(iNA,iarrival)/(n*t_data));
                        vuround_temp = floor(vutimeslot_temp/path(vusendauvrun_node_random(iNA,iarrival)).count);
                        vuroundplus_temp = mod(vutimeslot_temp,path(vusendauvrun_node_random(iNA,iarrival)).count);
                        vupacketerrorate(iNA,iarrival) = vupacketerrorate(iNA,iarrival)*packeterroratefunction(dis_table(vusendauvrun_node_random(iNA,iarrival),path(vusendauvrun_node_random(iNA,iarrival)).hop),L_packet)*packeterroratefunction(dis_table(path(vusendauvrun_node_random(iNA,iarrival)).hop,path(path(vusendauvrun_node_random(iNA,iarrival)).hop).hop),L_packet); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if node(vusendauvrun_node_random(iNA,iarrival)).seq <=  vuroundplus_temp
                            vuhop1time = vuqueuesendtime(iNA,iarrival)+(path(vusendauvrun_node_random(iNA,iarrival)).count-vuroundplus_temp)*n*t_data+(node(vusendauvrun_node_random(iNA,iarrival)).seq-1)*n*t_data+dis_table(vusendauvrun_node_random(iNA,iarrival),path(vusendauvrun_node_random(iNA,iarrival)).hop)/vs;
                            vuhop2timeslot = floor(vuhop1time/(n*t_data));
                            vuhop2round = floor(vuhop2timeslot/path(path(vusendauvrun_node_random(iNA,iarrival)).hop).count);
                            vuhop2roundplus = mod(vuhop2timeslot,path(path(vusendauvrun_node_random(iNA,iarrival)).hop).count);
                            if node(path(vusendauvrun_node_random(iNA,iarrival)).hop).seq <=  vuhop2roundplus
                                vuhop2time = (path(path(vusendauvrun_node_random(iNA,iarrival)).hop).count-vuhop2roundplus)*n*t_data+(node(path(vusendauvrun_node_random(iNA,iarrival)).hop).seq-1)*n*t_data+dis_table(path(vusendauvrun_node_random(iNA,iarrival)).hop,path(path(vusendauvrun_node_random(iNA,iarrival)).hop).hop)/vs;
                                vutotaldelay(iNA,iarrival) = vuhop1time+vuhop2time-vuqueuearrivaltime(iNA,iarrival);
                            elseif node(path(vusendauvrun_node_random(iNA,iarrival)).hop).seq > vuhop2roundplus
                                vuhop2time = (node(path(vusendauvrun_node_random(iNA,iarrival)).hop).seq-vuhop2roundplus)*n*t_data+dis_table(path(vusendauvrun_node_random(iNA,iarrival)).hop,path(path(vusendauvrun_node_random(iNA,iarrival)).hop).hop)/vs;
                                vutotaldelay(iNA,iarrival) = vuhop1time+vuhop2time-vuqueuearrivaltime(iNA,iarrival);
                            end
                        elseif node(vusendauvrun_node_random(iNA,iarrival)).seq > vuroundplus_temp
                            vuhop1time = vuqueuesendtime(iNA,iarrival)+(node(vusendauvrun_node_random(iNA,iarrival)).seq-vuroundplus_temp)*n*t_data+dis_table(vusendauvrun_node_random(iNA,iarrival),path(vusendauvrun_node_random(iNA,iarrival)).hop)/vs;
                            vuhop2timeslot = floor(vuhop1time/(n*t_data));
                            vuhop2round = floor(vuhop2timeslot/path(path(vusendauvrun_node_random(iNA,iarrival)).hop).count);
                            vuhop2roundplus = mod(vuhop2timeslot,path(path(vusendauvrun_node_random(iNA,iarrival)).hop).count);
                            if node(path(vusendauvrun_node_random(iNA,iarrival)).hop).seq <=  vuhop2roundplus
                                vuhop2time = (path(path(vusendauvrun_node_random(iNA,iarrival)).hop).count-vuhop2roundplus)*n*t_data+(node(path(vusendauvrun_node_random(iNA,iarrival)).hop).seq-1)*n*t_data+dis_table(path(vusendauvrun_node_random(iNA,iarrival)).hop,path(path(vusendauvrun_node_random(iNA,iarrival)).hop).hop)/vs;
                                vutotaldelay(iNA,iarrival) = vuhop1time+vuhop2time-vuqueuearrivaltime(iNA,iarrival);
                            elseif node(path(vusendauvrun_node_random(iNA,iarrival)).hop).seq > vuhop2roundplus
                                vuhop2time = (node(path(vusendauvrun_node_random(iNA,iarrival)).hop).seq-vuhop2roundplus)*n*t_data+dis_table(path(vusendauvrun_node_random(iNA,iarrival)).hop,path(path(vusendauvrun_node_random(iNA,iarrival)).hop).hop)/vs;
                                vutotaldelay(iNA,iarrival) = vuhop1time+vuhop2time-vuqueuearrivaltime(iNA,iarrival);
                            end
                        end
                    end
                end
            end
        end
    end
end


%% Calculate VoI
vuvoi = zeros(NA,max(max(vutotaltimearrival)));
vut = zeros(NA,max(max(vutotaltimearrival)));
for iNA = 1:NA
    for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload) %每次数据到达的时候时间
        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)
                if vuqueuesendtime(iNA,iarrival)>0&& vuqueuesendtime(iNA,iarrival)<len_xyz(iNA)*timeplus-timestart-timeend
                    vut(iNA,iarrival)  =  valueofinformation(vusenddelay(iNA,iarrival),betau,vu0,alphau,Tu);
                    vuvoi(iNA,iarrival)  =  vut(iNA,iarrival)*(1-vupacketerrorate(iNA,iarrival));
                end
            
        end
    end
end










%% The collision for the unknown static nodes 
vuAUVrun_undistancenode = zeros(NA,max(max(vutotaltimearrival)),N);
for iNA = 1:NA
    for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload)

        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)-1 &&vusendauvrun_node_random(iNA,iarrival)~= 0
            
            for jn = 1:vuunknown_nodesincount(iNA,iarrival)
                vuAUVrun_undistancenode(iNA,iarrival,jn) =  sqrt((vuAx_packetarrivaltime(iNA,iarrival)-Nx(vuAUVrun_unknownnodesin(iNA,iarrival,jn)))^2+(vuAy_packetarrivaltime(iNA,iarrival)-Ny(vuAUVrun_unknownnodesin(iNA,iarrival,jn)))^2+(vuAz_packetarrivaltime(iNA,iarrival)-Nz(vuAUVrun_unknownnodesin(iNA,iarrival,jn)))^2);
                
                % Case1: Receiving-Sending Collisions at  the next-hop node
           
                vuunreceiveauvrun_timecase1 = vusendauvrun_time_random(iNA,iarrival)+vuAUVrun_undistancenode(iNA,iarrival,jn)/vs;
                if abs(vureceiveauvrun_time(iNA,iarrival)-vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)*n*t_data)< t_data  && ismember(vusendauvrun_node_random(iNA,iarrival), timeslot(vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)+1).nodeset)==1
                    vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                    break;
              
                end
                
                % Case 2: Receiving-Receiving Collisions at the  next-hop node
            
                
                
                % Case 3: Receiving-Sending Collisions at  the interference  nodes
            
                % Case 4: Receiving-Receiving Collisions at the interference  nodes
                
                if vusendauvrun_node_random(iNA,iarrival)<= N0+N1  
                    for sn = 1:length(path(vusendauvrun_node_random(iNA,iarrival)).sonnode) 
                        vusendnoderun_temp = path(vusendauvrun_node_random(iNA,iarrival)).sonnode(sn);
                        if ismember(vusendnoderun_temp, timeslot(vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)).nodeset)==1
                            vuunreceivesnrun_timecase2 = vureceiveauvrun_timecase1timeslot_number(iNA,iarrival)*n*t_data+dis_table(vusendnoderun_temp,vusendauvrun_node_random(iNA,iarrival))/vs;
                            vuunreceiveauvrun_timecase2 = vusendauvrun_time_random(iNA,iarrival)+vuAUVrun_distancenode(iNA,iarrival,jn)/vs;
                            if abs(vuunreceivesnrun_timecase2-vuunreceiveauvrun_timecase2)<t_data
                                vucollisioncount(iNA,iarrival) = vucollisioncount(iNA,iarrival)+1;
                                break;
                            end
                            
                        end
                    end
                end
            end
        end
    end
end

vuvoi_collision = vuvoi;
for iNA = 1:NA
    for iarrival = 1:vutotaltimearrival(iNA,ivutrafficload) 
        if vuAUVtimebefore(iNA,iarrival)>0 && vuAUVtimebefore(iNA,iarrival)<len_xyz(iNA)-1 &&vusendauvrun_node_random(iNA,iarrival)~= 0
            
            if vucollisioncount(iNA,iarrival)~= 0 || vuauvcollisioncount(iNA,iarrival)~= 0
                vuvoi_collision(iNA,iarrival)  =  0;
            end
        end
    end
end 


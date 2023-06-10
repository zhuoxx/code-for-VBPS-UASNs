function [delay_table,dis_table,TM,IM,timeslot,sequence,node,path,Ts] = generate_TDMA(N,Nx,Ny,Nz,dc,vs,N0,N1,N2,path1,n,t_data,Round,totaltimeslot,iL_packet)

 
% % %%Gnerate the network topology for static nodes
% % %The item in the topology matrix is 1 when two nodes has the communication
% % %connections, and 0, otherwise;
% % %The item in the interference matrix is 1 when the distance between two nodes is
% % %smaller than the communication range;
% % %The item in the delay table is the propagation delay between nodes;
% % %The item in the distance table is the propagation delay between nodes;
% % 
%%---------------------Initialization-------------------------------------- 
TM = zeros(N,N);            %The topology matrix for static nodes
IM = zeros(N,N);            %The interference matrix for static nodes
delay_table = zeros(N,N);   %The delay table for static nodes 
dis_table = zeros(N,N);     %The distance table for static nodes 

for i=1:N
    node(i).interference=[];
end

%%-------Calculate the distance and interference matrix between nodes------
for i = 1:N
    for j = i+1:N
        distance = sqrt(power((Nx(i)-Nx(j)),2)+power((Ny(i)-Ny(j)),2)+power((Nz(i)-Nz(j)),2));
        dis_table(i,j) = distance;
        dis_table(j,i) = distance;
        if distance<=dc
            IM(i,j) = 1;
            IM(j,i) = 1;
            
            node(i).interference = union(node(i).interference,j);
            node(j).interference = union(node(j).interference,i);
            
            delay_table(i,j) = distance/vs;
            delay_table(j,i) = delay_table(i,j);
        else
            delay_table(i,j) = inf;
            delay_table(j,i) = delay_table(i,j);
        end
    end
end


%%---------------------Calculate the topology matrix--------------- 
for i = 1:N0
    TM(i,i) = 1;
end
for i = 1:N
    path(i).sonnode = [];
end
for i = N0+1:N0+N1
    min_delay = delay_table(i,1);
    path(i).source = i;
    path(i).hop = 1;
    index = 1;
    for j = 2:N0
        if delay_table(i,j)<min_delay
            path(i).hop = j;
            index = j;
        end
    end
    TM(i,index) = 1;
    TM(index,i) = 1;
    path(index).sonnode = union(path(index).sonnode,i);
end

for i = N0+N1+1:N0+N1+N2
    min_delay = delay_table(i,N0+1);
    path(i).source = i;
    path(i).hop = N0+1;
    index = N0+1;
    for j = 2+N0:N0+N1
        if IM(i,j) == 1
            
            if delay_table(i,j)<min_delay
                path(i).hop = j;
                index = j;
            end
        end
    end
    TM(i,index) = 1;
    TM(index,i) = 1;
    path(index).sonnode = union(path(index).sonnode,i);
end

%%Generate the sending time for each static nodes based on TDMA protocol
 
%The interference nodes should not send data packets in the same time slots.
% So the determination of sending sequence of static nodes should consider
% the collisions between  interference nodes. 

%%--------------------Calculate the collisons between sensor nodes---------
 
 
for i = N0+N1+1:N0+N1+N2
    for j = N0+1:N0+N1
        if IM(i,j)==1 && TM(i,j)~=1
            temp = find(TM(:,j)==1);
            for k = 1:length(temp)
                if temp(k)>N0  
                    if abs(delay_table(temp(k),j)-delay_table(i,j))<t_data
                        node(i).conflict = temp(k);
                        node(temp(k)).conflict = union(node(temp(k)).conflict,i);
                    end
                end
            end
        end
    end
end

%%----------------Calculate the collisons between relay nodes---------  
for i = N0+1:N0+N1
    for j = 1:N0
        if IM(i,j)==1 && TM(i,j)~=1
            temp = find(TM(:,j)==1);
            for k = 1:length(temp)
                if abs(delay_table(temp(k),j)-delay_table(i,j))<t_data
                    node(i).conflict = temp(k);
                    node(temp(k)).conflict = union(node(temp(k)).conflict,i);
                end
            end
        end
    end
end

%%-----Calculate the sending sequence of sensor nodes and relay nodes------
sequence = zeros(N,N);
for i = 1:N
    node(i).seq = 0;
end
 
for i = N0+1:N0+N1+N2
    if node(i).seq==0
        if ~isempty(node(i).conflict)
            for j = 1:length(node(i).conflict) %If the collision set is larger than 1, the sending sequence should be calculated for all collision nodes iteratively.  
                if node(node(i).conflict(j)).seq==0
                    path(node(i).conflict(j)).count = sum(TM(:,path(node(i).conflict(j)).hop))-1;
                    path(i).count = sum(TM(:,path(i).hop))-1;
                    %The least common multiple of two nodes min_GBS
                    if path(i).count> path(node(i).conflict(j)).count
                        count_max = path(i).count;
                    else
                        count_max = path(node(i).conflict(j)).count;
                    end
                    for b = count_max:path(i).count*path(node(i).conflict(j)).count
                        if mod(b,path(i).count)==0 && mod(b,path(node(i).conflict(j)).count)==0
                            min_GBS = b;
                            break;
                        end
                    end
                    if j==1  %Determine the first collision nodes
                        %Calculate the sequence of collision-free  
                        for s = 1:path(i).count
                            for a = 1:path(node(i).conflict(j)).count
                                for k = 0:min_GBS/path(i).count
                                    if mod(path(i).count*k+s-a,path(node(i).conflict(j)).count)==0 %Judge if two collision nodes send with the same time slot判断是否会在同一个slot
                                        %   node(i).seq=0;
                                        %   node(node(i).conflict(j)).seq=0;
                                        %   sequence(path(i).hop,s)=0;
                                        %   sequence(path(node(i).conflict(j)).hop,a)=0;
                                        break;
                                    end
                                end
                                %If two collision nodes is not in the same slot and if the sequence have not been assigend, then assign the sequence for this two nodes. 
                                if  k==min_GBS/path(i).count && sequence(path(i).hop,s)==0 && sequence(path(node(i).conflict(j)).hop,a)==0
                                    node(i).seq = s;
                                    node(node(i).conflict(j)).seq = a;
                                    sequence(path(i).hop,s) = i;
                                    sequence(path(node(i).conflict(j)).hop,a) = node(i).conflict(j);
                                    break;
                                end
                            end
                            
                            if node(i).seq~=0
                                break;
                            end
                        end
                    else  %Determine the other collision nodes
                        for a = 1:path(node(i).conflict(j)).count
                            for k = 0:min_GBS/path(i).count
                                if mod(path(i).count*k+node(i).seq-a,path(node(i).conflict(j)).count)==0  
                                    %                                 node(i).seq=0;
                                    %                                 node(node(i).conflict(j)).seq=0;
                                    %                                 sequence(path(i).hop,s)=0;
                                    %                                 sequence(path(node(i).conflict(j)).hop,a)=0;
                                    break;
                                end
                            end
                             
                            if  k==min_GBS/path(i).count && sequence(path(node(i).conflict(j)).hop,a)==0
                                node(node(i).conflict(j)).seq = a;
                                sequence(path(node(i).conflict(j)).hop,a) = node(i).conflict(j);
                                break;
                            end
                        end
                    end
                end
            end
            
        end
    end
end


for i = N0+1:N0+N1+N2
    if node(i).seq==0
        path(i).count = sum(TM(:,path(i).hop))-1;
        for s = 1:path(i).count
            if sequence(path(i).hop,s)==0
                node(i).seq = s;
                sequence(path(i).hop,s) = i;
                break;
            end
        end
    end
end


%%--------------Calculate TDMA sending time for static nodes--------------- 

Ts = zeros(N,Round);
for iround = 1:Round
    for i = N0+1:N0+N1+N2
        Ts(i,iround) = ((iround-1)*path(i).count+(node(i).seq-1))*n*t_data;  %Sending time
    end
end

for itimeslot = 1:totaltimeslot
    timeslot(itimeslot).nodeset = [];
    for iN = N0+1:N0+N1+N2
        if mod(itimeslot-node(iN).seq,path(iN).count)==0
            timeslot(itimeslot).nodeset = union(timeslot(itimeslot).nodeset,iN);
        end
    end
end

%%------Save and load gnerated network topology for static nodes-----------

% If  the locations of static nodes is changed, you can save the new
% generated topology matrix.
% save(fullfile(path1,'TM3'),'TM');
% save(fullfile(path1,'IM3'),'IM');
% save(fullfile(path1,'delay_table3'),'delay_table');
% save(fullfile(path1,'dis_table3'),'dis_table');
 
 

% save(fullfile(path1,strcat('timeslot',num2str(iL_packet))),'timeslot');
% save(fullfile(path1,strcat('sequence',num2str(iL_packet))),'sequence');
% save(fullfile(path1,strcat('node',num2str(iL_packet))),'node');
% save(fullfile(path1,strcat('path',num2str(iL_packet))),'path');
% save(fullfile(path1,strcat('Ts',num2str(iL_packet))),'Ts');



delay_table = cell2mat(struct2cell(load(fullfile(path1,'delay_table3'))));
dis_table = cell2mat(struct2cell(load(fullfile(path1,'dis_table3'))));
TM = cell2mat(struct2cell(load(fullfile(path1,'TM3'))));
IM = cell2mat(struct2cell(load(fullfile(path1,'IM3'))));
 


timeslot = cell2mat(struct2cell(load(fullfile(path1,strcat('timeslot',num2str(iL_packet))))));
sequence = cell2mat(struct2cell(load(fullfile(path1,strcat('sequence',num2str(iL_packet))))));
node = cell2mat(struct2cell(load(fullfile(path1,strcat('node',num2str(iL_packet))))));
path = cell2mat(struct2cell(load(fullfile(path1,strcat('path',num2str(iL_packet))))));
Ts = cell2mat(struct2cell(load(fullfile(path1,strcat('Ts',num2str(iL_packet))))));
%%-------------------Plot the locations for static nodes-------------------
% figure
% cmap=colormap;
% scatter3(Nx(1:N0),Ny(1:N0),Nz(1:N0),'r');
% hold on
% scatter3(Nx(N0+1:N0+N1),Ny(N0+1:N0+N1),Nz(N0+1:N0+N1),'g');
% hold on
% scatter3(Nx(N0+N1+1:N0+N1+N2),Ny(N0+N1+1:N0+N1+N2),Nz(N0+N1+1:N0+N1+N2),'b');
% hold on
% % gplot3(TM,Local_N,'-');
% 
% for i = 1:N
%     for j = i+1:N
%         if TM(i,j) == 1
%             line([Nx(i),Nx(j,1)],[Ny(i),Ny(j)],[Nz(i),Nz(j)]);
%             %     text(x(i),y(i),num2str(i),'fontsize',8,'HorizontalAlignment','center');
%         end
%     end
% end
end
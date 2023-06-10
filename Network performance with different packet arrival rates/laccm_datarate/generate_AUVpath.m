function [Ax,Ay,Az,len_xyz,timeplus] = generate_AUVpath(Nx,Ny,Nz, NA,path1)
% %Generate the path of AUV
%  
% nx_max =  max(Nx)-4000;
% nx_min =  min(Nx)+2000;
% ny_max =  max(Ny)-2000;
% ny_min =  min(Ny)+2000;
% nz_max =  max(Nz)-1000;
% nz_min =  min(Nz);
% 
%  
%  
%  
 timeplus =  5;
% V_A =  timeplus*v_A;
% len_xyz =  zeros(1,NA);
% 
% 
% for iNA =  1:NA
%     Ax(iNA,1) =  Nx(mod(iNA,2)+1);
%     Ay(iNA,1) =  Ny(mod(iNA,2)+1);
%     Az(iNA,1) =  nz_max/(NA*5)*iNA+500;
%     jtime =  2;
%     flag =  0;
%     while  Ax(iNA,jtime-1)<=nx_max && Ax(iNA,jtime-1)>=nx_min &&  Ay(iNA,jtime-1)<=ny_max && Ay(iNA,jtime-1)>=ny_min &&  Az(iNA,jtime-1)<=nz_max && Az(iNA,jtime-1)>=nz_min
% 
%         A_v =  [0,0,0];%The moving direction  of the AUV in x, y, and z coordination
%         temp =  randi([1,3],1,1);
% 
% 
%         A_v(temp) =  V_A;%The randome direction
% 
%         Ax(iNA,jtime) =  Ax(iNA,jtime-1)+(-1)^(iNA)*A_v(1);
%         Ay(iNA,jtime) =  Ay(iNA,jtime-1)+A_v(2);
%         Az(iNA,jtime) =  Az(iNA,jtime-1)+A_v(3);
% 
%         jtime =  jtime+1;
% 
%     end
% 
% 
%     if  (Ax(iNA,jtime-1)>=nx_max || Ax(iNA,jtime-1)<=nx_min) &&  (Ay(iNA,jtime-1)>=ny_max || Ay(iNA,jtime-1)<=ny_min)
%         flag =  4; 
%     elseif  (Ax(iNA,jtime-1)>=nx_max || Ax(iNA,jtime-1)<=nx_min) && (Az(iNA,jtime-1)>=nz_max || Az(iNA,jtime-1)<=nz_min)
%         flag =  5;
%     elseif   (Ay(iNA,jtime-1)>=ny_max || Ay(iNA,jtime-1)<=ny_min) &&  (Az(iNA,jtime-1)>=nz_max || Az(iNA,jtime-1)<=nz_min)
%         flag =  6;
%     elseif   (Ax(iNA,jtime-1)>=nx_max || Ax(iNA,jtime-1)<=nx_min)
%         flag =  1;
%     elseif  (Ay(iNA,jtime-1)>=ny_max || Ay(iNA,jtime-1)<=ny_min)
%         flag =  2;
%     elseif  (Az(iNA,jtime-1)>=nz_max || Az(iNA,jtime-1)<=nz_min)
%         flag =  3;
%     end
% 
%     if flag==1
%         while   Ay(iNA,jtime-1)<=ny_max && Ay(iNA,jtime-1)>=ny_min &&  Az(iNA,jtime-1)<=nz_max && Az(iNA,jtime-1)>=nz_min
%             A_v =  [0,0,0]; 
%             temp =  randi([2,3],1,1);
%             A_v(temp) =  V_A; 
% 
%             Ax(iNA,jtime) =  Ax(iNA,jtime-1)+(-1)^(iNA)*A_v(1);
%             Ay(iNA,jtime) =  Ay(iNA,jtime-1)+A_v(2);
%             Az(iNA,jtime) =  Az(iNA,jtime-1)+A_v(3);
% 
%             jtime =  jtime+1;
%         end
%         if  (Ay(iNA,jtime-1)>=ny_max || Ay(iNA,jtime-1)<=ny_min)
%             flag =  4;
%         elseif (Az(iNA,jtime-1)>=nz_max || Az(iNA,jtime-1)<=nz_min)
%             flag =  5;
%         end
%     end
% 
% 
%     if flag==2
%         while   Ax(iNA,jtime-1)<=nx_max && Ax(iNA,jtime-1)>=nx_min &&  Az(iNA,jtime-1)<=nz_max && Az(iNA,jtime-1)>=nz_min
%             A_v =  [0,0,0]; 
%             a =  [1,3];
%             b =  randi([1 2],1,1);
%             temp =  a(b);
%             A_v(temp) =  V_A; 
% 
%             Ax(iNA,jtime) =  Ax(iNA,jtime-1)+(-1)^(iNA)*A_v(1);
%             Ay(iNA,jtime) =  Ay(iNA,jtime-1)+A_v(2);
%             Az(iNA,jtime) =  Az(iNA,jtime-1)+A_v(3);
% 
%             jtime =  jtime+1;
%         end
%         if  (Ax(iNA,jtime-1)>=nx_max || Ax(iNA,jtime-1)<=nx_min)
%             flag =  4;
%         elseif  (Az(iNA,jtime-1)>=nz_max || Az(iNA,jtime-1)<=nz_min)
%             flag =  6;
%         end
%     end
% 
% 
%     if flag==3
%         while    Ax(iNA,jtime-1)<=nx_max && Ax(iNA,jtime-1)>=nx_min &&  Ay(iNA,jtime-1)<=ny_max && Ay(iNA,jtime-1)>=ny_min
%             A_v =  [0,0,0]; 
%             temp =  randi([1,2],1,1);
%             A_v(temp) =  V_A ;
% 
%             Ax(iNA,jtime) =  Ax(iNA,jtime-1)+(-1)^(iNA)*A_v(1);
%             Ay(iNA,jtime) =  Ay(iNA,jtime-1)+A_v(2);
%             Az(iNA,jtime) =  Az(iNA,jtime-1)+A_v(3);
% 
%             jtime =  jtime+1;
%         end
%         if  (Ax(iNA,jtime-1)>=nx_max || Ax(iNA,jtime-1)<=nx_min)
%             flag =  5;
%         elseif (Ay(iNA,jtime-1)>=ny_max || Ay(iNA,jtime-1)<=ny_min)
%             flag =  6;
%         end
%     end
% 
%     if flag==4
%         while   Az(iNA,jtime-1)<=nz_max && Az(iNA,jtime-1)>=nz_min
%             Ax(iNA,jtime) =  Ax(iNA,jtime-1);
%             Ay(iNA,jtime) =  Ay(iNA,jtime-1);
%             Az(iNA,jtime) =  Az(iNA,jtime-1)+(-1)^(iNA)*V_A;
% 
%             jtime =  jtime+1;
%         end
%     end
% 
%     if flag ==5
%         while  Ay(iNA,jtime-1)<=ny_max && Ay(iNA,jtime-1)>=ny_min
%             Ax(iNA,jtime) =  Ax(iNA,jtime-1);
%             Ay(iNA,jtime) =  Ay(iNA,jtime-1)+(-1)^(iNA)*V_A;
%             Az(iNA,jtime) =  Az(iNA,jtime-1);
% 
%             jtime =  jtime+1;
%         end
%     end
% 
%     if flag==6
%         while   Ax(iNA,jtime-1)<=nx_max && Ax(iNA,jtime-1)>=nx_min &&   Ay(iNA,jtime-1)>=ny_min &&   Az(iNA,jtime-1)>=nz_min %Ax(iNA,jtime-1)<=nx_max && Ax(iNA,jtime-1)>=nx_min
%             
%                 
%                 A_v =  [0,0,0];%指示A_vx,A_vy, A_vz
%                 temp =  randi([1,3],1,1);
%                 
%                 
%                 A_v(temp) =  V_A;%往三个方向中的任意一个方向走
%                 
%                 Ax(iNA,jtime) =  Ax(iNA,jtime-1)+(-1)^(iNA)*A_v(1);
%                 Ay(iNA,jtime) =  Ay(iNA,jtime-1)-A_v(2);
%                 Az(iNA,jtime) =  Az(iNA,jtime-1)-A_v(3);
%                 
%                 jtime =  jtime+1;
%                 
%             
%             
%         end
%     end
%     len_xyz(iNA) =  jtime-1;
% end
%%-------------------Save and load the locations--------------------------
%Since the path of AUVs are predefined, we directly load a predefined 
% locations in the simulation.

% If you want to change the path of AUV, you can save the new
% generated path.
% save((fullfile(path1,'Ax'),'Ax');
% save((fullfile(path1,'Ay'),'Ay');
% save((fullfile(path1,'Az'),'Az');
% save((fullfile(path1,'len_xyz'),'len_xyz');

% path1='D:\software\Matlab\install\work\LAB\zxxpaper\2021auv';
Ax = cell2mat(struct2cell(load(fullfile(path1,'Ax'))));
Ay = cell2mat(struct2cell(load(fullfile(path1,'Ay'))));
Az = cell2mat(struct2cell(load(fullfile(path1,'Az'))));
len_xyz = cell2mat(struct2cell(load(fullfile(path1,'len_xyz'))));

% % Plot the AUV path
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
% for i=1:N
%     for j=i+1:N
%         if TM(i,j)==1
%             line([Nx(i),Nx(j,1)],[Ny(i),Ny(j)],[Nz(i),Nz(j)]);
%             %     text(x(i),y(i),num2str(i),'fontsize',8,'HorizontalAlignment','center');
%         end
%     end
% end
% clr(1,:)=[0 113 188]/255;
% clr(2,:)=[216 82 24]/255;
% clr(3,:)=[236 176 31]/255;
% clr(4,:)=[125 46 141]/255;
% clr(5,:)=[118 171 47]/255;
% clr(6,:)=[76 189 237]/255;
% clr(7,:)=[161 19 46]/255;
% clr(8,:)=[255 255 255]/255;
% hold on
% for iNA=1:NA
%     plot3(Ax(iNA,1:len_xyz(iNA)),Ay(iNA,1:len_xyz(iNA)),Az(iNA,1:len_xyz(iNA)),'color',clr(iNA,:));
% end


%%-------Another alternative algorithm for AUV path planning based on the  ant colony algorithm
% %some nodes are randomly generated and find the shortest path
% NR=10;
% for iNA=1:NA
%     for i=1:NR
%         Nr_x(iNA,i) = randi([nx_min,nx_max],1);
%         Nr_y(iNA,i) = randi([ny_min,ny_max],1);
%         Nr_z(iNA,i) = randi([nz_min,nz_max],1);
%     end
% end
%
% for iNA=1:NA
%     Nr(iNA).distance=zeros(NR,NR);
% end
%
% for iNA=1:NA
%     for i=1:NR
%         for j=i:NR
%             Nr(iNA).distance(i,j)=sqrt(power((Nr_x(i)-Nr_x(j)),2)+power((Nr_y(i)-Nr_y(j)),2)+power((Nr_z(i)-Nr_z(j)),2));
%             Nr(iNA).distance(j,i)= Nr(iNA).distance(i,j);
%         end
%     end
% end
%
% for iNA=1:NA
%     [Nr(iNA).Shortest_Length,Nr(iNA).Shortest_Route] = ACOTSP(Nr(iNA).distance,NR);
% end
%
% % 
% figure
% cmap=colormap;
% scatter3(Nx(1:N0),Ny(1:N0),Nz(1:N0),'r');
% hold on
% scatter3(Nx(N0+1:N0+N1),Ny(N0+1:N0+N1),Nz(N0+1:N0+N1),'g');
% hold on
% scatter3(Nx(N0+N1+1:N0+N1+N2),Ny(N0+N1+1:N0+N1+N2),Nz(N0+N1+1:N0+N1+N2),'b');
% hold on
%
% % gplot3(TM,Local_N,'-');
%
% for i=1:N
%     for j=i+1:N
%         if TM(i,j)==1
%             line([Nx(i),Nx(j,1)],[Ny(i),Ny(j)],[Nz(i),Nz(j)]);
%             %     text(x(i),y(i),num2str(i),'fontsize',8,'HorizontalAlignment','center');
%         end
%     end
% end
% clr(1,:)=[0 113 188]/255;
% clr(2,:)=[216 82 24]/255;
% clr(3,:)=[236 176 31]/255;
% clr(4,:)=[125 46 141]/255;
% clr(5,:)=[118 171 47]/255;
% clr(6,:)=[76 189 237]/255;
% clr(7,:)=[161 19 46]/255;
% clr(8,:)=[255 255 255]/255;
%
% hold on
% for iNA=1:NA
%     clr_nr=repmat(clr(iNA+3,:),NR,1);
%     scatter3(Nr_x(iNA,:),Nr_y(iNA,:),Nr_z(iNA,:),8,clr_nr);%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     hold on
%     for i=2:NR
%     plot3([Nr_x(iNA,Nr(iNA).Shortest_Route(i-1)),Nr_x(iNA,Nr(iNA).Shortest_Route(i))],[Nr_y(iNA,Nr(iNA).Shortest_Route(i-1)),Nr_y(iNA,Nr(iNA).Shortest_Route(i))],[Nr_z(iNA,Nr(iNA).Shortest_Route(i-1)),Nr_z(iNA,Nr(iNA).Shortest_Route(i))],'color',clr(iNA+3,:))
%     hold on
%     end
%
% end
end
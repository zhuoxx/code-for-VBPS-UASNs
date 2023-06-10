function [Nx,Ny,Nz] = generate_staticnode_location(N0,N1,N2,dc,path1)
% Generate the locations for static nodes
% 
% %N0
% N0_x = zeros(N0,1);
% N0_y = zeros(N0,1); 
% N0_z = zeros(N0,1); 
% 
% 
% N0_x(1) = 1500;
% N0_y(1) = 1500;
% N0_z(1) = 0;
% i = 2;
% while N0>=i    
%     N0_x(i) = randi([1500,7000],1);
%     N0_y(i) = 1500;
%     N0_z(i) = 0;
%     if (N0_x(i)-N0_x(i-1))^2+(N0_y(i)-N0_y(i-1))^2+(N0_z(i)-N0_z(i-1))^2> dc^2 
%         i=i+1;
%     end
% end
% 
% %N1
% N1_x = zeros(N1,1);
% N1_y = zeros(N1,1);
% N1_z = zeros(N1,1);
% 
% i = 1;
% while N1>=i    
% 
% 
%     for j=1:N0
%         k=1;
%         while k<=N1/N0
%             N1_x(i) = 2*randi([N0_x(j)-dc,N0_x(j)+dc],1);
%             N1_y(i) = 2*randi([N0_y(j)-dc,N0_y(j)+dc],1);
%             N1_z(i) = 2*randi([500,1500],1);
% 
%             if (N1_x(i)-N0_x(j))^2+(N1_y(i)-N0_y(j))^2+(N1_z(i)-N0_z(j))^2< dc^2 
%                 i = i+1;
%                 k = k+1;
%             end
%             if i>N1
%                 break;
%             end
%         end
%     end
% end
% 
% %N2
% N2_x = zeros(N2,1);
% N2_y = zeros(N2,1);
% N2_z = zeros(N2,1);
% 
% i = 1;
% while N2>=i     
%     for j=1:N1
% 
%         k=1;
%         while k<=N2/N1
%             N2_x(i) = 2*randi([N1_x(j)-dc,N1_x(j)+dc],1);
%             N2_y(i) = 2*randi([N1_y(j)-dc,N1_y(j)+dc],1);
%             N2_z(i) = 2*randi([2000,3000],1);
%             if (N2_x(i)-N1_x(j))^2+(N2_y(i)-N1_y(j))^2+(N2_z(i)-N1_z(j))^2< dc^2 
%                 i=i+1;
%                 k=k+1;
%             end
%             if i>N2
%                 break;
%             end
%         end
% 
%     end
% end
% 
% % figure
% % scatter3(N0_x,N0_y,N0_z,'r');
% % hold on
% % scatter3(N1_x,N1_y,N1_z,'g');
% % hold on
% % scatter3(N2_x,N2_y,N2_z,'b');
% 
% 
% 
% Nx = [N0_x;N1_x;N2_x];
% Ny = [N0_y;N1_y;N2_y];
% Nz = [N0_z;N1_z;N2_z];
% 
% Local_N = [Nx,Ny,Nz];

%%-------------------Save and load the locations--------------------------
%Since the locations of static nodes are predefined, we directly load a predefined 
% static nodes in the simulation.

% If you want to change the locations of static nodes, you can save the new
% generated locations of static nodes.
% save(fullfile(path1,'Nx3'),'Nx');
% save(fullfile(path1,'Ny3'),'Ny');
% save(fullfile(path1,'Nz3'),'Nz');


%path1='D:\software\Matlab\install\work\LAB\zxxpaper\2021auv\results_compare\datarate\vbps_datarate';
Nx = cell2mat(struct2cell(load(fullfile(path1,'Nx3'))));
Ny = cell2mat(struct2cell(load(fullfile(path1,'Ny3'))));
Nz = cell2mat(struct2cell(load(fullfile(path1,'Nz3'))));

end
function [Shortest_Length,Shortest_Route] = ACOTSP(D,n)
%% ��Ⱥ�㷨
% tic;           %��ʱ��ʼ
%%1.��ʼ��
% global D
% global ND
% n=ND; %n��ʾ����Ĺ�ģ�����и�����

% [i,j,k]=size(D_ind);
% D=zeros(n,n);
%     for jj=1:n
%         for kk=1:n
%             D(jj,kk)=D_ind(1,jj,kk);
%         end
%     end
m=floor(0.9*n);          %m ���ϸ���
Alpha=1;       %Alpha ������Ϣ����Ҫ�̶ȵĲ���
Beta=5;        %Beta ��������ʽ������Ҫ�̶ȵĲ���
Rho=0.5;       %Rho ��Ϣ������ϵ��
NC_max=180;    %��������������������Ϊ180���ο����Ŵ��˻��㷨��ִ�е������������Ʊ���
Q=100;         %��Ϣ������ǿ��ϵ��




Eta=1./D;                     %EtaΪ�������ӣ�������Ϊ����ĵ���
Tau=ones(n,n);                %TauΪ��Ϣ�ؾ���
Tabu=zeros(m,n);              %�洢����¼·��������
NC=1;                         %��������������¼��������
R_best=zeros(NC_max,n);       %�������·��
L_best=inf.*ones(NC_max,1);   %�������·�ߵĳ���
L_ave=zeros(NC_max,1);        %����·�ߵ�ƽ������
 
 
 
while NC<=NC_max        %ֹͣ����֮һ���ﵽ������������ֹͣ
    %%2.��mֻ���Ϸŵ�n��������
    Randpos=[];   %�漴��ȡ
    for i=1:(ceil(m/n))
        Randpos=[Randpos,randperm(n)];
    end
    Tabu(:,1)=(Randpos(1,1:m))';   
    %%3.mֻ���ϰ����ʺ���ѡ����һ�����У���ɸ��Ե�����
    for j=2:n     %���ڳ��в�����
        for i=1:m
            visited=Tabu(i,1:(j-1)); %��¼�ѷ��ʵĳ��У������ظ�����
            J=zeros(1,(n-j+1));       %�����ʵĳ���
            P=J;                      %�����ʳ��е�ѡ����ʷֲ�
            Jc=1;
            for k=1:n
                if length(find(visited==k))==0   %��ʼʱ��0
                    J(Jc)=k;
                    Jc=Jc+1;                         %���ʵĳ��и����Լ�1
                end
            end
            %��������ѡ���еĸ��ʷֲ�
            for k=1:length(J)
                P(k)=(Tau(visited(end),J(k))^Alpha)*(Eta(visited(end),J(k))^Beta);
            end
            P=P/(sum(P));
            %������ԭ��ѡȡ��һ������
            Pcum=cumsum(P);     %cumsum��Ԫ���ۼӼ����
            Select=find(Pcum>=rand); %������ĸ��ʴ���ԭ���ľ�ѡ������·��
            to_visit=J(Select(1));
            Tabu(i,j)=to_visit;
        end
    end
    if NC>=2
        Tabu(1,:)=R_best(NC-1,:);
    end
    %%4.��¼���ε������·��
    L=zeros(m,1);     %��ʼ����Ϊ0��m*1��������
    for i=1:m
        R=Tabu(i,:);
        for j=1:(n-1)
            L(i)=L(i)+D(R(j),R(j+1));    %ԭ������ϵ�j�����е���j+1�����еľ���
        end
        L(i)=L(i)+D(R(1),R(n));      %һ���������߹��ľ���
    end
    L_best(NC)=min(L);           %��Ѿ���ȡ��С
    pos=find(L==L_best(NC));
    R_best(NC,:)=Tabu(pos(1),:); %���ֵ���������·��
    L_ave(NC)=mean(L);           %���ֵ������ƽ������
    NC=NC+1;                      %��������
 
 
    %%5.������Ϣ��
    Delta_Tau=zeros(n,n);        %��ʼʱ��Ϣ��Ϊn*n��0����
    for i=1:m
        for j=1:(n-1)
            Delta_Tau(Tabu(i,j),Tabu(i,j+1))=Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
            %�˴�ѭ����·����i��j���ϵ���Ϣ������
        end
        Delta_Tau(Tabu(i,n),Tabu(i,1))=Delta_Tau(Tabu(i,n),Tabu(i,1))+Q/L(i);
        %�˴�ѭ��������·���ϵ���Ϣ������
    end
    Tau=(1-Rho).*Tau+Delta_Tau; %������Ϣ�ػӷ������º����Ϣ��
    %%6.���ɱ�����
    Tabu=zeros(m,n);             %ֱ������������
end
Pos=find(L_best==min(L_best)) ;  %�ҵ����·������0Ϊ�棩
Shortest_Route=R_best(Pos(1),:) ; %���������������·��
Shortest_Length=L_best(Pos(1));   %��������������̾��� 
% toc %��ʱ����
end






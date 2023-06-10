function [Shortest_Length,Shortest_Route] = ACOTSP(D,n)
%% 蚁群算法
% tic;           %计时开始
%%1.初始化
% global D
% global ND
% n=ND; %n表示问题的规模（城市个数）

% [i,j,k]=size(D_ind);
% D=zeros(n,n);
%     for jj=1:n
%         for kk=1:n
%             D(jj,kk)=D_ind(1,jj,kk);
%         end
%     end
m=floor(0.9*n);          %m 蚂蚁个数
Alpha=1;       %Alpha 表征信息素重要程度的参数
Beta=5;        %Beta 表征启发式因子重要程度的参数
Rho=0.5;       %Rho 信息素蒸发系数
NC_max=180;    %最大迭代次数，这里设置为180，参考了遗传退火算法的执行迭代次数，控制变量
Q=100;         %信息素增加强度系数




Eta=1./D;                     %Eta为启发因子，这里设为距离的倒数
Tau=ones(n,n);                %Tau为信息素矩阵
Tabu=zeros(m,n);              %存储并记录路径的生成
NC=1;                         %迭代计数器，记录迭代次数
R_best=zeros(NC_max,n);       %各代最佳路线
L_best=inf.*ones(NC_max,1);   %各代最佳路线的长度
L_ave=zeros(NC_max,1);        %各代路线的平均长度
 
 
 
while NC<=NC_max        %停止条件之一：达到最大迭代次数，停止
    %%2.将m只蚂蚁放到n个城市上
    Randpos=[];   %随即存取
    for i=1:(ceil(m/n))
        Randpos=[Randpos,randperm(n)];
    end
    Tabu(:,1)=(Randpos(1,1:m))';   
    %%3.m只蚂蚁按概率函数选择下一座城市，完成各自的周游
    for j=2:n     %所在城市不计算
        for i=1:m
            visited=Tabu(i,1:(j-1)); %记录已访问的城市，避免重复访问
            J=zeros(1,(n-j+1));       %待访问的城市
            P=J;                      %待访问城市的选择概率分布
            Jc=1;
            for k=1:n
                if length(find(visited==k))==0   %开始时置0
                    J(Jc)=k;
                    Jc=Jc+1;                         %访问的城市个数自加1
                end
            end
            %下面计算待选城市的概率分布
            for k=1:length(J)
                P(k)=(Tau(visited(end),J(k))^Alpha)*(Eta(visited(end),J(k))^Beta);
            end
            P=P/(sum(P));
            %按概率原则选取下一个城市
            Pcum=cumsum(P);     %cumsum，元素累加即求和
            Select=find(Pcum>=rand); %若计算的概率大于原来的就选择这条路线
            to_visit=J(Select(1));
            Tabu(i,j)=to_visit;
        end
    end
    if NC>=2
        Tabu(1,:)=R_best(NC-1,:);
    end
    %%4.记录本次迭代最佳路线
    L=zeros(m,1);     %开始距离为0，m*1的列向量
    for i=1:m
        R=Tabu(i,:);
        for j=1:(n-1)
            L(i)=L(i)+D(R(j),R(j+1));    %原距离加上第j个城市到第j+1个城市的距离
        end
        L(i)=L(i)+D(R(1),R(n));      %一轮下来后走过的距离
    end
    L_best(NC)=min(L);           %最佳距离取最小
    pos=find(L==L_best(NC));
    R_best(NC,:)=Tabu(pos(1),:); %此轮迭代后的最佳路线
    L_ave(NC)=mean(L);           %此轮迭代后的平均距离
    NC=NC+1;                      %迭代继续
 
 
    %%5.更新信息素
    Delta_Tau=zeros(n,n);        %开始时信息素为n*n的0矩阵
    for i=1:m
        for j=1:(n-1)
            Delta_Tau(Tabu(i,j),Tabu(i,j+1))=Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
            %此次循环在路径（i，j）上的信息素增量
        end
        Delta_Tau(Tabu(i,n),Tabu(i,1))=Delta_Tau(Tabu(i,n),Tabu(i,1))+Q/L(i);
        %此次循环在整个路径上的信息素增量
    end
    Tau=(1-Rho).*Tau+Delta_Tau; %考虑信息素挥发，更新后的信息素
    %%6.禁忌表清零
    Tabu=zeros(m,n);             %直到最大迭代次数
end
Pos=find(L_best==min(L_best)) ;  %找到最佳路径（非0为真）
Shortest_Route=R_best(Pos(1),:) ; %最大迭代次数后最佳路径
Shortest_Length=L_best(Pos(1));   %最大迭代次数后最短距离 
% toc %计时结束
end






%%
clc;clear;close all;

%% parameter
% f=30*10^3;
%% NL
s=0.5;
w=0;
 r=0:100:3500;
% r=100000;
% f=18000:100:34000;
f=26000;
% f=0:100:1000000;
B=16000;
% P0=1;
% Prx=0.75;
DI=0;
% BER=10^(-6);
% SNR=10*log10((qfuncinv(BER))^2/2);
%%%%%%%SNR=10;
pt=65;
H=100;
for i=1:length(f)
    %% noise
    n1(i)=17-30*log10(f(i)*10^(-3));%湍流
    n2(i)=40+20*(s-0.5)+26*log10(f(i)*10^(-3))-60*log10(f(i)*10^(-3)+0.03);%
    n3(i)=50+7.5*sqrt(w)+20*log10(f(i)*10^(-3))-40*log10(f(i)*10^(-3)+0.4);
    n4(i)=-15+20*log10(f(i)*10^(-3));
    NL(i)=power(10,0.1*n1(i))+power(10,0.1*n2(i))+power(10,0.1*n3(i))+power(10,0.1*n4(i));
    NLDB(i)=10*log10( NL(i));
    
     NL1DB(i)=50-18*log10(f(i)*10^(-3));
     NL1(i)=10^(NL1DB(i)/10);
     
% %     NL_w(i)=10^((NL(i)-170.77)/10);
    % TL
   
    %n=1;%柱面波
% %     SL=10*log10(pt)+170.77-10*log10(B);
% %      P=2*pi*H*10^(SL/10)*0.67*10^(-18);
    
    %% SNR
    
    %SNR=SL-TL-NL+DI;
     a(i)=0.11*(f(i)*10^(-3))^2/(1+(f(i)*10^(-3))^2)+44*(f(i)*10^(-3))^2/(4100+(f(i)*10^(-3))^2)+2.75*10^(-4)*(f(i)*10^(-3))^2+0.003;
    for j=1:length(r)
% %         TL(i,j)=1.5*10*log10(r(j))+r(j)*10^(-3)*a(i);%dB
% %         A(i,j)=10^(TL(i,j)/10);
% %          A1(i,j)=power((r(j)),1.5)*power(a2(i),r(j)*10^(-3));
% %         Gain(i,j)=1/A(i,j);
% %         Gaindb(i,j)=10*log10(1/A(i,j));
         A1(i,j)=(r(j)^1.5)*(10^(a(i)/10))^(r(j)/1000);
        %%%%%SL(i,j)=SNR-DI+NL(i)+TL(i,j);
        %p(i,j)=10^((A(i,j)-170.77)/10);
        
  
      Prx(i,j)=pt/A1(i,j);
      Prxdb(i,j)=10*log10(Prx(i,j));
% %            snr0(i,j)=1/A(i,j)/NL(i);
% %            snr0db(i,j)=10*log10(snr0(i,j));
          snr(i,j)= Prx(i,j)/(NL1(i)*10^(-12)*B);
        snrdb(i,j)=10*log10(snr(i,j));
        
% %         SNRDB1(i,j)=SL-TL(i,j)-NLDB(i);
% %          SNR1(i,j)=10^(SNRDB1(i,j)/10);
% %          SNR(i,j)=10*log10(P)-10*log(2*pi*H*0.67*10^(-18))-TL(i,j)-NL1(i);
% %          SNRDB(i,j)=10*log10(SNR(i,j));
          pe(i,j)=qfunc(sqrt(2*snr(i,j)));
% %         pe(i,j)= erfc(sqrt(2*snr(i,j))/sqrt(2))/2;
%             pe(i,j)=qfunc(sqrt(4*snr(i,j))*sin(pi/8));
% %          pb(i,j)=1/2-1/2*sqrt(10^(snr(i,j)/10)/1+10^(snr(i,j)/10));
% %     pe(i,j)=1/(4*snr(i,j));
       pb(i,j)=1-(1-pe(i,j))^100;
    end
end
% figure;
% surf(Ptx);
% figure;
% axis ([21 27 0 5000 0 25]);
% surf(r,f/1000,pe);

% x1=xlabel('Frequency (kHz)','Fontname', 'Times New Roman','FontSize',12);        %x轴标题'Fontname', 'Times New Roman',,'FontSize',12
% x2=ylabel('Distance (m)','Fontname', 'Times New Roman','FontSize',12);        %y轴标题,'FontSize',12
% x3=zlabel('Power (w)','Fontname', 'Times New Roman','FontSize',12);        %z轴标题
% set(x1,'Rotation',30);    %x轴名称旋转
% set(x2,'Rotation',-30);    %y轴名称旋转
clr(1,:)=[0 113 188]/255;
clr(2,:)=[216 82 24]/255;
clr(3,:)=[236 176 31]/255;
clr(4,:)=[125 46 141]/255;
clr(5,:)=[118 171 47]/255;
clr(6,:)=[76 189 237]/255;
clr(7,:)=[161 19 46]/255;
clr(8,:)=[255 255 255]/255;
clr(9,:)=[0 0 0]/255;
% 
% figure
% plot(r,snrdb);
% grid on;
% figure
% semilogy(snrdb,pe1);



figure
semilogy(snrdb,pe);
% % hold on
figure
semilogy(r,pb);

figure
 yyaxis left
 plot(r,pe,'*-','LineWidth',1,'MarkerSize',6);
 xlabel('\fontsize{12} distance(m)');
ylabel('\fontsize{12} p_e ');
 yyaxis right
 plot(r,snrdb,'o--','LineWidth',1,'MarkerSize',6); 
backColor = [245 249 253]/255;
set(gca, 'color', backColor);
grid on;
grid minor
% set(gca,'xtick',1:1:10);
xlabel('\fontsize{12} distance(m)');
ylabel('\fontsize{12} SNR ');

hold on;
axes('Position',[0.15,0.2,0.31,0.2]); % 生成子图   左右  上下 宽窄
plot(r,pe,'*-','LineWidth',1,'MarkerSize',6);                                                                                                         
xlim([0,2200]); % 设置坐标轴范围 
grid on;
grid minor
annotation('rectangle',[0.131,0.1,0.5,0.03],'LineStyle','-','Color','g','LineWidth',0.7)
% annotation('line',[0.1,0.13],[0.35,0.15],'LineStyle','-','Color','g','LineWidth',0.7)
annotation('arrow',[0.15,0.13],[0.2,0.13],'LineStyle','-','Color','g','LineWidth',0.7)
% annotation('line',[0.58,0.46],[0.35,0.15],'LineStyle','-','Color','g','LineWidth',0.7)
annotation('arrow',[0.46,0.63],[0.2,0.13],'LineStyle','-','Color','g','LineWidth',0.7)
%axis([0 0.1 0 35]);

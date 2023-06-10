function [packeterrorate] = packeterroratefunction(r,L_packet)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

%% parameter
% % f=30*10^3;
% % NL
s=0.5;
w=0;

%%distance
%r=0:100:3500;
% % r=100000;
% % f=18000:100:34000;
f=26000;
% % f=0:100:1000000;
B=16000;
% % P0=1;
% % Prx=0.75;
DI=0;
% % BER=10^(-6);
% % SNR=10*log10((qfuncinv(BER))^2/2);
% %%%%%%%SNR=10;
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
    
    % SNR
    
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
        pb(i,j)=1-(1-pe(i,j))^L_packet;
        packeterrorate=pb;
    end
end
end


%% Cooperative probability Vs Network size
clear;clc;close all;set(0,'DefaultAxesFontSize',13);
N=1024; M=round(0.25*N); a=0.7; bk=0.65; 
Psi=fft(eye(N))/sqrt(N);     rng(0)
Phi=Phi_selection(N,M,2,2e7);
[x,~,~,~,K,index,~] = WFG_Re(9e8,2e7,N,20); close figure 1
A=Phi/Psi; At=A';  Y=Phi*x;  
L=1:20; SNR=-5:5:5;  MC=10;   % parameters
l=length(L);    snr=length(SNR);    
Xf_mat=zeros(N,2*l+1); 
for s=1:snr
    tic;
disp(['SNR(dB)', ' = ', num2str(SNR(s))]);
Pd=zeros(2,l);Pfa=Pd;
for j=1:l
for frame=1:MC
y = awgn(repmat(Y, 1, L(j)), SNR(s), 'measured');
YL = Virtual(y);
for i=1:2*L(j)+1
	Xf_mat(:,i)=abs(fct_RMP(A,At,YL(:,i),K,bk,a));
end
xf=min(Xf_mat(:,L(j)+1:2*L(j)), [], 2);xf=(fftshift(xf)>0);
Pd(1,j)=Pd(1,j)+sum(xf(index)>0)/K;
Pfa(1,j)=Pfa(1,j)+(sum(xf>0)-sum(xf(index)>0))/(N-K);
xf=min(Xf_mat(:,1:L(j)), [], 2);xf=(fftshift(xf)>0);
Pd(2,j)=Pd(2,j)+sum(xf(index)>0)/K;
Pfa(2,j)=Pfa(2,j)+(sum(xf>0)-sum(xf(index)>0))/(N-K);
end
end
Pfa=Pfa/MC;Pd=Pd/MC;
subplot(snr,2,2*s-1);hold on; grid on;
xlabel('Network Size');ylabel('Q_d');
plot(L,Pd(1,:),'b','LineWidth',1.5);
[a1,b1,c1] = fit_logistic_inc(L,Pd(1,:));%fit Qd Virt
plot(L,Pd(2,:),'b--','LineWidth',1.5);
[a2,b2,c2] = fit_logistic_dec(L,Pd(2,:));
legend('Q_d Virt',sprintf('Logistic: Q_d = %.2f / (1 + exp(-%.2f(L - %.2f)))',a1, b1, c1),'Q_{fa} Real',...
sprintf('Logistic inv: Q_d = %.2f / (1 + exp(%.2f(%.2f+L )))', a2, b2, -c2),...
'Q_d Real','Location', 'best')
subplot(snr,2,2*s);hold on; grid on;xlabel('Network Size');
ylabel('Q_{fa}');plot(L,Pfa(1,:),'r','LineWidth',1.5);plot(L,Pfa(2,:),'r--','LineWidth',1.5);
legend('Q_{fa} Virt','Location', 'best')
titleStr = sprintf('SNR(dB) = %d', SNR(s));title(titleStr);toc
disp(['Date & Time', ' = ', num2str(datestr(now))]);
end
%global_title = sprintf('MC: %d, SNR(dB): %d, N: %d, M: %d, K: %d',MC,SNR(s),N,M,K);
%title(global_title);
%tightfig;
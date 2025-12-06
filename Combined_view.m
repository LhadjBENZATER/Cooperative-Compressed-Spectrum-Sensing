clear;clc;close all;set(0,'DefaultAxesFontSize',13);
N=1024; M=round(0.25*N); fs =20e6; fc=9e8;  pfa=0.07;
frames=50;  rng(0); a_cor=0.7;   bk=0.65; %a=0.5-0.7; bk=0.15-0.75
Phi=Phi_selection(N,M,2,fs);
SigMatr=zeros(N,9);
Pfa=zeros(13,1); Pd=Pfa;
Psi=fft(eye(N));    A=Phi/Psi; At=A';
[x, xf,t,f,K,ind] = WFG_Re(fc,fs,N,20);
%[x,xf,t,f,K,ind] = WFG_BPSK(f1,f2,N0,fc,fs,N);
close figure 1
V8=0.25;  V9=0.25;    V4=0.25;  SNR_y=5;  th=0.02;

%% Save Doub_to_bin_DAT(real(x),'x.dat'); Already saved
 %% Original spectrum
[threshold_1,~]=CA_CFAR1(pfa,32,16, xf);offset=0.2;
Pd(1,1)=sum(xf(ind)>threshold_1(ind)+offset)/K;
Pfa(1,1)=(sum(xf>threshold_1+offset)-sum(xf(ind)>threshold_1(ind))+offset)/(N-K);
for frame=1:frames
Y = awgn(repmat(Phi*x, 1, 4), SNR_y, 'measured');
y = Virtual(Y);
%% %%%%%%%%%%%%%% Recovery   %%%%%%%%%%%%%
    for i=1:9
    SigMatr(:,i)=abs(fct_RMP(A,At,y(:,i),K,bk,a_cor));
    end
%% Averaged channel Y=mean(Y_i)
xf_2=fftshift(SigMatr(:,9));% xf_2 =medfilt1(xf_2,2);
xf_2=xf_2/max(xf_2);    thr_AV = 0.2*ones(N,1);
Pd(2,1)=Pd(2,1)+sum(xf_2(ind)>thr_AV(ind))/K;
Pfa(2,1)=Pfa(2,1)+(sum(xf_2>thr_AV)-sum(xf_2(ind)>thr_AV(ind)))/(N-K);
%%  AND Rule
% 4 Virtual Channels
xf_3=min(SigMatr(:,5:8), [], 2);xf_3=(fftshift(xf_3)>0);
Pd(3,1)=Pd(3,1)+sum(xf_3(ind))/K;
Pfa(3,1)=Pfa(3,1)+(sum(xf_3)-sum(xf_3(ind)))/(N-K);
% 4 Real Channels
xf_4=min(SigMatr(:,1:4), [], 2);xf_4=(fftshift(xf_4)>0);
Pd(4,1)=Pd(4,1)+sum(xf_4(ind))/K;
Pfa(4,1)=Pfa(4,1)+(sum(xf_4)-sum(xf_4(ind)))/(N-K);

%% Voting 
% Virtual
xf_5=SigMatr(:,5:8); Nz = xf_5 ~= 0;
xf_5(Nz) = 1; xf_5=0.25*sum(xf_5, 2);
xf_5 =medfilt1(fftshift(xf_5),2);  
Pd(5,1)=Pd(5,1)+sum(xf_5(ind)>V4)/K;
Pfa(5,1)=Pfa(5,1)+(sum(xf_5>V4)-sum(xf_5(ind)>V4))/(N-K);
% Real
xf_6=SigMatr(:,1:4); Nz = xf_6 ~= 0;xf_6(Nz) = 1; 
xf_6=0.25*sum(xf_6, 2);xf_6 =medfilt1(fftshift(xf_6),2);
Pd(6,1)=Pd(6,1)+sum(xf_6(ind)>V4)/K;
Pfa(6,1)=Pfa(6,1)+(sum(xf_6>V4)-sum(xf_6(ind)>V4))/(N-K);
% 4 Real +4 Virt Channels
xf_10=SigMatr(:,1:8); Nz = xf_10 ~= 0;
xf_10(Nz) = 1; xf_10=sum(xf_10, 2)/8;
xf_10 =medfilt1(fftshift(xf_10),2);
Pd(10,1)=Pd(10,1)+sum(xf_10(ind)>V8)/K;
Pfa(10,1)=Pfa(10,1)+(sum(xf_10>V8)-sum(xf_10(ind)>V8))/(N-K);
% 4 Real +5 Virt Channels
xf_12=SigMatr; Nz = xf_12 ~= 0;xf_12(Nz) = 1; xf_12=sum(xf_12, 2)/9;
xf_12 =medfilt1(fftshift(xf_12),2);
Pd(12,1)=Pd(12,1)+sum(xf_12(ind)>V9)/K;
Pfa(12,1)=Pfa(12,1)+(sum(xf_12>V9)-sum(xf_12(ind)>V9))/(N-K);
%%  Average spectrum
% Virtual Channels
xf_7=mean(SigMatr(:,5:8),2);xf_7=abs(xf_7/max(xf_7));
xf_7 =medfilt1(fftshift(xf_7),2); xf_7=xf_7.^2;
Pd(7,1)=Pd(7,1)+sum(xf_7(ind)>th)/K;
Pfa(7,1)=Pfa(7,1)+(sum(xf_7>th)-sum(xf_7(ind)>th))/(N-K);
%4 Real Channels
xf_8=mean(SigMatr(:,1:4),2);xf_8=abs(xf_8/max(xf_8));
xf_8 =medfilt1(fftshift(xf_8),2);   xf_8=xf_8.^2;
Pd(8,1)=Pd(8,1)+sum(xf_8(ind)>th)/K;
Pfa(8,1)=Pfa(8,1)+(sum(xf_8>th)-sum(xf_8(ind)>th))/(N-K);
% 4 Real +4 Virt Channels
xf_11=mean(SigMatr(:,1:8),2);xf_11=abs(xf_11/max(xf_11));
xf_11 =medfilt1(fftshift(xf_11),2); xf_11=xf_11.^2;
Pd(11,1)=Pd(11,1)+sum(xf_11(ind)>th)/K;
Pfa(11,1)=Pfa(11,1)+(sum(xf_11>th)-sum(xf_11(ind)>th))/(N-K);
% 4 Real +5 Virt Channels
xf_13=mean(SigMatr,2);xf_13=abs(xf_13/max(xf_13));
xf_13 =medfilt1(fftshift(xf_13),2); xf_13=xf_13.^2;
Pd(13,1)=Pd(13,1)+sum(xf_13(ind)>th)/K;
Pfa(13,1)=Pfa(13,1)+(sum(xf_13>th)-sum(xf_13(ind)>th))/(N-K);
end
Pfa(2:end,1)=Pfa(2:end,1)/frames;Pd(2:end,1)=Pd(2:end,1)/frames;
subplot(4,1,1); grid on;hold on;
plot(f,xf,'LineWidth',1.5); plot(f,threshold_1, 'LineWidth',1.5);
plot(f,xf_2, 'LineWidth',1.5); plot(f,thr_AV, 'LineWidth',1.5);
l1 = sprintf('Pd=%.2f%%, Pfa=%.2f%%', 100*Pd(1,1), 100*Pfa(1,1));
l2 = sprintf('Pd=%.2f%%, Pfa=%.2f%%', 100*Pd(2,1), 100*Pfa(2,1));
legend({'Original Spectrum',l1,'Avg y_i',l2});
xlabel('Freq (MHz)');ylabel('Mag (V/Hz)');
titleText = sprintf(['X(f)',', N= ', num2str(N), ', M=', num2str(M), ', K=', num2str(K),', f_s=', num2str(fs/(1e6)),'MHz, f_c=', num2str(fc/(1e6)),'MHz, Frames=', num2str(frames),'   SNR(dB)=', num2str(SNR_y) ]);
title(titleText,'FontSize',14);

subplot(4,2,4);grid on;hold on
plot(f,xf_3); plot(f,th*ones(N,1), 'm--','LineWidth',1);
legend(sprintf('Pd=%.2f%%, Pfa=%.2f%%', 100*Pd(3,1), 100*Pfa(3,1)))
set(gca, 'xticklabels', []);
title('$|\mathbf{\hat{X}(f)}|=AND(|\mathbf{\hat{X}_i}|) \ Virt$', 'Interpreter','latex','FontSize',13);

subplot(4,2,3);grid on;hold on
plot(f,xf_4); plot(f,th*ones(N,1), 'm--','LineWidth',1);
legend(sprintf('Pd=%.2f%%, Pfa=%.2f%%', 100*Pd(4,1), 100*Pfa(4,1)),...
'Threshold');set(gca, 'xticklabels', []);
title('$|\mathbf{\hat{X}(f)}|=AND(|\mathbf{\hat{X}_i}|) \ Real$', 'Interpreter','latex','FontSize',13);



subplot(4,4,10); grid on;hold on;
plot(f,xf_5,'LineWidth',1); plot(f,V4*ones(N,1),'m--','LineWidth',1);
legend(sprintf('Pd=%.2f%%, Pfa=%.2f%%', 100*Pd(5,1), 100*Pfa(5,1)),...
'Threshold');set(gca, 'xticklabels', []);
title('$|\mathbf{\hat{X}(f)}|=Voting(|\mathbf{\hat{X}_i}|) \ Virt$', 'Interpreter','latex','FontSize',13);

subplot(4,4,9); grid on;hold on;
plot(f,xf_6,'LineWidth',1);plot(f,V4*ones(N,1), 'm--','LineWidth',1);
legend(sprintf('Pd=%.2f%%, Pfa=%.2f%%', 100*Pd(6,1), 100*Pfa(6,1)),...
'Threshold');set(gca, 'xticklabels', []);
title('$|\mathbf{\hat{X}(f)}|=Voting(|\mathbf{\hat{X}_i}|) \Re$', 'Interpreter','latex','FontSize',13);

subplot(4,4,11); grid on;hold on;
plot(f,xf_10,'LineWidth',1);plot(f,V8*ones(N,1), 'm--','LineWidth',1);
legend(sprintf('Pd=%.2f%%, Pfa=%.2f%%', 100*Pd(10,1), 100*Pfa(10,1)),...
'Threshold');set(gca, 'xticklabels', []);
title('$|\mathbf{\hat{X}(f)}|=Voting(|\mathbf{\hat{X}_i}|) \ 8ch$', 'Interpreter','latex','FontSize',13);

subplot(4,4,12); grid on;hold on;
plot(f,xf_12,'LineWidth',1);plot(f,V9*ones(N,1), 'm--','LineWidth',1);
legend(sprintf('Pd=%.2f%%, Pfa=%.2f%%', 100*Pd(12,1), 100*Pfa(12,1)),...
'Threshold');set(gca, 'xticklabels', []);
title('$|\mathbf{\hat{X}(f)}|=Voting(|\mathbf{\hat{X}_i}|) \ 9ch$', 'Interpreter','latex','FontSize',13);

subplot(4,4,14);grid on;hold on
plot(f,xf_7); plot(f,th*ones(N,1), 'm--','LineWidth',1);
legend(sprintf('Pd=%.2f%%, Pfa=%.2f%%', 100*Pd(7,1), 100*Pfa(7,1)),...
'Threshold');%set(gca, 'xticklabels', []);
title('$|\mathbf{\hat{X}(f)}|=Average (|\mathbf{\hat{X}_i}|) \ Virt$', 'Interpreter','latex','FontSize',13);

subplot(4,4,13);grid on;hold on
plot(f,xf_8); plot(f,th*ones(N,1), 'm--','LineWidth',1);
legend(sprintf('Pd=%.2f%%, Pfa=%.2f%%', 100*Pd(8,1), 100*Pfa(8,1)),...
'Threshold');
title('$|\mathbf{\hat{X}(f)}|=Average (|\mathbf{\hat{X}_i}|) \Re$', 'Interpreter','latex','FontSize',13);

subplot(4,4,15);grid on;hold on
plot(f,xf_11); plot(f,th*ones(N,1), 'm--','LineWidth',1);
legend(sprintf('Pd=%.2f%%, Pfa=%.2f%%', 100*Pd(11,1), 100*Pfa(11,1)),...
'Threshold');
title('$|\mathbf{\hat{X}(f)}|=Average (|\mathbf{\hat{X}_i}|) \ 8ch$', 'Interpreter','latex','FontSize',13);

subplot(4,4,16);grid on;hold on
plot(f,xf_13); plot(f,th*ones(N,1), 'm--','LineWidth',1);
legend(sprintf('Pd=%.2f%%, Pfa=%.2f%%', 100*Pd(13,1), 100*Pfa(13,1)),...
'Threshold');
title('$|\mathbf{\hat{X}(f)}|=Average (|\mathbf{\hat{X}_i}|) \ 9ch$', 'Interpreter','latex','FontSize',13);
%tightfig;
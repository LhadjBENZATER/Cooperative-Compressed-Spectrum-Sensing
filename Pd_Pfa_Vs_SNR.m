% Greddy CSS performance
clear;clc;close all;set(0,'DefaultAxesFontSize',13);
N=1024; M=round(0.25*N); a=0.7;   bk=0.65; 
SNR=-36:3:30;   MC=50; thr=0.04; rng(0)
Psi=fft(eye(N))/sqrt(N);Phi=Phi_selection(N,M,2,2e7);
[x,~,~,~,K,index,~] = WFG_Re(9e8,2e7,N,20); close figure 1
A=Phi/Psi; At=A';  Y=Phi*x;       
L=[1 4 8];V=[1;1/4;5/16]; Xf_mat=zeros(N,length(SNR));

%% CSS with empty Support
disp('CSS with RMP');
for l=1:length(L)
disp(['Network size', ' = ', num2str(L(l))]);
Pfa=zeros(10,length(SNR));Pd=Pfa;   
Sig=zeros(N,2*L(l)+1); tic
for snr=1:length(SNR)
    disp(['SNR = ', num2str(SNR(snr))]);
for frame=1:MC
y = awgn(repmat(Y, 1, L(l)), SNR(snr), 'measured');
YL = Virtual(y);
for i=1:2*L(l)+1
	Sig(:,i)=abs(fct_RMP(A,At,YL(:,i),K,bk,a));
end
[Pd,Pfa] = fct_Cooperat(Sig,V,index,N,K,L,Pd,Pfa,thr,snr,l);
end
end
Pfa=Pfa/MC;Pd=Pd/MC;    fct_plot(SNR,Pd,Pfa,L,l);
disp(['Time(min) = ', num2str(toc/60)]);
end
global_title = sprintf('Without partial Supp, MC: %d, N: %d, M: %d, K: %d',MC,N,M,K);
title(global_title)
legend('Proposed Supp','Supp ','Vote Virt n','Reference','Vote Virt n+1','Vote All','Avg X_i(f) Virt n','Avg X_i(f)','Avg X_i(f) Virt n+1','Avg X_i(f) All','P_{fa}=0.1','Location', 'best')
tightfig;
clear;clc;close all;set(0,'DefaultAxesFontSize',13);
N=1024;  MC=20;SNR=-37:0.5:30; M=256;th=0.3;
Pfa=zeros(6,length(SNR));Pd=Pfa; threshold=zeros(N,length(SNR)); 
SigMatr=zeros(N,9); Psi=fft(eye(N))/sqrt(N);
[x,~,~,~,K,index,~] = WFG_Re(9e8,2e7,N,30); close figure 1
Phi=Phi_selection(N,M,2,2e7);A=Phi/Psi;At=A';
alpha=6;
for snr=1:length(SNR)
disp(['SNR(dB)=', num2str(SNR(snr))]);    
    for mc=1:MC
xL = Virtual(awgn(repmat(x, 1, 4), SNR(snr), 'measured'));
y=Phi*xL;

for i=1:4:9
	SigMatr(:,i)=abs(fct_RMP(A,At, y(:,i), K, 0.65, 0.7));
end
xf=fftshift(SigMatr(:,1));xf =medfilt1(xf,2);xf=xf/max(xf);    
Pd(4,snr)=Pd(4,snr)+sum(xf(index)>th)/K;
Pfa(4,snr)=Pfa(4,snr)+(sum(xf>th)-sum(xf(index)>th))/(N-K);

xf=fftshift(SigMatr(:,9));xf =medfilt1(xf,2);xf=xf/max(xf);    
Pd(6,snr)=Pd(6,snr)+sum(xf(index)>th)/K;
Pfa(5,snr)=Pfa(6,snr)+(sum(xf>th)-sum(xf(index)>th))/(N-K);
    end
end
Pd=Pd/MC; Pfa=Pfa/MC;
subplot(1,2,1);hold on; grid on;
plot(SNR+alpha,smooth(Pd(1,:)),'b','LineWidth',1.5);
plot(SNR+alpha,smooth(Pd(2,:)),'g','LineWidth',1.5);
plot(SNR+alpha,smooth(Pd(3,:)),'m','LineWidth',1.5);
plot(SNR,0.9*ones(1,length(SNR)),'r','LineWidth',2)
plot(SNR,smooth(Pd(4,:)),'b--','LineWidth',1.5);
plot(SNR,smooth(Pd(5,:)),'g--','LineWidth',1.5);
plot(SNR,smooth(Pd(6,:)),'m--','LineWidth',1.5);
xlabel('SNR (dB)');ylabel('P_{d}');
legend('No CSS Re','No CSS Virt','No CSS Avg','P_{d}=0.9','CSS Re','CSS Virt','CSS Avg')
legend('No CSS','No CSS Avg','P_{fa}=0.1','CSS','CSS Avg')
subplot(1,2,2);grid on; hold on;
plot(SNR+alpha,smooth(Pfa(1,:)),'b','LineWidth',1.5);
plot(SNR+alpha,smooth(Pfa(2,:)),'g','LineWidth',1.5);
plot(SNR+alpha,smooth(Pfa(3,:)),'m','LineWidth',1.5);
plot(SNR,0.1*ones(1,length(SNR)),'r','LineWidth',2)
plot(SNR,smooth(Pfa(4,:)),'b:','LineWidth',1.5);
plot(SNR,smooth(Pfa(5,:)),'g:','LineWidth',1.5);
plot(SNR,smooth(Pfa(6,:)),'m:','LineWidth',1.5);
xlabel('SNR (dB)');ylabel('P_{fa}');
legend('No CSS Re','No CSS Virt','No CSS Avg','P_{fa}=0.1','CSS Re','CSS Virt','CSS Avg')

Pfa(4,:)=Pfa(4,:)/2;
hold on
plot(SNR,smooth(Pfa(4,:)),'b--','LineWidth',1.5);
tightfig;
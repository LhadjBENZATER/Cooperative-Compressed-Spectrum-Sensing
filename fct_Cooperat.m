function [Pd,Pfa] = fct_Cooperat(Sig,V,index,N,K,L,Pd,Pfa,thr,snr,l)
%% Voting 
% Virtual n
xf=Sig(:,L(l)+1:2*L(l)); Nz=xf~=0;xf(Nz)=1; xf=sum(xf, 2);
xf=abs(xf/max(xf)); xf =medfilt1(fftshift(xf),2);  
Pd(3,snr)=Pd(3,snr)+sum(xf(index)>V(l,1))/K;
Pfa(3,snr)=Pfa(3,snr)+(sum(xf>V(l,1))-sum(xf(index)>V(l,1)))/(N-K);
% % Real
xf=Sig(:,1:L(l)); Nz=xf~= 0;xf(Nz)=1; xf=sum(xf,2);
xf=abs(xf/max(xf));xf =medfilt1(fftshift(xf),2);
Pd(4,snr)=Pd(4,snr)+sum(xf(index)>V(l,1))/K;
Pfa(4,snr)=Pfa(4,snr)+(sum(xf>V(l,1))-sum(xf(index)>V(l,1)))/(N-K);
% Virtual n+1
xf=Sig(:,L(l)+1:2*L(l)+1); Nz=xf~=0;xf(Nz)=1; xf=sum(xf, 2);
xf=abs(xf/max(xf)); xf =medfilt1(fftshift(xf),2);  
Pd(5,snr)=Pd(5,snr)+sum(xf(index)>V(l,1))/K;
Pfa(5,snr)=Pfa(5,snr)+(sum(xf>V(l,1))-sum(xf(index)>V(l,1)))/(N-K);
%  All
xf=Sig; Nz=xf~=0;xf(Nz)=1; xf=sum(xf, 2);
xf=abs(xf/max(xf)); xf =medfilt1(fftshift(xf),2);  
Pd(6,snr)=Pd(6,snr)+sum(xf(index)>V(l,1))/K;
Pfa(6,snr)=Pfa(6,snr)+(sum(xf>V(l,1))-sum(xf(index)>V(l,1)))/(N-K);
%% Average spectrum
% Virtual Channels n
xf=mean(Sig(:,L(l)+1:2*L(l)),2);xf=abs(xf/max(xf));
xf =medfilt1(fftshift(xf),2);   xf=xf.^2;
Pd(7,snr)=Pd(7,snr)+sum(xf(index)>thr)/K;
Pfa(7,snr)=Pfa(7,snr)+(sum(xf>thr)-sum(xf(index)>thr))/(N-K);
%4 Real Channels
xf=mean(Sig(:,1:L(l)),2);xf=abs(xf/max(xf));
xf =medfilt1(fftshift(xf),2);   xf=xf.^2;
Pd(8,snr)=Pd(8,snr)+sum(xf(index)>thr)/K;
Pfa(8,snr)=Pfa(8,snr)+(sum(xf>thr)-sum(xf(index)>thr))/(N-K);
% Virtual Channels n+1
xf=mean(Sig(:,L(l)+1:2*L(l)+1),2);xf=abs(xf/max(xf));
xf =medfilt1(fftshift(xf),2);   xf=xf.^2;
Pd(9,snr)=Pd(9,snr)+sum(xf(index)>thr)/K;
Pfa(9,snr)=Pfa(9,snr)+(sum(xf>thr)-sum(xf(index)>thr))/(N-K);
% All Channels
xf=mean(Sig,2);xf=abs(xf/max(xf));
xf =medfilt1(fftshift(xf),2);   xf=xf.^2;
Pd(10,snr)=Pd(10,snr)+sum(xf(index)>thr)/K;
Pfa(10,snr)=Pfa(10,snr)+(sum(xf>thr)-sum(xf(index)>thr))/(N-K);
%% AND Rule
%4 Virtual Channels
xf=min(Sig(:,L(l)+1:2*L(l)), [], 2);xf=(fftshift(xf)>0);
Pd(1,snr)=Pd(1,snr)+sum(xf(index))/K;
Pfa(1,snr)=Pfa(1,snr)+(sum(xf)-sum(xf(index)))/(N-K);
% 4 Real Channels
xf=min(Sig(:,1:L(l)), [], 2);xf=(fftshift(xf)>0);
Pd(2,snr)=Pd(2,snr)+sum(xf(index))/K;
Pfa(2,snr)=Pfa(2,snr)+(sum(xf)-sum(xf(index)))/(N-K);

end
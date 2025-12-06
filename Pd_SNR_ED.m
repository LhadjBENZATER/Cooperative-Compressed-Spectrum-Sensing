function Pd_simg = Pd_SNR_ED(SNR,Pfa)
var=1;sig=sqrt(var); n=1e6; 
Vt=sqrt((2*sig^2)*log(1/Pfa)); % threshold according to Mahefza chapter single pulse Detection
Pd_simg=zeros(length(SNR),1);
for i=1:length(SNR)
    A=sqrt(2)*sig*10^(SNR(i)/20);
    vi=sig*randn(n,1);
    vq=sig*randn(n,1);
    signal=sqrt((A+vi).^2+vq.^2); % signal according to Mahefza chapter single Pulse Detection
    Pd_simg(i)=mean(signal>Vt);
end


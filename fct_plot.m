function  fct_plot(SNR,Pd,Pfa,L,l)
subplot(length(L),2,2*l-1);hold on; grid on;
plot(SNR,Pd(1,:),'b','LineWidth',1.5);plot(SNR,Pd(2,:),'b--','LineWidth',1.5);
plot(SNR,Pd(3,:),'k','LineWidth',1.5);plot(SNR,Pd(4,:),'k--','LineWidth',1.5);
plot(SNR,Pd(5,:),'k:','LineWidth',1.5);plot(SNR,Pd(6,:),'k-.','LineWidth',1.5);
plot(SNR,Pd(7,:),'m','LineWidth',1.5);plot(SNR,Pd(8,:),'m--','LineWidth',1.5);
plot(SNR,Pd(9,:),'m:','LineWidth',1.5);plot(SNR,Pd(10,:),'m-.','LineWidth',1.5);
plot(SNR,0.9*ones(1,length(SNR)),'r','LineWidth',2)
xlabel('SNR (dB)');ylabel('P_{d}');
title([num2str(L(l)),' Channels']);

subplot(length(L),2,2*l);grid on; hold on;
plot(SNR,Pfa(1,:),'b','LineWidth',1.5);plot(SNR,Pfa(2,:),'b--','LineWidth',1.5);
plot(SNR,Pfa(3,:),'k','LineWidth',1.5);plot(SNR,Pfa(4,:),'k--','LineWidth',1.5);
plot(SNR,Pfa(5,:),'k:','LineWidth',1.5);plot(SNR,Pfa(6,:),'k-.','LineWidth',1.5);
plot(SNR,Pfa(7,:),'m','LineWidth',1.5);plot(SNR,Pfa(8,:),'m--','LineWidth',1.5);
plot(SNR,Pfa(9,:),'m:','LineWidth',1.5);plot(SNR,Pfa(10,:),'m-.','LineWidth',1.5);
plot(SNR,0.1*ones(1,length(SNR)),'r','LineWidth',2)
xlabel('SNR (dB)');ylabel('P_{fa}');
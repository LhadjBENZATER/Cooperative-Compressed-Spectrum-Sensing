function result = demo_ccss(makePlot)
% Synthetic coherent-channel cooperation example; no SDR files required.
if nargin<1, makePlot=true; end
root=fileparts(mfilename('fullpath')); addpath(root,fullfile(root,'examples'));
rng(20261121,'twister'); N=128; M=64; K=6; L=4; MC=20;
snrDb=[-15 -10 -5 0 5 10]; A=randn(M,N); A=A./sqrt(sum(A.^2,1));
truth=false(N,1); truth([12:14 80:82])=true;
x=zeros(N,1); x(truth)=1; clean=repmat(A*x,1,L); power=mean(abs(clean(:)).^2);
scores=zeros(numel(snrDb),2,3);
for s=1:numel(snrDb)
 for trial=1:MC
  Y=clean+sqrt(power/10^(snrDb(s)/10)/2)*(randn(M,L)+1i*randn(M,L));
  Yv=Y*(ones(L)-eye(L))/(L-1); local=false(N,L); transformed=local;
  for l=1:L
   local(:,l)=any(abs(somp_baseline(A,Y(:,l),K))>0,2);
   transformed(:,l)=any(abs(somp_baseline(A,Yv(:,l),K))>0,2);
  end
  decisions=[sum(local,2)>=ceil(L/2),all(transformed,2)];
  for method=1:2
   scores(s,method,:)=scores(s,method,:)+reshape(support_metrics(truth,decisions(:,method)),1,1,3)/MC;
  end
 end
end
result=struct('snrDb',snrDb,'metrics',scores,'seed',20261121,'trials',MC);
out=fullfile(root,'results'); if ~exist(out,'dir'),mkdir(out);end
rows=[];for m=1:2,rows=[rows;snrDb(:),repmat(m,numel(snrDb),1),squeeze(scores(:,m,:))];end
write_demo_csv(fullfile(out,'demo_metrics.csv'),'snr_db,method,pd,pf,f1',rows);
save(fullfile(out,'demo_results.mat'),'result'); fprintf('CCSS demo completed: %d paired trials/SNR.\n',MC);
if makePlot
 fig=figure('Visible','off','Color','w','Position',[100 100 1100 350]);
 for q=1:3
  subplot(1,3,q);plot(snrDb,scores(:,:,q),'-o','LineWidth',1.5);grid on;ylim([0 1]);
  xlabel('Measurement SNR (dB)');names={'P_d','P_f','F1'};ylabel(names{q});
 end
 subplot(1,3,1);legend('Physical channels: 2-of-4 voting','VDS channels: AND','Location','best');
 sgtitle('Synthetic coherent-channel example; no hardware measurements');
 print(fig,fullfile(root,'docs','preview.svg'),'-dsvg');close(fig);
end
end

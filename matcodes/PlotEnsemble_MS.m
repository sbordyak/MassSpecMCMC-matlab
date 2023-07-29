%%

figfolder = ['./results/' runname '/figs/'];
if ~exist(figfolder,'dir')
    mkdir(figfolder)
end


%load(datamat)
Ntb = d0.Ntb;







H=figure(1);clf
set(H,'color','w')
set(H,'units','inches','position',[6 2 6 11])
set(H,'invertHardcopy','off')
set(H,'units','inches','PaperSize', [6 11]);
set(H,'PaperPositionMode','auto')

fs = 10;


for ii = 1:4
    
    subplot(5,2,ii);
    set(gca,'FontSize',fs)
    [rhist,bins]=hist(ens_rat(ii,burn:cnt),50);
    dbins = bins(2)-bins(1);
    %bar(bins,rhist/(cnt-burn)/dbins,'Barwidth', 1.05,'EdgeColor','k','facecolor','c')
    aa=bar(bins,rhist,'Barwidth', 1.01,'EdgeColor','none','facecolor','b');hold on;
    xlabel(sprintf('Log %s/%s',Iso_Name{ii},Iso_Name{5}));ylabel('Frequency','FontSize',fs)
    %title(sprintf('Posterior for %s/%s',Iso_Name{m},Iso_Name{2}))
    spos=get(gca,'position');
    set(gca,'position',spos+[0 .02 0 0])
    %legend(aa,sprintf('%s/%s',Iso_Name{ii},Iso_Name{5}),'FontSize',fs)
    
end
% yy = get(gca,'ylim');
% xx=log(truedata.ratio208206^-1);
% plot(xx*[1 1],yy,'r--','LineWidth',1.5)

%%

subplot(5,2,5);
set(gca,'FontSize',fs)
mcol = cool(8);
for ii = 1:length(BLmean)
    [rhist,bins]=hist(ens_BL(ii,burn:cnt),50);
    dbins = bins(2)-bins(1);
    %bar(bins,rhist/(cnt-burn)/dbins,'Barwidth', 1.05,'EdgeColor','k','facecolor','c')
    bar(bins,rhist,'Barwidth', 1.01,'EdgeColor','none','facecolor',mcol(ii,:));
    hold on
end
xlabel('Baseline (counts)');ylabel('Frequency','FontSize',fs)
%title(sprintf('Posterior for %s/%s',Iso_Name{m},Iso_Name{2}))
%legend({'L1','H1'},'FontSize',fs)

% yy = get(gca,'ylim');
% yy=yy*1.6;
% set(gca,'ylim',yy)
% xx=massspec.H1baseline;
% plot(xx*[1 1],yy,'r--','LineWidth',1.5)
% xx=massspec.L1baseline;
% plot(xx*[1 1],yy,'r--','LineWidth',1.5)

for ii = 1:length(BLmean)
    [rhist,bins]=hist(ens_BL(ii,burn:cnt),50);
    dbins = bins(2)-bins(1);
    %bar(bins,rhist/(cnt-burn)/dbins,'Barwidth', 1.05,'EdgeColor','k','facecolor','c')
    bar(bins,rhist,'Barwidth', 1.01,'EdgeColor',mcol(ii,:),'facecolor',mcol(ii,:));
    hold on
end
%%

subplot(5,2,6);
set(gca,'FontSize',fs)
[rhist,bins]=hist(ens_DF(1,burn:cnt),50);
dbins = bins(2)-bins(1);
%bar(bins,rhist/(cnt-burn)/dbins,'Barwidth', 1.05,'EdgeColor','k','facecolor','c')
bar(bins,rhist,'Barwidth', 1.01,'EdgeColor','none','facecolor','b');hold on
xlabel('Daly/Faraday Gain');ylabel('Frequency','FontSize',fs)
%title(sprintf('Posterior for %s/%s',Iso_Name{m},Iso_Name{2}))
%legend(sprintf('D/F'))

% yy = get(gca,'ylim');
% xx=massspec.dalyFaradayGain;
% plot(xx*[1 1],yy,'r--','LineWidth',1.5)


%%
subplot(5,2,[7 8]);cla
set(gca,'FontSize',fs)

for m=1:d0.Nblock
 %m=1;   
    sampo = burn:1:cnt;
    [~,hi]=max(InterpMat{m});
    
    %plot(truedata.iPb208(301:end),'r','linewidth',1.5)    ;hold on
    plot(InterpMat{m}*Imean{m}/DFmean,'k--','linewidth',1);hold on
    plot(hi',Imean{m}/DFmean,'kv','markersize',7,'linewidth',1)
    
    xlabel('Time index','FontSize',fs);ylabel('Intensity (counts)','FontSize',fs)

    legend({'I_{mean}','knots'})
    
end

% %%
% subplot(5,2,[7 8])
% set(gca,'FontSize',fs)
% 
% %for m=1:d0.Nblock
%     
%     sampo = burn:1:cnt;
%     
%     %Iprint = (InterpMat{m}*ens_I{m}(:,sampo))/DFmean-repmat(Imean{m},1,length(sampo));
%     Iprint = (InterpMat{m}*ens_I{m}(:,sampo))/.9;
%     
%     
%     xedges=(1:Ntb(m))';yedges = linspace(min(Iprint(:)),max(Iprint(:)),100)';
%     for kk=1:Ntb(m)
%         Ihist(kk,:) = hist(Iprint(kk,:),yedges);
%     end
%     imagesc(xedges,yedges,Ihist')
%     hold on
%     plot(xedges,zeros(size(xedges)),'r--','linewidth',1.5)
%     
%     xlabel('Time index','FontSize',fs);ylabel('Residual intensity','FontSize',fs)
%     colormap([[1 1 1];paruly])
%     cb=colorbar('FontSize',fs);
%     
% %end


%%


% subplot(5,2,[9 10])
% set(gca,'FontSize',fs)
% 
% mcol = cool(8);
% for m=1:d0.Ndet-1
%     [rhist,bins]=hist(ens_sig(m,burn:cnt),30);
%     dbins = bins(2)-bins(1);
%     %bar(bins,rhist/(cnt-burn)/dbins,'facecolor',mcol(m,:),'edgecolor','none');hold on;
%     bar(bins,rhist./max(rhist),'facecolor',mcol(m,:),'edgecolor','none');hold on;
%     %bar(bins,rhist/(cnt-burn)/dbins,'Barwidth', 1.05,'EdgeColor','none');hold on;
%     xlabel('Noise hyperparmeter');ylabel('Frequency','FontSize',fs)
% end
% %legend({'L1','H1','Daly'})
% 
% spos=get(gca,'position');
% set(gca,'position',spos+[0 -.02 0 0])    
%     
% yy = get(gca,'ylim');
% yy=yy*1.6;
% set(gca,'ylim',yy)
% xx = get(gca,'xlim');
% xx(1)=xx(1) - 0.1*(xx(2)-xx(1));
% xx(2)=xx(2) + 0.2*(xx(2)-xx(1));
% set(gca,'xlim',xx)
% 
% % 
% % xx=sqrt(massspec.FnoiseVarianceCPS/integrationTime.reportinterval);
% % plot(xx*[1 1],yy,'r--','LineWidth',1.5)
% 
% %plot(0*[1 1],yy,'r--','LineWidth',1.5)

%%
set(H,'renderer','painters');
%print(H,sprintf('%s/Ensemble_run%02d.eps',figfolder,iset),'-depsc','-r200'); 




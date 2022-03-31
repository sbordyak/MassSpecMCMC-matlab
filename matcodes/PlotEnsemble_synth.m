%%

figfolder = ['./results/' runname '/figs/'];
if ~exist(figfolder,'dir')
    mkdir(figfolder)
end


load(datamat)
Ntb = d0.Ntb;



vir10 = viridis(10);



H=figure(1);clf
set(H,'color','w')
set(H,'units','inches','position',[12 2 6 11])
set(H,'invertHardcopy','off')
set(H,'units','inches','PaperSize', [6 11]);
set(H,'PaperPositionMode','auto')

fs = 10;




subplot(5,2,[1 2]);
set(gca,'FontSize',fs)
[rhist,bins]=hist(exp(ens_rat(1,burn:cnt)),50);
dbins = bins(2)-bins(1);
%bar(bins,rhist/(cnt-burn)/dbins,'Barwidth', 1.05,'EdgeColor','k','facecolor','c')
aa=bar(bins,rhist,'Barwidth', 1.01,'EdgeColor','none','facecolor',vir10(3,:));hold on;
xlabel('Ratio');ylabel('Frequency','FontSize',fs)
%title(sprintf('Posterior for %s/%s',Iso_Name{m},Iso_Name{2}))
spos=get(gca,'position');
set(gca,'position',spos+[0 .02 0 0]) 

yy = get(gca,'ylim');
xx=exp(log(truedata.ratio208206^-1));
plot(xx*[1 1],yy,'r--','LineWidth',1.5)

legend(aa,sprintf('%s/%s',Iso_Name{1},Iso_Name{2}));%,'FontSize',fs)

%%

subplot(5,2,3);
set(gca,'FontSize',fs)


vir10 = viridis(10);
mcol = vir10([4 8],:);
%mcol = viridis(2);
for ii = 1:length(BLmean)
    [rhist,bins]=hist(ens_BL(ii,burn:cnt)/6.24e7*1e6,50);
    dbins = bins(2)-bins(1);
    %bar(bins,rhist/(cnt-burn)/dbins,'Barwidth', 1.05,'EdgeColor','k','facecolor','c')
    bb(ii) = bar(bins,rhist,'Barwidth', 1.01,'EdgeColor','non','facecolor',mcol(ii,:));
    hold on
end
xlabel('Baseline (\muV)');ylabel('Frequency','FontSize',fs)
%title(sprintf('Posterior for %s/%s',Iso_Name{m},Iso_Name{2}))

yy = get(gca,'ylim');
yy=yy*1.6;
set(gca,'ylim',yy)
xx = get(gca,'xlim');
set(gca,'xlim',[xx(1)-(xx(2)-xx(1))/20, xx(2)+(xx(2)-xx(1))/20])

xx=massspec.H1baseline/6.24e7*1e6;
plot(xx*[1 1],yy,'r--','LineWidth',1.5)
hold on
xx=massspec.L1baseline/6.24e7*1e6;
plot(xx*[1 1],yy,'r--','LineWidth',1.5)



for ii = 1:length(BLmean)
    [rhist,bins]=hist(ens_BL(ii,burn:cnt)/6.24e7*1e6,50);
    dbins = bins(2)-bins(1);
    %bar(bins,rhist/(cnt-burn)/dbins,'Barwidth', 1.05,'EdgeColor','k','facecolor','c')
    bar(bins,rhist,'Barwidth', 1.01,'EdgeColor',mcol(ii,:),'facecolor',mcol(ii,:));
    hold on
end



legend(bb,{'L1','H1'})%,'FontSize',fs)

%%

subplot(5,2,4);
set(gca,'FontSize',fs)
[rhist,bins]=hist(ens_DF(1,burn:cnt),50);
dbins = bins(2)-bins(1);
%bar(bins,rhist/(cnt-burn)/dbins,'Barwidth', 1.05,'EdgeColor','k','facecolor','c')
bar(bins,rhist,'Barwidth', 1.01,'EdgeColor','none','facecolor',vir10(1,:));hold on
xlabel('Daly/Faraday Gain');ylabel('Frequency','FontSize',fs)
%title(sprintf('Posterior for %s/%s',Iso_Name{m},Iso_Name{2}))
%legend(sprintf('D/F'))

yy = get(gca,'ylim');
xx=massspec.dalyFaradayGain;
plot(xx*[1 1],yy,'r--','LineWidth',1.5)


%%
subplot(5,2,[5 6]);cla
set(gca,'FontSize',fs)

%for m=1:d0.Nblock
 m=1;   
    sampo = burn:1:cnt;
    [~,hi]=max(InterpMat{m});
    
    plot((0:1499)/5,truedata.iPb208(301:end)/6.24e7*1e6,'r','linewidth',1.5)    ;hold on
    plot((0:1499)/5,InterpMat{m}*Imean{m}/DFmean/6.24e7*1e6,'k--','linewidth',1)
    plot(hi'/5,Imean{m}/DFmean/6.24e7*1e6,'kv','markersize',7,'linewidth',1)
    
    xlabel('Time (s)','FontSize',fs);ylabel('Intensity (\muV)','FontSize',fs)

    legend({'I_{true}','I_{mean}','knots'})
    
%end

%%
subplot(5,2,[7 8])
set(gca,'FontSize',fs)

%for m=1:d0.Nblock
    
    sampo = burn:1:cnt;
    
    Iprint = ((InterpMat{m}*ens_I{m}(:,sampo))/DFmean-repmat(truedata.iPb208(301:end),1,length(sampo)))/6.24e7*1e6;
    %Iprint = (InterpMat{m}*ens_I{m}(:,sampo))/.9;
    
    
    xedges=(1:Ntb(m))';yedges = linspace(min(Iprint(:)),max(Iprint(:)),100)';
    for kk=1:Ntb(m)
        Ihist(kk,:) = hist(Iprint(kk,:),yedges);
    end
    imagesc(xedges/5,yedges,Ihist')
    hold on
    plot(xedges/5,zeros(size(xedges)),'r--','linewidth',1.5)
    
    xlabel('Time (s)','FontSize',fs);ylabel('Residual intensity (\muV)','FontSize',fs)
    %colormap([[1 1 1];paruly])
    colormap([[1 1 1];viridis(64)])
    cb=colorbar('FontSize',fs);
    
%end


%%


subplot(5,2,[9 10])
set(gca,'FontSize',fs)

mcol = vir10([4 8],:);

%mcol = cool(3);
for m=1:d0.Ndet-1
    [rhist,bins]=hist(ens_sig(m,burn:cnt),30);
    dbins = bins(2)-bins(1);
    %bar(bins,rhist/(cnt-burn)/dbins,'facecolor',mcol(m,:),'edgecolor','none');hold on;
    aa(m)=bar(bins,rhist,'facecolor',mcol(m,:),'edgecolor','none');hold on;
    %bar(bins,rhist/(cnt-burn)/dbins,'Barwidth', 1.05,'EdgeColor','none');hold on;
    xlabel('Noise hyperparmeter');ylabel('Frequency','FontSize',fs)
end

spos=get(gca,'position');
set(gca,'position',spos+[0 -.02 0 0])    
    
yy = get(gca,'ylim');
yy=yy*1.6;
set(gca,'ylim',yy)
xx = get(gca,'xlim');
xx(1)=xx(1) - 0.1*(xx(2)-xx(1));
xx(2)=xx(2) + 0.2*(xx(2)-xx(1));
set(gca,'xlim',xx)


xx=sqrt(massspec.FnoiseVarianceCPS/integrationTime.reportinterval);
plot(xx*[1 1],yy,'r--','LineWidth',1.5)

legend([aa(1) aa(2)],{'L1','H1'});%,'Daly'})

%plot(0*[1 1],yy,'r--','LineWidth',1.5)

%%
set(H,'renderer','painters');
print(H,sprintf('%s/Ensemble_run%02d.eps',figfolder,iset),'-depsc','-r200'); 




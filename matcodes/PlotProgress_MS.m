%%

figfolder = ['./results/' runname '/figs/'];
if ~exist(figfolder,'dir')
    mkdir(figfolder)
end


%load(datamat)
Ntb = d0.Ntb;


fs = 10;





H=figure(2);clf
set(H,'color','w')
set(H,'units','inches','position',[6 5 4 2.5*(d0.Niso-1)])
set(H,'invertHardcopy','off')
set(H,'units','inches','PaperSize', [4 2.5*(d0.Niso-1)]);
set(H,'PaperPositionMode','auto')

for ii=1:d0.Niso-1
    subplot(d0.Niso-1,1,ii)
    aa=plot(ens_rat(ii,1:cnt),'k');set(gca,'xscale','log')
    xlabel('Saved iteration','FontSize',fs)
    ylabel(sprintf('log %s/%s',Iso_Name{ii},Iso_Name{end}),'FontSize',fs)
end


set(H,'renderer','painters');
print(H,sprintf('%s/Converge_Ratio_run%02d.eps',figfolder,iset),'-depsc','-r200'); 

%%


H=figure(3);clf
set(H,'color','w')
set(H,'units','inches','position',[6 5 4 5])
set(H,'invertHardcopy','off')
set(H,'units','inches','PaperSize', [4 5]);
set(H,'PaperPositionMode','auto')


subplot(2,1,1)
aa=plot(sqrt(ens_E(1:cnt)),'k');set(gca,'xscale','log','yscale','log')
xlabel('Saved iteration','FontSize',fs)
ylabel(sprintf('Weighted misfit'),'FontSize',fs)

subplot(2,1,2)
aa=plot(sqrt(ens_E0(1:cnt)),'k');set(gca,'xscale','log','yscale','log')
xlabel('Saved iteration','FontSize',fs)
ylabel(sprintf('Raw misfit'),'FontSize',fs)


set(H,'renderer','painters');
print(H,sprintf('%s/Converge_Error_run%02d.eps',figfolder,iset),'-depsc','-r200'); 


%%


H=figure(4);clf
set(H,'color','w')
set(H,'units','inches','position',[6 5 4 2.5*(d0.Nfar)])
set(H,'invertHardcopy','off')
set(H,'units','inches','PaperSize', [4 2.5*(d0.Nfar)]);
set(H,'PaperPositionMode','auto')

Farname = {'L4','L3','L2','L1','H1','B2','B3','B4'};

for ii=1:d0.Nfar
    subplot(d0.Nfar,1,ii)
    aa=plot(ens_BL(ii,1:cnt),'k');set(gca,'xscale','log')
    xlabel('Saved iteration','FontSize',fs)
    ylabel(sprintf('%s baseline',Farname{ii}),'FontSize',fs)
end


set(H,'renderer','painters');
print(H,sprintf('%s/Converge_Baseline_run%02d.eps',figfolder,iset),'-depsc','-r200'); 


%% Noise Hyperparameter plot (not used currently)


% H=figure(5);clf
% set(H,'color','w')
% set(H,'units','inches','position',[18 5 4 2.5*(d0.Ndet)])
% set(H,'invertHardcopy','off')
% set(H,'units','inches','PaperSize', [4 2.5*(d0.Ndet)]);
% set(H,'PaperPositionMode','auto')
% 
% hypername = {'L4','L3','L2','L1','H1','B2','B3','B4','Daly'};
% 
% for ii=1:d0.Ndet
%     subplot(d0.Ndet,1,ii)
%     aa=plot(ens_sig(ii,1:cnt),'k');set(gca,'xscale','log')
%     xlabel('Saved iteration','FontSize',fs)
%     ylabel(sprintf('%s baseline',hypername{ii}),'FontSize',fs)
% end
% 
% 
% set(H,'renderer','painters');
% print(H,sprintf('%s/Converge_Noise_run%02d.eps',figfolder,iset),'-depsc','-r200'); 
% 



%%


H=figure(6);clf
set(H,'color','w')
set(H,'units','inches','position',[6 5 4 2.5*(d0.Nblock)])
set(H,'invertHardcopy','off')
set(H,'units','inches','PaperSize', [4 2.5*(d0.Nblock)]);
set(H,'PaperPositionMode','auto')


for ii=1:d0.Nblock
    subplot(d0.Nblock,1,ii)
    Icolor = cool(d0.Ncycle(ii));

    for jj=1:d0.Ncycle(ii)
    aa=plot(ens_I{ii}(jj,1:cnt),'color',Icolor(jj,:));hold on;set(gca,'xscale','log')
    end
    xlabel('Saved iteration','FontSize',fs)
    ylabel(sprintf('Block %d',ii),'FontSize',fs)
    set(gca,'xlim',[0 cnt])
    
end


set(H,'renderer','painters');
print(H,sprintf('%s/Converge_Intensity_run%02d.eps',figfolder,iset),'-depsc','-r200'); 





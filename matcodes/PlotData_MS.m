
%%

figfolder = ['./results/' runname '/figs/'];
if ~exist(figfolder,'dir')
    mkdir(figfolder)
end


% load(datamat)
Ntb = d0.Ntb;


fs = 10;



% for n = 1:d0.Nblock
%     
%     Intensity{n} = InterpMat{n}*x.I{n};
%     
%     for mm=1:d0.Niso;
%         itmp = d0.iso_ind(:,mm) & d0.axflag & d0.block(:,n);
%         %d(itmp,1) = (x.lograt(mm))*Intensity{n}(d0.time_ind(itmp)); %debug
%         %dnobl(itmp,1) = (x.lograt(mm))*Intensity{n}(d0.time_ind(itmp)); %debug
%         d(itmp,1) = exp(x.lograt(mm))*Intensity{n}(d0.time_ind(itmp));
%         dnobl(itmp,1) = exp(x.lograt(mm))*Intensity{n}(d0.time_ind(itmp));
%         
%         itmp = d0.iso_ind(:,mm) & ~d0.axflag & d0.block(:,n);
%         %d(itmp,1) = (x.lograt(mm))*x.DFgain^-1 *Intensity{n}(d0.time_ind(itmp)) + x.BL(d0.det_vec(itmp)); %debug
%         %dnobl(itmp,1) = (x.lograt(mm))*x.DFgain^-1 *Intensity{n}(d0.time_ind(itmp)); %debug        
%         d(itmp,1) = exp(x.lograt(mm))*x.DFgain^-1 *Intensity{n}(d0.time_ind(itmp)) + x.BL(d0.det_vec(itmp)); 
%         dnobl(itmp,1) = exp(x.lograt(mm))*x.DFgain^-1 *Intensity{n}(d0.time_ind(itmp));
%     end
% end

    d = ModelMSData(x0,d0);


%Dsig = sqrt(x.sig(d0.det_vec).^2 + x.sig(d0.iso_vec+d0.Ndet).*dnobl); % New data covar vector
Dsig = sqrt(x0.Dsig); %sqrt(x.sig(d0.det_vec).^2 + x.sig(end).*dnobl); % New data covar vector   


rnr = 10;
reduce = (1:rnr:d0.Ndata)';



restmp = (d0.data-d).^2;
        

E=sum(restmp./Dsig.^2);
E0=sum(restmp);  % Unweighted error func (for visualization)


xx=1:length(d);


H=figure(7);clf
set(H,'color','w')
set(H,'units','inches','position',[6 5 7.4 6.4])
set(H,'invertHardcopy','off')
set(H,'units','inches','PaperSize', [7.5 6.5]);
set(H,'PaperPositionMode','auto')


isotype = d0.iso_vec + 1;
isocol = [0 0 0; flipud(jet(d0.Niso))];

meastype = 1+d0.blflag*2 + d0.axflag;
measmark = {'v','+','o'};


subplot(2,1,1)
spos = get(gca,'position');set(gca,'position',spos+[0 .05 0 -.05])
%plot(d0.data,'k.');hold on;
scatter(xx(reduce),d0.data(reduce),10,isotype(reduce));hold on
plot(d,'r');
set(gca,'xlim',[1,length(d)])
% yl=get(gca,'ylim');set(gca,'ylim',[yl(1) yl(2)*1.2])
xlabel('Data index','FontSize',fs)
ylabel(sprintf('Data (counts)'),'FontSize',fs)
box on
colormap(isocol)

%%

reduce = zeros(d0.Ndata,1); reduce(1:rnr:end)=1;

subplot(2,1,2)
spos = get(gca,'position');set(gca,'position',spos+[0 .15 0 -.05])
for ii = 1:3
    ind = meastype==ii & reduce;
    sc(ii)=scatter(xx(ind),d0.data(ind)-d(ind),10,isotype(ind,:),'marker',measmark{ii});hold on;
    colormap(isocol)
end

plot(zeros(size(d)),'r--','linewidth',1.2);
plot(1*(Dsig),'r','linewidth',1.2);plot(-1*(Dsig),'r','linewidth',1.2)

set(gca,'xlim',[1,length(d)])
% yl=get(gca,'ylim');set(gca,'ylim',[yl(1)*1 yl(2)*2])


xlabel('Data index','FontSize',fs)
ylabel(sprintf('Residual counts'),'FontSize',fs)
box on

%%
axes
spos = get(gca,'position');set(gca,'position',spos+[0 -.1 0 0])

sy(1) = plot(1,1,'kv','linewidth',2);hold on
sy(2) = plot(1,1,'k+','linewidth',2);
sy(3) = plot(1,1,'ko','linewidth',2);

legend(sy,{'Faraday','Daly','Baseline'},'Location','southeast','fontsize',fs)

set(gca,'visible','off')
set(sy,'visible','off')

%%
axes
spos = get(gca,'position');set(gca,'position',spos+[0 -.1 0 0])

for ii = 1:d0.Niso+1;
    c(ii) = plot(1,1,'o','color',isocol(ii,:),'linewidth',2);hold on
end
cleg{1} = 'Baseline';

for ii=1:d0.Niso;
    cleg{ii+1} = sprintf('%s',Iso_Name{ii});
end
    
legend(c,cleg,'Location','southwest','fontsize',fs)


set(c,'visible','off')
set(gca,'visible','off')

%%
set(H,'renderer','painters');
print(H,sprintf('%s/DataFit__run%02d.pdf',figfolder,iset),'-dpdf','-r200'); 





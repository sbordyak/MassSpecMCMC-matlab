function [x0,C0,d] = InitializeModelCov_synth(d0,user_DFgain)



%% CREATE INITIAL MODEL

for m=1:d0.Nfar%+1
    x0.BL(m,1) = mean(d0.data(d0.blflag & d0.det_ind(:,m)));
    x0.BLstd(m,1) = std(d0.data(d0.blflag & d0.det_ind(:,m)));
end
%x0.BL(d0.Nfar+1,1) = 0;


iden = d0.Niso; % Denominator Isotope

x0.DFgain = user_DFgain;

for m=1:d0.Nblock
    II = d0.InterpMat{m};
    
    dind = (d0.iso_ind(:,iden) &   d0.block(:,m));
    dd=d0.data(dind);
    tmpdetvec = d0.det_vec(dind);
    
    dd(~d0.axflag(dind)) = (dd(~d0.axflag(dind)) - x0.BL(tmpdetvec(~d0.axflag(dind))))*x0.DFgain ;
    
    tmptime = d0.time_ind(dind);
    [~,dsort]=sort(tmptime);
    
    dd=dd(dsort);
    
    %I0=(II'*II)^-1*II'*dd;
    I0=II\dd;  % Solve least squares problem for initial intensity values
    
    %x0.I{m} =  tmpDenIso*ones(Nknots(m),1);
    %x0.I{m} = linspace(maxtmpCounts(m,iden),mintmpCounts(m,iden),d0.Nknots(m))';
    x0.I{m} = I0;
end





for m=1:d0.Nblock
    II = d0.InterpMat{m};
    IntensityFn = II*x0.I{m};
    
    for ii = 1:d0.Niso
        dind = (d0.iso_ind(:,ii) &   d0.block(:,m));
        dd=d0.data(dind);
        tmpdetvec = d0.det_vec(dind);
        
        dd(~d0.axflag(dind)) = (dd(~d0.axflag(dind)) - x0.BL(tmpdetvec(~d0.axflag(dind))))*x0.DFgain ;
        
        tmptime = d0.time_ind(dind);
        [~,dsort]=sort(tmptime);
        
        dd=dd(dsort);
        
        IsoBlockMean(ii,m) = log(mean(dd./IntensityFn));

    end
end

for ii=1:d0.Niso
    x0.lograt(ii,1) = mean(IsoBlockMean(ii,:));
end

%%
Dsig = zeros(size(d0.data));

ReportInterval = 0.1;


II = d0.InterpMat;
for n = 1:d0.Nblock
    Intensity{n} = II{n}*x0.I{n};
    for m=1:d0.Niso;
        itmp = d0.iso_ind(:,m) & d0.axflag & d0.block(:,n);
        ddd(itmp,1) = exp(x0.lograt(m))*Intensity{n}(d0.time_ind(itmp));
        Dsig(itmp,1) = ddd(itmp)/ReportInterval;
        
        itmp = d0.iso_ind(:,m) & ~d0.axflag & d0.block(:,n);
        ddd(itmp) = exp(x0.lograt(m))*x0.DFgain^-1 *Intensity{n}(d0.time_ind(itmp));
    end
    for m=1:d0.Nfar%+1
        itmp = d0.det_ind(:,m) & ~d0.axflag & d0.block(:,n);
        Dsig(itmp,1) = ddd(itmp)/ReportInterval + x0.BLstd(m).^2;
    end
    
end
for m=1:d0.Nfar%+1
    Dsig(d0.blflag & d0.det_ind(:,m),1) = x0.BLstd(m).^2;
end

x0.Dsig = Dsig;


%%

for m = 1:d0.Niso
    testLR = linspace(-.5,.5,1001);
    Etmp = zeros(size(testLR));
    for ii = 1:length(testLR) 
        testx0 = x0;
        testx0.lograt(m) = x0.lograt(m) + testLR(ii);
        d = ModelMSData(testx0,d0);
        Etmp(ii) = sum((d0.data-d).^2./Dsig);
    end    
    EE=(Etmp-min(Etmp));
    p = exp(-EE/2)/sum(exp(-EE/2)); 
    x0.logratVar(m) = sum(p.*(testLR-0).^2);
    
end

for n = 1:d0.Nblock
    for m = 1:length(x0.I{n})
        testI = linspace(-mean(x0.BLstd),mean(x0.BLstd),101);
        Etmp = zeros(size(testI));
        for ii = 1:length(testI)
            testx0 = x0;
            testx0.I{n}(m) = x0.I{n}(m) + testI(ii);
            d = ModelMSData(testx0,d0);
            Etmp(ii) = sum((d0.data-d).^2./Dsig);
        end
        EE=(Etmp-min(Etmp));
        p = exp(-EE/2)/sum(exp(-EE/2));
        x0.IVar{n}(m) = sum(p.*(testI-0).^2);
        
    end
end

%%

testDF = linspace(-.1,.1,1001);
Etmp = zeros(size(testDF));
for ii = 1:length(testDF)
    testx0 = x0;
    testx0.DFgain = x0.DFgain + testDF(ii);
    d = ModelMSData(testx0,d0);
    Etmp(ii) = sum((d0.data-d).^2./Dsig);
end
EE=(Etmp-min(Etmp));
p = exp(-EE/2)/sum(exp(-EE/2)); 
x0.DFgainVar = sum(p.*(testDF-0).^2);



for m = 1:d0.Nfar	
    testBL = linspace(-x0.BLstd(m),x0.BLstd(m),1001);
    Etmp = zeros(size(testBL));
    for ii = 1:length(testBL) 
        testx0 = x0;
        testx0.BL(m) = x0.BL(m) + testBL(ii);
        d = ModelMSData(testx0,d0);
        Etmp(ii) = sum((d0.data-d).^2./Dsig);
    end  
    EE=(Etmp-min(Etmp));
    p = exp(-EE/2)/sum(exp(-EE/2)); 
    x0.BLVar(m) = sum(p.*(testBL-0).^2);
    
end



d = ModelMSData(x0,d0);


C0diag =  [sqrt(x0.logratVar) sqrt([x0.IVar{:}]) sqrt(x0.BLVar) sqrt(x0.DFgainVar)];
%C0diag =  [(x0.logratVar) ([x0.IVar{:}]) (x0.BLVar) (x0.DFgainVar)];
d0.Nmod = length(C0diag);
C0 = (0.1)^2*d0.Nmod^-1*diag(C0diag);


%return

% %% MODEL DATA WITH INITIAL MODEL
% II = d0.InterpMat;
% 
% for m=1:d0.Nfar%+1
%     d(d0.blflag & d0.det_ind(:,m),1) = x0.BL(m);
% end
% 
% for n = 1:d0.Nblock
%     Intensity{n} = II{n}*x0.I{n};
%     for m=1:d0.Niso;
%         itmp = d0.iso_ind(:,m) & d0.axflag & d0.block(:,n);
%         d(itmp) = exp(x0.lograt(m))*Intensity{n}(d0.time_ind(itmp));
%         
%         itmp = d0.iso_ind(:,m) & ~d0.axflag & d0.block(:,n);
%         d(itmp) = exp(x0.lograt(m))*x0.DFgain^-1 *Intensity{n}(d0.time_ind(itmp)) + x0.BL(d0.det_vec(itmp));
%     end
% end
% 
% for m = 1:d0.Nfar%+1
%     itmp = d0.det_vec==m & d0.blflag==1;
%     x0.sig(m,1) = 1*std(d0.data(itmp));
% end
% 
% x0.sig(d0.Nfar+1,1) = 0;
% 
% for m = 1: d0.Niso;
%     itmp = d0.iso_vec==m ;
%     x0.sig(d0.Ndet + m,1) = 1.1*10;
% end
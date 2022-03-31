function [x0,d,Intensity] = InitializeModel_synth(d0)



%% CREATE INITIAL MODEL

for m=1:d0.Nfar%+1
    x0.BL(m,1) = mean(d0.data(d0.blflag & d0.det_ind(:,m)));
    x0.BLstd(m,1) = std(d0.data(d0.blflag & d0.det_ind(:,m)));
end
%x0.BL(d0.Nfar+1,1) = 0;


for m=1:d0.Niso;
    tmpCounts(m,1) = mean(d0.data( (d0.iso_ind(:,m) & d0.axflag)));
    
    for n = 1:d0.Nblock
        maxtmpCounts(n,m) = max(d0.data( (d0.iso_ind(:,m) & d0.axflag & d0.block(:,n))));
        mintmpCounts(n,m) = min(d0.data( (d0.iso_ind(:,m) & d0.axflag & d0.block(:,n))));
    end
    
    itmp = (d0.iso_ind(:,m) & ~d0.axflag);
    tmpFar(m,1)  = mean(d0.data(itmp)-x0.BL(d0.det_vec(itmp)));
end
[~,imaxC] = max(tmpCounts);

iden = d0.Niso;

x0.DFgain = tmpCounts(imaxC)/tmpFar(imaxC);


for m=1:d0.Niso
    x0.lograt(m,1) = log(tmpCounts(m)/tmpCounts(iden));
end



%%




for m=1:d0.Nblock
    II = d0.InterpMat{m};
    
    dind = ( d0.axflag & d0.block(:,m));
    dd=d0.data(dind)./exp(x0.lograt(d0.iso_vec(dind)));
    [~,dsort]=sort(d0.time_ind(dind));
    
    dd=dd(dsort);
    
    I0=(II'*II)^-1*II'*dd;
    
    %x0.I{m} =  tmpDenIso*ones(Nknots(m),1);
    %x0.I{m} = linspace(maxtmpCounts(m,iden),mintmpCounts(m,iden),d0.Nknots(m))';
    x0.I{m} = I0;
end





%%% MODEL DATA WITH INITIAL MODEL
II = d0.InterpMat;

for m=1:d0.Nfar%+1
    d(d0.blflag & d0.det_ind(:,m),1) = x0.BL(m);
end

for n = 1:d0.Nblock
    Intensity{n} = II{n}*x0.I{n};
    for m=1:d0.Niso;
        itmp = d0.iso_ind(:,m) & d0.axflag & d0.block(:,n);
        d(itmp) = exp(x0.lograt(m))*Intensity{n}(d0.time_ind(itmp));
        
        itmp = d0.iso_ind(:,m) & ~d0.axflag & d0.block(:,n);
        d(itmp) = exp(x0.lograt(m))*x0.DFgain^-1 *Intensity{n}(d0.time_ind(itmp)) + x0.BL(d0.det_vec(itmp));
    end
end



% % Define initial sigmas based on current error
% x0.sig = 2*std(d-d0.data)*ones(d0.Nsig,1);
% for m = 1:d0.Nsig
%     itmp = d0.sig_ind==m;
%     x0.sig(m) = 1*std(d(itmp)-d0.data(itmp));
% end

% Define initial sigmas based on baseline
%x0.sig = 2*std(d-d0.data)*ones(d0.Nsig,1);
for m = 1:d0.Nfar%+1
    itmp = d0.det_vec==m & d0.blflag==1;
    x0.sig(m,1) = 1*std(d0.data(itmp));
end

x0.sig(d0.Nfar+1,1) = 0;

for m = 1: d0.Niso;
    itmp = d0.iso_vec==m ;
    x0.sig(d0.Ndet + m,1) = 1.1*10;
end
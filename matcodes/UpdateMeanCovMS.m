function [xmean,xcov] = UpdateMeanCovMS(x,xcov,xmean,ensemble,m,iterflag)



Niso = length(x.lograt);
Nblock = length(x.I);
for ii=1:Nblock;
    Ncycle(ii) = length(x.I{ii});
end
Nfar = length(x.BL);
Ndf = 1;

Nmod = Niso + sum(Ncycle) + Nfar + Ndf;


if iterflag
    
    
    xx = x.lograt;
    
    for ii=1:Nblock
        xx = [xx; x.I{ii}];
    end
    
    xx = [xx; x.BL(1:Nfar)];
    
    xx = [xx; x.DFgain];
    
    
    
    xmean = (xmean*(m-1) + xx)/m;
    
    
    xctmp = (xx-xmean)*(xx-xmean)';
    xctmp = (xctmp+xctmp')/2;
    
    
    xcov = (xcov*(m-1) + (m-1)/m*xctmp)/m;
    
    
end




if ~iterflag
    
    cnt = length(ensemble);
    
    enso = [ensemble.lograt];
    
    for ii = 1:Nblock
        for n = 1:cnt;
            ens_I{ii}(:,n) =[ensemble(n).I{ii}];
        end
        enso = [enso; ens_I{ii}];
        
    end
    enso = [enso; [ensemble.BL]];
    enso = [enso; [ensemble.DFgain]];
    
    %xcov = cov(enso(:,ceil(end/2):end)');
    xmean = mean(enso(:,m:end)');
    xcov = cov(enso(:,m:end)');
    
    
    
end

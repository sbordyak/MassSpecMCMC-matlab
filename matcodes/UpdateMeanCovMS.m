function [xmean,xcov] = UpdateMeanCovMS(x,xmean,xcov,m)


Niso = length(x.lograt);
Nblock = length(x.I);
for ii=1:Nblock;
    Ncycle(ii) = length(x.I{ii});
end
Nfar = length(x.BL);
Ndf = 1;

Nmod = Niso + sum(Ncycle) + Nfar + Ndf;


xx = x.lograt;
for ii=1:Nblock
    xx = [xx; x.I{ii}];
end
xx = [xx; x.BL(1:Nfar)];
xx = [xx; x.DFgain];

xmeantmp = xmean;
xmean = xmeantmp + (xx-xmeantmp)/m;

%xmean = (xmean*(m-1) + xx)/m;

%xctmp = (xx-xmean)*(xx-xmean)';
%xctmp = (xctmp+xctmp')/2;


%xcov = (xcov*(m-1) + (m-1)/m*xctmp)/m;

xcov = xcov*(m-1)/m + (m-1)/m^2*(xx-xmean)*(xx-xmeantmp)';

    
    


% 
%         xmeantmp = xmean;
%         xmean = xmeantmp + (x-xmeantmp)/(cnt);
%         xcov = xcov*(cnt)/(cnt+1) + (cnt)/(cnt+1)^2*(x-xmean)'*(x-xmeantmp);
%         mxa(ii,:)=xmean;
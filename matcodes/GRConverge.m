    
function Rexit = GRConverge(x,ensemble);

Niso = length(x.lograt)-1;
Nblock = length(x.I);
for ii=1:Nblock;
    Ncycle(ii) = length(x.I{ii});
end
Nfar = length(x.BL);
Ndf = 1;

Nmod = Niso + sum(Ncycle) + Nfar + Ndf;


cnt = length(ensemble);

xall = [ensemble.lograt];
xall = xall(1:Niso,:);

for ii = 1:Nblock
    for n = 1:cnt;
        ens_I{ii}(:,n) =[ensemble(n).I{ii}];
    end
    xall = [xall; ens_I{ii}];
    
end
xall = [xall; [ensemble.BL]];
xall = [xall; [ensemble.DFgain]];
xall = xall';


    ngroup = round(sqrt(cnt));
    gsize = round(sqrt(cnt));
    for jj = 1:ngroup % Iterate over ngroups groups of size gsize
        tmpxs(:,:,jj) = cov(xall(1+(jj-1)*gsize:jj*gsize,:));
        tmpxm(:,jj) = mean(xall(1+(jj-1)*gsize:jj*gsize,:));
    end
    MeanofVar = sum(tmpxs(:,:,1:ngroup),3)/ngroup; % Mean of variances
    VarofMean = diag(std(tmpxm(:,1:ngroup),[],2).^2); % Variance of means
    
    %Rexit(ngroup) = sqrt((ngroup-1)/ngroup+(det(VarofMean)/det(MeanofVar))^(1/N)/ngroup);
    Rexit = sqrt((ngroup-1)/ngroup+(det(VarofMean)/det(MeanofVar))^(1/Nmod)/ngroup);
    
%     if Rexit(ngroup)<=ExitCrit
%         disp(sprintf('MCMC exiting after %d iters with R of %0.6f',cnt,Rexit(ngroup)))
%         break
%     end
    

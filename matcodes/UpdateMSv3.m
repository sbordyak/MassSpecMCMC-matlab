function  [x2,delx,xcov] = UpdateMSv3(oper,x,psig,prior,ensemble,xcov,delx_adapt,adaptflag,allflag)

%%
cnt = length(ensemble);

Niso = length(x.lograt);
Nblock = length(x.I);
for ii=1:Nblock;
    Ncycle(ii) = length(x.I{ii});
end
Nfar = length(x.BL);
Ndf = 1;



%xind = [ones(Niso,1); psig.I*ones(sum(Ncycle),1); psig.BL*ones(Nfar,1); psig.DFgain*ones(Ndf,1)];

ps0diag =  [psig.lograt*ones(Niso,1);          psig.I*ones(sum(Ncycle),1);     psig.BL*ones(Nfar,1);     psig.DFgain*ones(Ndf,1)];
priormin = [prior.lograt(1)*ones(Niso-1,1); 0; prior.I(1)*ones(sum(Ncycle),1); prior.BL(1)*ones(Nfar,1); prior.DFgain(1)*ones(Ndf,1)];
priormax = [prior.lograt(2)*ones(Niso-1,1); 0; prior.I(2)*ones(sum(Ncycle),1); prior.BL(2)*ones(Nfar,1); prior.DFgain(2)*ones(Ndf,1)];


%xx0 = [x.lograt; x.I{1}; x.I{2}; x.BL; x.DFgain];


xx0 = x.lograt;
xind = ones(Niso,1);

for ii=1:Nblock
    xx0 = [xx0; x.I{ii}];
    xind = [xind; 1+ii*ones(Ncycle(ii),1)];
end

xx0 = [xx0; x.BL];
xind = [xind; (2+Nblock)*ones(Nfar,1)];

xx0 = [xx0; x.DFgain];
xind = [xind; (3+Nblock)*ones(Ndf,1)];



if strcmp(oper(1:3),'cha')
    
    
    
    if ~allflag
        

        if strcmp(oper,'changer')
            nind = randi(Niso-1);
        elseif strcmp(oper,'changeI')
            nind = Niso+randi(sum(Ncycle));
        elseif strcmp(oper,'changebl')
            nind = Niso + sum(Ncycle) + randi(Nfar);
        elseif strcmp(oper,'changedfg')
            nind = Niso + sum(Ncycle) + Nfar + randi(Ndf);
        end
        
        
        if adaptflag
            delx = sqrt(xcov(nind,nind))*randn(1); %mvnrnd(zeros(1),xcov(nind,nind));
        else
            delx = ps0diag(nind)*randn(1);
        end
        
        
        xx =  xx0;
        xx(nind) = xx(nind) + delx;
        
        inprior = xx<=priormax & xx>=priormin;
        
        xx(~inprior) = xx0(~inprior);
        
        
    else
        
        %
        %delx(:,1) = mvnrnd(zeros(size(xx0)),2.38^2*diag(ps0diag).^2)/10/length(xind);
        %
        %     %%
        %         enso = [ensemble.lograt];
        %
        %         for ii = 1:Nblock
        %             for n = 1:cnt;
        %                 ens_I{ii}(:,n) =[ensemble(n).I{ii}];
        %             end
        %             enso = [enso; ens_I{ii}];
        %
        %         end
        %         enso = [enso; [ensemble.BL]];
        %         enso = [enso; [ensemble.DFgain]];
        %
        %         xcov = cov(enso(:,ceil(end/2):end)');
        %         %%
        %         delx(:,1) = mvnrnd(zeros(size(xx0)),2.38^2*xcov/length(xind));
        
        
        
        
        %VARY ALL AT A TIME
        
        %%delx(:,1) = mvnrnd(zeros(size(xx0)),xcov);
        %delx(:,1) = mvnrnd(zeros(size(xx0)),2.38^2*xcov/length(xind));
        delx = delx_adapt;
        
        xx =  xx0 + delx;
        
        inprior = xx<=priormax & xx>=priormin;
        
        xx(~inprior) = xx0(~inprior);
        
        
    end
    
    
    
    x2.lograt = xx(xind==1);
    for ii=1:Nblock
        x2.I{ii} = xx(xind==(1+ii));
    end
    x2.BL = xx(xind==(2+Nblock));
    x2.DFgain = xx(xind==(3+Nblock));
    
    x2.sig = x.sig;
    
    
    
    
    
    
elseif strcmp(oper(1:3),'noi')    %CHANGE NOISE
    
    %nind = randi(length(x.sig));
    nind = randi(length(x.BL)); %Just for the faradays
    
    
    x2=x;
    
    
    delx=psig.sig*randn(1);
    
    if x2.sig(nind) + delx >= prior.sig(1) && x2.sig(nind) + delx <= prior.sig(2)
        x2.sig(nind) = x2.sig(nind)+delx;
    else
        delx=0;
    end
    
    
else
    display('Thats not a thing')
    
    
    
end

%

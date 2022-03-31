function keep = AcceptItMS(oper,dE,psig,delx,prior,Dsig,Dsig2,d0)

% Determine probability of acceptance based on updated error function,
% Number of model parameters (may add later)

if strcmp(oper(1:3),'cha')
    
    X = 1;

    keep = min(1,X*exp(-(dE)/2));
    
    
    
elseif  strcmp(oper,'noise')
    
    %X = sum(-log(Dsig2(d0.blflag)))-sum(-log(Dsig(d0.blflag)));
    X = sum(-log(Dsig2))-sum(-log(Dsig));
    keep = min(1,exp(X/2-(dE)/2));

    %display([X dE/2 keep])
    
else
    display('Nope, come on now')
end
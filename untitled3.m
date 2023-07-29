
%%


for m = d0.Niso-1%:d0.Niso-1
    testLR = linspace(-.001,.001,101);
    %sb629 
    delta_testLR = testLR(2)-testLR(1);
    minvarLR = (delta_testLR/2)^2;
    Etmp = zeros(size(testLR));
    for ii = 1:length(testLR) 
        testx0 = x0;
        testx0.lograt(m) = x0.lograt(m) + testLR(ii);
        d = ModelMSData(testx0,d0);
        Etmp(ii) = sum((d0.data-d).^2./Dsig);
    end    
    EE=(Etmp-min(Etmp));
    p = exp(-EE/2)/sum(exp(-EE/2)); 
    %x0.logratVar(m) = max(sum(p.*(testLR-0).^2),minvarLR);%sb629 
    
end
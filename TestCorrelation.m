% Test correlation between model parameters that might be causing difficulty 
% In MCMC algorithm.

p1 = 1000; % Amount of Isotope 1
p2 = 1; % Amount of Isotope 2

% Noise on isotope samples as proportion
s1 = p1/10;  
s2 = p2/100;

% Generate N random samples of each
N=100;
d1 = p1 + s1*randn(N,1);
d2 = p2 + s2*randn(N,1);

% MCMC iterations
M = 100000;

% Proposal distributions
proplograt = 0.002;
prop1 = p1/400;
propnoise = 0.001;

% Allocate parameters vectors
xlograt = zeros(M,1); % Log ratio log(p2/p1)
x1 = zeros(M,1); % Intensity of Isotope 1
xnoise =  zeros(M,2); 

% Set equal to correct values to start
x1(1) = p1;% + randn*prop1;
xlograt(1) = log(p2/p1);% + randn*proplograt;
xnoise(1,:) = log([s1 s2]);

kept=0;

% Begin MCMC
for ii = 1:M-1
    nind = randi(2); % Choose param to update
    
    if nind==1 % If intensity of 1
        x1tmp = x1(ii) + randn*prop1; % Update this
        xlrtmp = xlograt(ii);
        xntmp = xnoise(ii,:);
    elseif nind==2  % If log ratio
        x1tmp = x1(ii);
        xlrtmp = xlograt(ii) + randn*proplograt; % Update that
        xntmp = xnoise(ii,:);
    else
        x1tmp = x1(ii);
        xlrtmp = xlograt(ii);
           
        nind2 = randi(2);
        xntmp = xnoise(ii,:);
        xntmp(1,nind2) = xntmp(1,nind2) + randn*propnoise;
    end

    % Model isotope amounts for old values
    mod01 = x1(ii);
    mod02 = x1(ii)*exp(xlograt(ii));
    mod0s = exp(xnoise(ii,:));
    
    % Model isotope amounts for new values
    mod1 = x1tmp;
    mod2 = x1tmp*exp(xlrtmp);
    mods = exp(xntmp);
    
    % Calculate misfit with both sets
    E =  sum((d1-mod1).^2/2/mods(1).^2+(d2-mod2).^2/2/mods(2).^2);
    E0 = sum((d1-mod01).^2/2/mod0s(1).^2+(d2-mod02).^2/2/mod0s(2).^2);
    
    % Determine acceptance probability
    dE = E-E0;
    keep = exp(-dE/2);
    
    if keep>rand(1) % If keep new values
        x1(ii+1) = x1tmp;
        xlograt(ii+1) = xlrtmp;
        xnoise(ii+1,:) = xntmp;

        kept=kept+1;
    else % If keep old values
        x1(ii+1) = x1(ii);
        xlograt(ii+1) = xlograt(ii);
        xnoise(ii+1,:) = xnoise(ii,:);
    end
    
end

burn = M/2; % Burn in period

% Trim vectors
x1 = x1(burn+1:end);
xlograt = xlograt(burn+1:end);

%%

figure
set(gcf,'Position',[360   480   522   218])
subplot(1,2,1)
plot(x1,xlograt,'k.',p1,log(p2/p1),'rv',mean(d1),log(mean(d2)/mean(d1)),'co')
xlabel('Isotope 1');ylabel('Log Ratio')
title(sprintf('S1=%.0e, S2=%.0e',s1/p1,s2/p2));
legend(sprintf('Corr = %.2f',corr(x1,xlograt)))

subplot(1,2,2)
plot(x1,x1.*exp(xlograt),'k.',p1,p2,'rv',mean(d1),mean(d2),'co')
xlabel('Isotope 1');ylabel('Isotope 2')
title(sprintf('%.1f%% kept',100*kept/M));
legend(sprintf('Corr = %.2f',corr(x1,x1.*exp(xlograt))))

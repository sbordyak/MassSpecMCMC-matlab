function [mean_of_logs, std_of_logs] = MeanofLogs(a,b,s2);


d = [a; b];

m0guess = real(2*log(mean(a./b))); % For initial guess, log of means

m0 = [m0guess; 1]; % [log(a/b); log(b)]
% a = exp( log(a/b) + log(b) )
% b = exp( log(b) )

% check that model parameterization works
dhat = evalg(a, b, m0);

% best fit to data using fminunc
opts = optimoptions("fminunc");
opts.StepTolerance = 1e-10;
opts.Display = 'off';
[x, ~] = fminunc(@(m) objfunc(d, a, b, s2, m), m0, opts);

mean_of_logs = x(1);

% %calculate uncertainty
% G = makeG(a, b, x);
% covx = inv( (G.*(1./s2))' * G ); % inv(G'*W*G), where W = diag(1./s2);
% std_of_logs = sqrt(diag(covx));
% std_of_logs = std_of_logs(1);

std_of_logs = NaN; % This is still giving me some problems %sb726

% x(1) = alrab -- the best fit is the alr mean!
% unctx(1) == logStde when all(a > 0 & b > 0).  Awesome!

% Note after coding: this all works out!  No logs of negative ratios
% due to negative measured intensities (a and b)!! Best fit value 
% of model parameters converges to geometric mean/alr results
% and uncertainty converges to standard error of logratios
% when all(a > 0 & b > 0)



%% functions

function G = makeG(a, b, m)

    arange = 1:length(a);
    brange = max(arange)+1 : max(arange) + length(b);
    
    lograbtrue = m(1);
    logbtrue = m(2);

    G(arange,1) = exp(lograbtrue + logbtrue);
    G(arange,2) = exp(lograbtrue + logbtrue);
    G(brange,1) = zeros(size(b));
    G(brange,2) = exp(logbtrue);

end % fuction makeG


% evaluate model function g(m), 
% return predicted data vector dhat
function dhat = evalg(a, b, m)

    arange = 1:length(a);
    brange = max(arange)+1 : max(arange) + length(b);
    
    dhat = zeros(size([a; b]));
    dhat(arange) = exp(m(1) + m(2));
    dhat(brange) = exp(m(2));

end % function evalg


% objective function - weighted sum of squares for d-dhat
function fval = objfunc(d, a, b, s2, m)

    dhat = evalg(a, b, m);
    r = d - dhat;
    fval = sum( r.^2 ./ s2);

end %function objfunc

end
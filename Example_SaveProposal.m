
N = 20; % Number of parameters

tmp = randn(N);
Sigma0 = tmp'*tmp+0.1*eye(N); 


filename = 'test.csv';

for ii = 1:100

    dx = mvnrnd(zeros(1,N),Sigma0,1);
    
    
    dlmwrite(filename,dx,'delimiter',',','-append');
    
end
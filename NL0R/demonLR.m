% demon sparse logistic regression problems 
clc; clear; close all;

n    = 10000;  
m    = ceil(n/5); 
s    = ceil(0.05*n);
test = 2;

switch test
    case 1  % You could input data including (data.A, data.At, data.b), e.g.,
         I0       = randperm(m);  
         I        = I0(1:ceil(m/2)); 
         b        = ones(m,1);    
         b(I)     = 0;
         data.A   = repmat(b.*rand(m,1),1,n) + randn(m,n);
         data.At  = data.A';    
         data.b   = b;
         pars.lam = 0.01;
    case 2  % Or you could input data by our data generation function 
         ExMat    = 1;
         MatType  = {'Correlated','Weakly-Indipendent'};
         data     = logistic_random_data(MatType{ExMat},m,n,s,0.5); 
         pars.lam = 1e-3;
end
 
pars.draw = 1;
pars.rate = (n>=1e3)*0.5+(n<1e3)*0.75;
out       = NL0R('LR',data,n,pars);

fprintf(' CPU time:      %.3fsec\n',  out.time);
fprintf(' Logistic Loss: %5.2e\n', out.obj);
fprintf(' Sample size:   %4dx%4d\n', m,n);

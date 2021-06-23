% demon sparse compressed sensing problems 
clc; clear; close all; 

n     = 10000; 
m     = ceil(0.25*n);
s     = ceil(0.05*n);     
test  = 2; 

switch test
  case 1       % Input data by our data generation function
       ExMat    = 1;
       MatType  = {'GaussianMat','PartialDCTMat'}; 
       nf       = 0.00;
       data     = CSdata(MatType{ExMat},m,n,s,nf);
       xopt     = data.xopt;
       if nf    > 0; pars.tau = 0.5; end
  case 2       % Input any data including (data.A, data.b), e.g.,
       I        = randperm(n); 
       Tx       = I(1:s);
       xopt     = zeros(n,1);  
       xopt(Tx) = randn(s,1); 
       data.A   = randn(m,n)/sqrt(m);
       data.b   = data.A*xopt;  
end  
pars.tol  = 1e-6; 
func      = @(x,T1,T2)CS(x,T1,T2,data);
pars.rate = (n>=1e3)*0.25+(n<1e3)*0.5;  
out       = NL0R(func,n,pars);

fprintf(' Sample size:       %dx%d\n', m,n);
fprintf(' CPU time:          %.3fsec\n',  out.time);
fprintf(' Sparsity:          %.2d\n', nnz(out.sol));
fprintf(' Objective:         %5.2e\n',  out.obj);
RecoverShow(xopt,out.sol,[900,500,500,250],1)




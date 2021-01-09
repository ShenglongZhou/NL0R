% demon sparse linear complementarity problems
clc; close all; clear all; warning off

n         = 5000;
r         = 0.01;
s         = ceil(r*n);
data      = lcp_data('sdp',n,s);
 
pars.rate = 0.25*(r<0.05)+0.5*(r>=0.05);
pars.tol  = 1e-15;
out       = NL0R('LCP',data,n,pars); 

fprintf(' Sample size:       %dx%d\n', n,n);
fprintf(' CPU time:          %.3fsec\n',  out.time);
fprintf(' Sparsity:          %.2d\n', nnz(out.sol));
fprintf(' Objective:         %5.2e\n',  out.obj);
if isfield(data,'xopt')
   fprintf(' Accuracy:          %5.2e\n',...
   norm(out.sol-data.xopt)/norm(data.xopt)); 
   if  s<=100
       ReoveryDisplay(data.xopt,out.sol,[900,500,500,250],1)
   end
end
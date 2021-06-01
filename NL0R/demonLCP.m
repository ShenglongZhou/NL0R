% demon sparse linear complementarity problems
clc; close all; clear all; warning off

n         = 5000;
r         = 0.01;
s         = ceil(r*n);
data      = LCPdata('sdp',n,s);
 
pars.rate = 0.25*(r<0.05)+0.5*(r>=0.05);
pars.tol  = 1e-10;
pars.xopt = data.xopt;
func      = @(x,T1,T2)LCP(x,T1,T2,data);
out       = NL0R(func,n,pars); 

fprintf(' Sample size:       %dx%d\n', n,n);
fprintf(' CPU time:          %.3fsec\n',  out.time);
fprintf(' Sparsity:          %.2d\n', nnz(out.sol));
fprintf(' Objective:         %5.2e\n',  out.obj);


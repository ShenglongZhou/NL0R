% demon sparse compressed sensing problems 
clc; clear; close all; 

n         = 500; 
m         = ceil(0.25*n);
s         = ceil(0.05*n);     

data      = compressed_sensing_data('GaussianMat',m,n,s,0); 
pars.rate = (n>=1e3)*0.25+(n<1e3)*0.5;  
out       = NL0R('CS',data,n,pars);

fprintf(' Sample size:       %dx%d\n', m,n);
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



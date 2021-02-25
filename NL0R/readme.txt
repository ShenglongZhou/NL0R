% To implement this solver, please 
% 
% [1] run startup.m first to add the path;
% [2] run demonXXXX.m to solve  different problems
% 
% This source code contains the algorithm described in:
% `S.Zhou, L.Pan and N.Xiu, 2020, Newton Method for l_0 Regularized Optimization'
% Please give credits to this paper if you use the code for your research.
% 
% ===================================================================
% The citation of the solver  NL0R takes the form of
%               
%                 out = NL0R(problem,data,n,pars)
% 
% It aims at solving the sparsity constrained optimization with form
% 
%         min_{x\in R^n} f(x) + \lambda \|x\|_0 
% 
% where \lambda>0 is updated iteratively.
% 
% Inputs:
%     problem:  A text string for different problems to be solved, (required)
%               = 'CS',  compressed sensing problems
%               = 'LCP', linear complementarity problems
%               = 'LR',  sparse logistic regression problems
%     data:     A triple structure (data.A, data.At, data.b) (required)
%               data.A, the measurement matrix, or a function handle @(x)A(x);
%               data.At = data.A',or a function handle @(x)At(x);
%               data.b, the observation vector 
%     n:        Dimension of the solution x, (required)             
%     pars:     Parameters are all OPTIONAL
%               pars.x0      --  Starting point of x (default, zeros(n,1))               
%               pars.lam     --  A positive scalar (default, 1e-3)
%               pars.rate    --  A positive scalar to adjust lam, (default, (n>=1e3)*0.25+(n<1e3)*0.75) 
%               pars.tau     --  A positive scalar (default, 1)
%               pars.disp    --  Display results for each iteration if pars.disp=1 (default)
%                                Don't display results for each iteration if pars.disp=1
%               pars.draw    --  A  graph will be drawn if pars.draw=0 (default) 
%                                No graph will be drawn if pars.draw=1 
%               pars.maxit   --  Maximum number of iterations, (default,2000) 
%               pars.tol     --  Tolerance of the halting condition, (default,1e-6)
%               pars.obj     --  An upper bound of f(x), (default,1e-10)
% Outputs:
%     out.sol:           The sparse solution x
%     out.sparsity:      Sparsity level of out.sol
%     out.error:         Error used to terminate this solver 
%     out.time           CPU time
%     out.iter:          Number of iterations
%     out.obj:           Objective function value at out.sol 
% 
% This code is programmed based on the algorithm proposed in 
% S. Zhou, L. Pan and N. Xiu, 2020, 
% Newton Method for l_0 Regularized Optimization.
% Send your comments and suggestions to <<< shenglong.zhou@soton.ac.uk >>> 
% Warning: Accuracy may not be guaranteed !!!!! 

% Here are some examples that you can run
% =================================================================
% Example I:  compressed sensing problem

n         = 2000; 
m         = ceil(0.25*n);
s         = ceil(0.01*n);     
x         = zeros(n,1);
I         = randperm(n);
x(I(1:s)) = randn(s,1);
data.A    = randn(m,n)/sqrt(n);
data.At   = data.A'
data.b    = data.A*x 
pars.rate = (n>=1e3)*0.25+(n<1e3)*0.75;  
out       = NL0R('CS',data,n,pars) 
ReoveryShow(out.sol,x,[900,500,500,250],1)

% =================================================================
% Example II:  linear complementarity problem 

n         = 2000; 
s         = ceil(0.01*n);     
x         = zeros(n,1);
I         = randperm(n); 
T         = I(1:s);
x(T)      = rand(s,1);
A         = randn(n,ceil(n/4));
data.A    = A*A'/n; 
data.At   = data.A;; 
Ax        = data.A*x;
data.b    = abs(Ax); 
data.b(T) = -Ax(T); 
data.n    = n;
pars.rate = (n>=1e3)*0.25+(n<1e3)*0.75;  
out       = NL0R('LCP',data,n,pars) 
ReoveryShow(out.sol,x,[900,500,500,250],1)

% =================================================================
% Example III:  Logistic regression problem

n         = 2000; 
m         = ceil(0.25*n);
s         = ceil(0.05*n);     
I         = randperm(n);
T         = I(1:s); 
data.A    = randn(m,n); 
data.At   = data.A'; 
q         = 1./(1+exp(-data.A(:,T)*randn(s,1)));
data.b    = zeros(m,1);
for i     = 1:m    
data.b(i) = randsrc(1,1,[0 1; 1-q(i) q(i)]);
end               
pars.lam  = 1e-3;
pars.draw = 1;
pars.rate = (n>=1e3)*0.5+(n<1e3)*0.75;
out       = NL0R('LR',data,n,pars) 

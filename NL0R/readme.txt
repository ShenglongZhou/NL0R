% To implement this solver, please 
% 
% [1] run startup.m first to add the path;
% [2] run demonXXXX.m to solve  different problems
% 
% This source code contains the algorithm described in:
% `S.Zhou, L.Pan and N.Xiu, Numerical Algorithm, 2021
%  Newton Method for l0 Regularized Optimization,'
% Please give credits to this paper if you use the code for your research.
% 
% =========================================================================
%The citation of the solver  NL0R takes the form of
%               
%                 out = NL0R(problem,data,n,pars)
% 
% It aims at solving the L_0-regularized sparse optimization with form
% 
%         min_{x\in R^n} f(x) + \lambda \|x\|_0 
% 
% where \lambda is updated iteratively.
% =========================================================================
% Inputs:
%     func:     A function handle defines (objective,gradient,sub-Hessain) (required)
%     n:        Dimension of the solution x (required)           
%     pars:     Parameters are all OPTIONAL
%               pars.x0      --  Starting point of x (default, zeros(n,1))
%               pars.tau     --  A positive scalar (default, (n<=1e3)+(n>1e3)/2)
%               pars.lam     --  An initial penalty parameter (default, maxlam/2)
%               pars.rate    --  A positive scalar to adjust lam, (default, rate0) 
%               pars.disp    --  Display results or not for each iteration (default, 1) 
%               pars.draw    --  Draw or not draw a graph (default, 0)  
%               pars.maxit   --  Maximum number of iterations, (default,2000) 
%               pars.tol     --  Tolerance of the halting condition, (default,1e-6)
%               pars.obj     --  A predefined lower bound of f(x), (default,1e-20)
% Outputs:
%     out.sol:           The sparse solution x
%     out.sparsity:      Sparsity level of out.sol
%     out.time           CPU time
%     out.iter:          Number of iterations
%     out.obj:           Objective function value at out.sol 
% ========================================================================= 
% This code is programmed based on the algorithm proposed in 
% S. Zhou, L. Pan and N. Xiu, Numerical Algorithms, 2021, 
% Newton Method for l_0 Regularized Optimization, Numerical Algorithm.
% Send your comments and suggestions to <<< slzhou2021@163.com >>> 
% Warning: Accuracy may not be guaranteed !!!!! 
% =========================================================================


% Here are some examples that you can run
% =================================================================
% Example I:  compressed sensing problem

n         = 10000; 
m         = ceil(0.25*n);
s         = ceil(0.01*n);     
x         = zeros(n,1);
I         = randperm(n);
x(I(1:s)) = randn(s,1);
data.A    = randn(m,n)/sqrt(n);
data.b    = data.A*x; 
func      = @(x,T1,T2)CS(x,T1,T2,data);
pars.rate = (n>=1e3)*0.25+(n<1e3)*0.5;  
out       = NL0R(func,n,pars);
RecoverShow(out.sol,x,[900,500,500,250],1)

% =================================================================
% Example II:  linear complementarity problem 

n         = 10000; 
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
pars.rate = 0.25*(s/n<0.05)+0.5*(s/n>=0.05);
func      = @(x,T1,T2)LCP(x,T1,T2,data);
out       = NL0R(func,n,pars); 
RecoverShow(out.sol,x,[900,500,500,250],1)

% =================================================================
% Example III:  Logistic regression problem

n         = 10000; 
m         = ceil(0.25*n);
s         = ceil(0.05*n);     
I         = randperm(n);
T         = I(1:s); 
data.A    = randn(m,n); 
q         = 1./(1+exp(-data.A(:,T)*randn(s,1)));
data.b    = zeros(m,1);
for i     = 1:m    
data.b(i) = randsrc(1,1,[0 1; 1-q(i) q(i)]);
end               
pars.draw = 1;
pars.rate = (n>=1e3)*0.5+(n<1e3)*0.75;
func      = @(x,T1,T2)LogitReg(x,T1,T2,data);
out       = NL0R(func,n,pars);

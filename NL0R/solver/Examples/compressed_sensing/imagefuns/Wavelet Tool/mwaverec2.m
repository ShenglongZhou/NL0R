function x = mwaverec2(c,s,rlo,rhi,md,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
% Multi-level 2-D wavelet reconstruction  of structure [c,s] by using     %
% filter rlo and rhi                                                      %
% (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)                          %
%  Created on Oct 17, 2011                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                   %
%                 c:  [ A(N)   | H(N)   | V(N)   | D(N) | ...             %
%                  H(N-1) | V(N-1) | D(N-1) | ...  | H(1) | V(1) | D(1) ].%
%                  where A, H, V, D, are row vectors such that:           %
%                  A = approximation coefficients, is bjll                %
%                  H = hori. detail coefficients,  is xjlh                %
%                  V = vert. detail coefficients,  is xjhl                %
%                  D = diag. detail coefficients,  is xjhh                %                                   
%             s(1,:): size of app. coef.(n)                               %
%             s(i,:): size of det. coef.(n-i+2) for i = 2,...,n+1         %
%           s(n+2,:): size(x).                                            %
%                 n:  0 <= n <= length(l)-2                               %
%               rlo:  the low frequency mask of  reconstruction           %   
%               rhi:  the high frequency mask of reconstruction           %
%                     using  [lo,hi,rlo,rhi] = wfilters('wavename')       %
%                     to get the filer of wavelet wavename                %
%                md:  the mehtod of boundry condition                     %    
%                     md=0  use mode 'ppd'                                %
%                     md=1  use mode 'per'                                %
%                     md=2  use mode 'sym'                                %
%                     md=3  use mode 'zpd'                                %
% Output                                                                  %
%                 x:  the signal to be reconstructed                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = mappcoef2(c,s,rlo,rhi,md,n);

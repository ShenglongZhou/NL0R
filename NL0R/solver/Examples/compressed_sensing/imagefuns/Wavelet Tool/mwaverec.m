function x = mwaverec(c,l,rlo,rhi,md,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
% Multi-level 1-D wavelet reconstruction  of structure [c,l] by using     %
% filter rlo and rhi                                                      %
% (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)                          %
%  Created on Oct 17, 2011                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                   %
%                  c:  [app. coef.(n)|det. coef.(n)|... |det. coef.(1)]   %
%               l(1):  length of app. coef.(n)                            %
%               l(i):  length of det. coef.(n-i+2) for i = 2,...,n+1      %
%             l(n+2):  length(x).                                         %
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
x = mappcoef(c,l,rlo,rhi,md,n);
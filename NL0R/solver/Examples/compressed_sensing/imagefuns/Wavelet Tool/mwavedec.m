
function [c,l] = mwavedec(x,n,lo,hi,md,ds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
% Multi-level 1-D wavelet decomposition of x by using filter lo and hi    %
% (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)                          %
%  Created on Oct 17, 2011                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                   %
%                 x:  the signal to be decomposed                         %
%                 n:  the level to be decomposed                          %
%                lo:  the low frequency mask of  decompostion             %   
%                hi:  the high frequency mask of decompostion             %
%                     using  [lo,hi,rlo,rhi] = wfilters('wavename')       %
%                     to get the filer of wavelet wavename                %
%                md:  the mehtod of boundry condition                     %    
%                     md=0  use mode 'ppd'                                %
%                     md=1  use mode 'per'                                %
%                     md=2  use mode 'sym'                                %
%                     md=3  use mode 'zpd'                                %
%                ds:  the down sample patten                              %
%                     ds = 1 down sample the even term                    %
%                     ds = 0 without down sample                          %
% Output                                                                  %
%              [c,l]                                                      %
%              The structure is organized as:                             %
%              c     = [app. coef.(n)|det. coef.(n)|... |det. coef.(1)]   %
%              l(1)   = length of app. coef.(n)                           %
%              l(i)   = length of det. coef.(n-i+2) for i = 2,...,n+1     %
%              l(n+2) = length(x).                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization.
s = size(x); x = x(:)'; % row vector
if s(1)>1
    disp('erro x shoud be a vector')
end
c = [];
l = zeros(1,n+2);
if isempty(x) , return; end

l(end) = length(x);
for k = 1:n
    [x,d] = dwt1(x,lo,hi,md,ds); % decomposition 
    c     = [d c];               % store detail
    l(n+2-k) = length(d);        % store length
end

% Last approximation.
c = [x c];
l(1) = length(x); % x now is the app.coef



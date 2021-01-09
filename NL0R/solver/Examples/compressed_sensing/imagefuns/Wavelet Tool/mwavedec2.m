function [c,s] = mwavedec2(x,n,lo,hi,md,ds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
% Multi-level 2-D wavelet decomposition of x by using filter lo and hi    %
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
%              [c,s]                                                      %
%                  the structure of  c:                                   %
%                  c = [ A(N)   | H(N)   | V(N)   | D(N) | ...            %
%                  H(N-1) | V(N-1) | D(N-1) | ...  | H(1) | V(1) | D(1) ].%
%                  where A, H, V, D, are row vectors such that:           %
%                  A = approximation coefficients, is bjll                %
%                  H = hori. detail coefficients,  is xjlh                %
%                  V = vert. detail coefficients,  is xjhl                %
%                  D = diag. detail coefficients,  is xjhh                %
%                  the structure of  s:                                   %
%                  s(1,:) = size of app. coef.(n)                         %
%                  s(i,:) = size of det. coef.(n-i+2) for i = 2,...,n+1   %
%                  and s(n+2,:) = size(x).                                %
%   NOTE: When x represents an indexed image, then x as well  as the      %
%   arrays CA, CH, CV, CD are m-by-n matrices.                            %
%   When x represents a truecolor image, then they become   m-by-n-by-3   %
%   arrays. These arrays consist of three m-by-n matrices                 %
%   So the sizes of the vector c and the matrix s depend on the           %
%   type of analyzed image.                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization.
c = [];
sx =  size(x);
s = zeros(n+2,length(sx));

if isempty(x) 
    return; 
end

s(end,:) = size(x);
for i=1:n
    [x,h,v,d] = mdwt2(x,lo,hi,md,ds); % decomposition
    c = [h(:)' v(:)' d(:)' c];        % store details
    s(n+2-i,:) = size(x);             % store size
end

% Last approximation.
c = [x(:)' c];
s(1,:) = size(x);


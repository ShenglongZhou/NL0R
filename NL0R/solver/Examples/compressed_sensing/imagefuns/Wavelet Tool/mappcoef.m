function a = mappcoef(c,l,rlo,rhi,md,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
% Computes the approximation  coefficients at level n using the wavelet   %
% decomposition    structure [c,l]          (see MWAVEDEC)                %
% (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)                          %
%  Created on Oct 17, 2011                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                   %
%             [c,l]:  the decomposition  structure                        %
%                       (see MWAVEDEC)                                    %
%                 x:  the signal to be decomposed                         %
%               rlo:  the low frequency mask of reconstructon             %   
%               rhi:  the high frequency mask of reconstruction           %
%                     using [lo,hi,rlo,rhi] = wfilters('wavename')        %
%                     to get the filer of wavelet wavename                %
%                md:  the mehtod of boundry condition                     %    
%                     md=0  use mode 'ppd'                                %
%                     md=1  use mode 'per'                                %
%                     md=2  use mode 'sym'                                %
%                     md=3  use mode 'zpd'                                %
%                n:   the level of decompotion of the structure   [c,l]   %
%                     must be an integer such that 0 <= n <= length(L)-2. %
% Output                                                                  %
%                a:  the approximate coefficient at level n               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  

rmax = length(l);
nmax = rmax-2;
if nargin==5
    n=nmax;
end

if (n < 0) || (n > nmax) || (n ~= fix(n))
    error( 'Invalid level value n ');
end

% Initialization.
a = c(1:l(1));

% Iterated reconstruction.
imax = rmax+1;
for p = nmax:-1:n+1
    d = detcoef(c,l,p);                % extract detail
    a = idwt1(a,d,rlo,rhi,md,l(imax-p));
end

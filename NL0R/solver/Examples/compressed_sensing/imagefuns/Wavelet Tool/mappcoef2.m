function a = mappcoef2(c,s,rlo,rhi,md,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
% Computes the approximation  coefficients at level n using the wavelet   %
% decomposition    structure [c,l]          (see MWAVEDEC)                %
% (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)                          %
%  Created on Oct 17, 2011                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                   %
%             [c,s]:  the decomposition  structure                        %
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
%                     must be an integer such that 0 <= n <= size(s,1)-2. %
% Output                                                                  %
%                a:  the approximate coefficient at level n               %
%   NOTE: If c and s are obtained from an indexed image analysis          %
%   (respectively a truecolor image analysis) then A is an                %
%   m-by-n matrix (respectively  an m-by-n-by-3 array).                   %
%   For more information on image formats, see the reference              %
%   pages of IMAGE and IMFINFO functions.                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmax = size(s,1);
nmax = rmax-2;

if (n<0) || (n > nmax) || (n~=fix(n))
    error( 'Invalid level value n ');
end

% Initialization.
nl   = s(1,1);
nc   = s(1,2);
if length(s(1,:))<3 , dimFactor = 1; else dimFactor = 3; end;
a    = zeros(nl,nc,dimFactor);
a(:) = c(1:nl*nc*dimFactor);

% Iterated reconstruction.
rm   = rmax+1;
for p=nmax:-1:n+1
    [h,v,d] = detcoef2('all',c,s,p);       % use matlab function detcoef2
    a = midwt2(a,h,v,d,rlo,rhi,md,s(rm-p,:));
end



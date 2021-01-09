function sx = idwt1(bj,xj,rlo,rhi,md,l) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
% One level 1-dimension   inverse wavelet transform of bj and xj          %
% by using the reconstruction  filter rlo and rhi get origial singal sx.  %
% (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)                          %
%  Created on Oct 17, 2011                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                   %
%                bj:  the low frequency part                              %
%                xj:  the high frequency part                             %
%               rlo:  the low frequency mask of reconstructon             %   
%               rhi:  the high frequency mask of reconstruction           %
%                     using [lo,hi,rlo,rhi] = wfilters('wavename')        %
%                     to get the filer of wavelet wavename                %
%                 l:  get the central l elements  of sx  when    md<>1,   %
%                     and get the first l elements of sx when md =1       %
%                md:  the mehtod of boundry condition                     %    
%                     md=0  use mode 'ppd'                                %
%                     md=1  use mode 'per'                                %
%                     md=2  use mode 'sym'                                %
%                     md=3  use mode 'zpd'                                %
%                ds:  the down sample patten                              %
%                     ds = 1 down sample the even term                    %
%                     ds = 0 without down sample                          %
% Output                                                                  %
%                sx: the reconstruciton signal                            %
%                    the length of sx is                                  %
%                    2*length(bj) when  md = 1;                           %
%                    2*length(bj)-length(rlo)+2 when  md<>1;               % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% get sx by the main function mupsconv1 
lx=[];
if nargin==6
    lx=l;
end
sx = mupsconv1(bj,rlo,md,lx) + mupsconv1(xj,rhi,md,lx);


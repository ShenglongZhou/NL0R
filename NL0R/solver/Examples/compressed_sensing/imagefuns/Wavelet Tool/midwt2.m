function sx = midwt2(bjll,xjlh,xjhl,xjhh,rlo,rhi,md,l)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% One level 2-dimension inverse wavelet transform of bjll,xjlh,xjhl,xjhh  %
% by using the reconstruction  filter rlo rhi. And 2d filter rll rlh rhl  %
% rhh and rlo are gotten by 1-d tensor of rlo (rhi) and rlo (rhi)         %
% (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)                          %
%  Created on Oct 17, 2011                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                   %
%              bjll:  approximation coeffiecient                          %
%              xjlh:  detail coeffiecient LH                              %
%              xjhl:  detail coeffiecient HL                              %
%              xjhh:  detail coeffiecient HH                              %
%               rlo:  the low frequency mask of reconstruction            %   
%               rhi:  the high frequency mask of reconstruction           %
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
%                 l:  get the central l elements  of sx  when    md<>1,   %
%                     and get the first l elements of sx when md =1       %
% Output                                                                  %
%                sx:  the reconstructed signal                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get x by the main function mupsconv2
sx=[];
if nargin==8
    sx=l;
end
sx =mupsconv2(bjll,rlo,rlo,md,sx)+ ... % Approximation.rbjll
    mupsconv2(xjlh,rhi,rlo,md,sx)+ ... % Horizontal Detail. rxjhl  which is inverse of direct transfrom  xjlh(lh)--rxjhl(hl)
    mupsconv2(xjhl,rlo,rhi,md,sx)+ ... % Vertical Detail.   rxjlh 
    mupsconv2(xjhh,rhi,rhi,md,sx);     % Diagonal Detail.   rxjhh
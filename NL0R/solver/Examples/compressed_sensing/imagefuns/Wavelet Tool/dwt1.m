function [bj,xj]=dwt1(x,lo,hi,md,ds) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
% One level 1-dimension  wavelet transform of x by using filter lo and hi %
% get low frequency part bj and high frequency part xj                    %
% (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)                          %
%  Created on Oct 17, 2011                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                   %
%                 x:  the signal to be decomposed                         %
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
%                bj:  approximation coeffiecient                          %
%                xj:  detail coeffiecient                                 %
%                   the length of bj and xj is                            %
%                   ceil(length(x)/2) when   ds = 1 and md=1;             %
%                   length(x) when ds = 0 and  md=1;                      %         
%                   floor((length(x)+length(lo)-1)/2) when ds =1 and md<>1;%
%                   floor(length(x)+length(lo)-1) when ds =0 and md<>1.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lf = length(lo);    % lf must be even
lx = length(x);

% Extend, Decompose &  Extract coefficients.
if md==0
    lenEXT = lf-1; last = lx+lf-1;
    y = wextend('1D','ppd',x,lenEXT);
elseif  md==1
    lenEXT = lf/2; 
    last = 2*ceil(lx/2);
    y = wextend('1D','per',x,lenEXT);
elseif md==2 
    lenEXT = lf-1; last = lx+lf-1;
    y = wextend('1D','sym',x,lenEXT); % this is the defout mode used in Matlab
else
    lenEXT = lf-1; last = lx+lf-1;
    y = wextend('1D','zpd',x,lenEXT);   % for md >=3 use 'zpd'
end

% Compute coefficients of approximation bj  and detail xj. 
z1 = wconv1(y,lo,'valid'); 
if (length(z1)-length(x))~=1
    disp('error')
end
z2 = wconv1(y,hi,'valid');
if ds==0  % with out down sample
    bj = z1(1:last);
    xj = z2(1:last);
elseif ds==1   % down sample 
         bj = z1(2:2:last);
         xj = z2(2:2:last);
else
    disp('erro ! ds must be 0 or 1')
end





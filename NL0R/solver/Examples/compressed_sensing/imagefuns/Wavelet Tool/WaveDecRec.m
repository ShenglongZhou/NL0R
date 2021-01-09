function [Y,strc] = WaveDecRec(X,level,lo,hi,rlo,rhi,md,ds,isdec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  Multi-level 2-D wavelet decomposition and reconstruction               %
% (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)                          %
%  Created on Oct 17, 2011                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                   %
%                 x:  the signal or coefficient                           %
%             level:  decompostion level                                  %
%                lo:  the low frequency mask of  decomposition            %   
%                hi:  the high frequency mask of decomposition            %
%               rlo:  the low frequency mask of  reconstruction           %   
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
%             isdec:  1 to do decomposition                               %
%                     0 to do reconstruction                              %
% Output                                                                  %
%                 Y:  the coefficient or signal                           %
%              strc:  the structure of the coefficient                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent s;
if isdec == 1
   [Y,s] = mwavedec2(X,level,lo,hi,md,ds);
    strc = s;
    Y = Y';
else
    Y  = mwaverec2(X(:),s,rlo,rhi,md,0);
    strc = s;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example 
% imgsize = 64;
% Img = phantom(imgsize);
% [my,mx] = size(Img);
% Imgbar = mean(Img(:));
% Imgcntr = Img -Imgbar;
% level = 2;                   % levels of wavelet transform
% wav = 'db1';                 % type of wavelet 
% [lo,hi,rlo,rhi] = wfilters(wav);
% W = @(z) WaveDecRec(z,level,lo,hi,rlo,rhi,1,1,0);  % Reconstuction operator
% Wt = @(z) WaveDecRec(z,level,lo,hi,rlo,rhi,1,1,1); % Decomopositon operator
% [xe,c] =  Wt(Imgcntr);
%  a = c(1,:);
% xt = zeros(size(xe));
% id = find(abs(xe)>1e-3);
% xt(id)=xe(id);
% Imr = W(xt);
% imshow(Imr+Imgbar)
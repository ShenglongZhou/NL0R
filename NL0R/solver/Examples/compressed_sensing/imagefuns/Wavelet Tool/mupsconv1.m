
function y = mupsconv1(x,f,md,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
% (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)                          %
%  Created on Oct 17, 2011                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lx = 2*length(x);
lf = length(f);
if isempty(s)
   if md==1 
        s = lx; 
   else
        s = lx-lf+2; 
   end
end
% Compute Upsampling and Convolution.
y = x;
if md==1
    y = dyadup(y,0,1);%%%
    y = wextend('1D','per',y,lf/2);    
    y = wconv1(y,f);
    y = y(lf:lf+s-1);
else
    y = wconv1(dyadup(y,0),f); % conv sape full
    y = wkeep1(y,s,'c');
end


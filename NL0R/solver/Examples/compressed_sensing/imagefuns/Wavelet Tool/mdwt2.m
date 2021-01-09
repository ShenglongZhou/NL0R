function [bjll,xjlh,xjhl,xjhh] = mdwt2(x,lo,hi,md,ds) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
% One level 2-dimension wavelet transform of x which is a matrix size m*n %
% or m*n*3 by using filter lo and hi to get low frequency part bjll and   %
% high  frequency part xjlh xjlh xjhh. And 2d filter ll lh hl hh lo are   %
% gotten by one dimension tensor of lo (hi) and lo (hi);                  %
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
%              bjll:  approximation coeffiecient                          %
%              xjlh:  detail coeffiecient LH                              %
%              xjhl:  detail coeffiecient HL                              %
%              xjhh:  detail coeffiecient HH                              %
%                   the length of bj and xj is                            %
%                   ceil(length(x)/2) when   ds = 1 and md=1;             %
%                   length(x) when ds = 0 and  md=1;                      %         
%                   floor((length(x)+length(lo)-1)/2) when ds =1 and md<>1;%
%                   floor(length(x)+length(lo)-1) when ds =0 and md<>1.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute sizes.
lf = length(lo); % lf must be even when using orthogal wavelet transform 
sx = size(x);

% Extend, Decompose &  Extract coefficients.
if md==1
    sizeEXT = lf/2; 
    last = 2*ceil(sx/2);
else 
    sizeEXT = lf-1; last = sx+lf-1;
end

x = double(x); % made x to be double 
if ds==1
          if length(sx)==2
              if md==0
                   y = wextend('addcol','ppd',x,sizeEXT);
              elseif md==1
                y = wextend('addcol','per',x,sizeEXT); 
               
              elseif md==2
                   y = wextend('addcol','sym',x,sizeEXT); 
              else
                   y = wextend('addcol','zpd',x,sizeEXT); % for md >=3 use 'zpd'
              end
                z = conv2(y,lo(:)','valid');         % lo(:)' which is a row vector and z is got by conv y with Lo_D(:)' one row by one row 
                bjll = convds(z,lo,sizeEXT,last,md); % this shoud do conv of z with lo(which is stored in row vector) downsapled and then  one colum by one colum    
                xjlh = convds(z,hi,sizeEXT,last,md); % horizen direction  low vertical direciton high 
                z = conv2(y,hi(:)','valid');
                xjhl = convds(z,lo,sizeEXT,last,md); % horizen direction  high vertical direciton low  
                xjhh = convds(z,hi,sizeEXT,last,md);
          else
                a = cell(0,3);
                h = cell(0,3);
                v = cell(0,3);
                d = cell(0,3);
                for k = 1:3
                    if md==0
                            y = wextend('addcol','ppd',x(:,:,k),sizeEXT);
                    elseif md==1
                            y = wextend('addcol','per',x(:,:,k),sizeEXT); 
               
                    elseif md==2
                            y = wextend('addcol','sym',x(:,:,k),sizeEXT); 
                    else
                            y = wextend('addcol','zpd',x(:,:,k),sizeEXT); % for md >=3 use 'zpd'
                    end
                           
                           z = conv2(y,lo(:)','valid');
                           a{k} = convds(z,lo,sizeEXT,last,md);%fu zhi yong de shi {}
                           h{k} = convds(z,hi,sizeEXT,last,md);
                           z = conv2(y,hi(:)','valid');
                           v{k} = convds(z,lo,sizeEXT,last,md);
                           d{k} = convds(z,hi,sizeEXT,last,md);
                 end
                 bjll = cat(3,a{:});
                 xjlh = cat(3,h{:});
                 xjhl = cat(3,v{:});
                 xjhh = cat(3,d{:});  
          end
elseif ds==0
          if length(sx)==2
              if md==0
                   y = wextend('addcol','ppd',x,sizeEXT);
              elseif md==1
                y = wextend('addcol','per',x,sizeEXT); 
               
              elseif md==2
                   y = wextend('addcol','sym',x,sizeEXT); 
              else
                   y = wextend('addcol','zpd',x,sizeEXT); % for md >=3 use 'zpd'
              end
                z = conv2(y,lo(:)','valid');            % lo(:)' which is a row vector and z is got by conv y with Lo_D(:)' one row by one row 
                bjll = convndown(z,lo,sizeEXT,last,md); % this shoud do conv of z with lo(which is stored in row vector)  one colum by one colum    
                xjlh = convndown(z,hi,sizeEXT,last,md); % horizen direction  low vertical direciton high  
                z = conv2(y,hi(:)','valid');
                xjhl = convndown(z,lo,sizeEXT,last,md); % horizen direction  high vertical direciton low 
                xjhh = convndown(z,hi,sizeEXT,last,md);
          else
                a = cell(0,3);
                h = cell(0,3);
                v = cell(0,3);
                d = cell(0,3);
                for k = 1:3
                    if md==0
                            y = wextend('addcol','ppd',x(:,:,k),sizeEXT);
                    elseif md==1
                            y = wextend('addcol','per',x(:,:,k),sizeEXT); 
               
                    elseif md==2
                            y = wextend('addcol','sym',x(:,:,k),sizeEXT); 
                    else
                            y = wextend('addcol','zpd',x(:,:,k),sizeEXT); % for md >=3 use 'zpd'
                    end
                           z = conv2(y,lo(:)','valid');
                           a{k} = convndown(z,lo,sizeEXT,last,md); % fu zhi yong de shi {}
                           h{k} = convndown(z,hi,sizeEXT,last,md);
                           z = conv2(y,hi(:)','valid');
                           v{k} = convndown(z,lo,sizeEXT,last,md);
                           d{k} = convndown(z,hi,sizeEXT,last,md);
                 end
                 bjll = cat(3,a{:});
                 xjlh = cat(3,h{:});
                 xjhl = cat(3,v{:});
                 xjhh = cat(3,d{:});  
          end
else
    disp('erro ! ds must be 0 or 1')
end

%%
% Internal Function(s)
function y = convds(x,F,sizeEXT,last,md) 
y = x(:,2:2:last(2));
              if md==0
                   y = wextend('addrow','ppd',y,sizeEXT);
              elseif md==1
                y = wextend('addrow','per',y,sizeEXT); 
               
              elseif md==2
                   y = wextend('addrow','sym',y,sizeEXT); 
              else
                   y = wextend('addrow','zpd',y,sizeEXT); % for md >=3 use 'zpd'
              end
y = conv2(y',F(:)','valid');
y = y';
y = y(2:2:last(1),:);
%%
function y = convndown(x,F,sizeEXT,last,md) %conv no downsample
y = x(:,1:last(2));
              if md==0
                   y = wextend('addrow','ppd',y,sizeEXT);
              elseif md==1
                y = wextend('addrow','per',y,sizeEXT); 
               
              elseif md==2
                   y = wextend('addrow','sym',y,sizeEXT); 
              else
                   y = wextend('addrow','zpd',y,sizeEXT); % for md >=3 use 'zpd'
              end
y = conv2(y',F(:)','valid');
y = y';
y = y(1:last(1),:);
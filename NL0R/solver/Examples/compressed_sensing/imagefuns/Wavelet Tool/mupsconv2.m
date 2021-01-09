function y = mupsconv2(x,F1_R,F2_R,md,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
% (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)                          %
%  Created on Oct 17, 2011                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lf = length(F1_R);
sx = 2*size(x);

ndimX = ndims(x);

   if ndimX>2 , sx = sx(1:2); end

if isempty(s)
    if md==1 
           s = sx; 
    else
            s = sx-lf+2; 
    end
end
if ndimX<3
    y = mupsconv2ONE(x,F1_R,F2_R,md,s);
else
    y = cell(0,3);
    for j = 1:3
        y{j} = mupsconv2ONE(x(:,:,j),F1_R,F2_R,md,s);
    end
     y = cat(3,y{:});
end

    function y = mupsconv2ONE(z,F1_R,F2_R,md,s)
        % Compute Upsampling and Convolution.
        if md==1
            lf = length(F1_R);
            y = dyadup(z,'row',0,1);  
            y = wextend('addrow','per',y,lf/2);
            y = conv2(y,F1_R(:),'full'); 
            y = y(lf:lf+s(1)-1,:);
            y = dyadup(y,'col',0,1);
            y = wextend('addcol','per',y,lf/2);
            y = conv2(y,F2_R(:)','full');
            y = y(:,lf:lf+s(2)-1);
        else
            y = dyadup(z,'row',0);        % add even zero row since the direct transform is add colum and conv by row first  
            y = conv2(y,F1_R(:),'full');  % then the inverse transform  is add row and conv by colum first 
            y = dyadup(y,'col',0);
            y = conv2(y ,F2_R(:)','full');
            y = wkeep2(y,s,'c');
        end



function [X,Xt,y,ye,xe,A,out] = getrealdata(Img,nsam,sigma,seedno)

% imgsize    = 512
% nsam       = 60 number of samples 
% sigma      = 0.01 noise ratio
% matrixtype = 'real1d' or 'real2d'
 
rand('seed',seedno);   % fix seed
randn('seed', seedno); % fix seed

%Img   = phantom(imgsize);
dW_L  = 4;
out.I = Img;

[~,mx] = size(Img);
Imgbar = mean(Img(:));
Imgcntr = Img -Imgbar;
%  Sampling operator
[~,~,~,mhi] = LineMask(nsam,mx);
OMEGA = mhi;
Phi   = @(z) A_fhp(z,OMEGA);      % z is same size of original image, out put is a vector.
Phit  = @(z) At_fhp(z,OMEGA,mx); % z is a vector, out put is same size of original image.
% taking measurements
ye = Phi(Imgcntr);
M  = length(ye(:));
% Effective sensing operator
wav = 'db1';                 % type of wavelet
[lo,hi,rlo,rhi] = wfilters(wav);
X  = @(z) HDWIdW(Phi,Phit,z,lo,hi,rlo,rhi,1,1,dW_L,0);
Xt = @(z) HDWIdW(Phi,Phit,z,lo,hi,rlo,rhi,1,1,dW_L,1);
W  = @(z) WaveDecRec(z,dW_L,lo,hi,rlo,rhi,1,1,0);  % Rec
Wt = @(z) WaveDecRec(z,dW_L,lo,hi,rlo,rhi,1,1,1); % Dec
xe =  Wt(Imgcntr);
xe = xe(:);
A  = find(xe);
randn('state',0)
noise    = sigma*randn(M,1);
y        = ye + noise;
out.Ibar = Imgbar;
out.W    = W;
test=Xt(zeros(M,1)); clear test; 
end

function [M,Mh,mi,mhi] = LineMask(L,N)
% Returns the indicator of the domain in 2D fourier space for the
% specified line geometry.
% Usage :  [M,Mh,mi,mhi] = LineMask(L,N)
%
% Written by : Justin Romberg
% Created : 1/26/2004
% Revised : 12/2/2004
thc = linspace(0, pi-pi/L, L);
%thc = linspace(pi/(2*L), pi-pi/(2*L), L);
M = zeros(N);
% full mask
for ll = 1:L
    
    if ((thc(ll) <= pi/4) || (thc(ll) > 3*pi/4))
        yr = round(tan(thc(ll))*(-N/2+1:N/2-1))+N/2+1;
        for nn = 1:N-1
            M(yr(nn),nn+1) = 1;
        end
    else
        xc = round(cot(thc(ll))*(-N/2+1:N/2-1))+N/2+1;
        for nn = 1:N-1
            M(nn+1,xc(nn)) = 1;
        end
    end
    
end
% upper half plane mask (not including origin)
Mh = M;
Mh(N/2+2:N,:) = 0;
Mh(N/2+1,N/2+1:N) = 0;
M = ifftshift(M);
mi = find(M);
Mh = ifftshift(Mh);
mhi = find(Mh);
end

function y = A_fhp(x, OMEGA)
% Takes measurements in the upper half-plane of the 2D Fourier transform.
% x - N vector
% y - K vector = [mean; real part(OMEGA); imag part(OMEGA)]
% OMEGA - K/2-1 vector denoting which Fourier coefficients to use
%         (the real and imag parts of each freq are kept).
% Written by: Justin Romberg, Caltech and Modified by Yuling Jiao
[s1,s2] = size(x);
n = round(sqrt(s1*s2));
yc = 1/n*fft2(x);
y = [yc(1,1); sqrt(2)*real(yc(OMEGA)); sqrt(2)*imag(yc(OMEGA))];
end

function x = At_fhp(y, OMEGA, n)
% Adjoint of At_fhp (2D Fourier half plane measurements).
% y - K vector = [mean; real part(OMEGA); imag part(OMEGA)]
% OMEGA - K/2-1 vector denoting which Fourier coefficients to use
%         (the real and imag parts of each freq are kept).
% n - Image is nxn pixels
% x - N vector
% Written by: Justin Romberg, Caltech and modified by Yuling Jiao

K = length(y);
fx = zeros(n,n);
fx(1,1) = y(1);
fx(OMEGA) = sqrt(2)*(y(2:(K+1)/2) + 1i*y((K+3)/2:K));
x = real(n*ifft2(fx));
end

function [Y,strc] = WaveDecRec(X,level,lo,hi,rlo,rhi,emd,ds,isdec)
persistent s;
if isdec == 1
    % Phi' * X
    [Y,s] = mwavedec2(X,level,lo,hi,emd,ds);
    strc = s;
    Y = Y';
else
    Y  = mwaverec2(X(:),s,rlo,rhi,emd,0);
    strc = s;
end
end

function Y = HDWIdW(A,At,X,lo,hi,rlo,rhi,emd,ds,L,isdec)
% X: sparse vector which is taken from Wavelet transform of an image
% A: random projection matrix M x N
% rlo,rhi: scaling filter
% L: level of decomposition
% m, n: size of image
% Return Y: vector M x 1
% Written by Kun Qiu, ISU, April 2009 and modified by Yuling Jiao
persistent s
if isdec == 1 % Ht
    % converting measurements into samples
    if ~isa(At, 'function_handle')
        Y=At*X;
    else
        Y=At(X);
    end
    % converting samples into wavelet coefficients (sparse representation)
    [Y,s]= mwavedec2(Y,L,lo,hi,emd,ds);
    Y = Y';
else
    Y  = mwaverec2(X(:),s,rlo,rhi,emd,0);
    if ~isa(A, 'function_handle')
        Y = A*Y;
    else
        Y = A(Y);
    end
end
end
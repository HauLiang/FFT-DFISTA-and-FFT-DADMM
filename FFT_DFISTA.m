function x = FFT_DFISTA(PSF, D, b, x0, lambda, tol)
% This code implements the FFT-FISTA
% -- Fast Fourier Transform Differential Fast Iterative Shrinkage Thresholding Algorithm
%
%   minimize 1/2*|| F-1[F(x).*F(PSF)]-b]||_F^2 +lambda*||Dx||_2
%   s.t. x>=0
%
%   The solution is returned in the matrix x.
% 
% More information about FISTA can be found in the paper linked at:
% https://epubs.siam.org/doi/abs/10.1137/080716542
%
%
% Inputs:
%    PSF:  corresponding point spread function 
%    D:  modified first-order difference matrix 
%    b:  beamforming map, obtained by DAS
%    x0: initialization matrix
%    lambda:  parameter controlling the smoothness
%    tol:     tolerance of convergence criterion
% Outputs:
%    x:  high-resolution beamforming map, obtained by FFT-DFISTA
% Author: Hao Liang (haoliang@stu.xmu.edu.cn) 
% Last modified by: 21/05/16
%

% Initialize variables
x = x0; xold = x; y = x; t = 1; 

% Precompute fft of PSF
Fps = fft2(PSF); 

% Compute Lipschitz constant
L = lipschitz1(PSF,Fps,D,lambda);

% For n = 1
[~,grady] = gradient(Fps,b,y,D,lambda);

% Start iteration
MAX_ITER = 100; n = 0; xDif = 1;
while ( xDif > tol &&  n <= MAX_ITER ) 
    % update x
    x = max(0,y - (1/L)*grady);
    
    % updata t
    tnew = (1+sqrt(1+4*t*t))/2;
    
    % update y
    y = x + ((t-1)/tnew)*(x-xold);
    
    % gradient
    [~,grady] = gradient(Fps,b,y,D,lambda);
    
    % convergence judgment
    xDif = (norm(x - xold)/norm(xold))^2;
    xold = x; n = n+1;
    t = tnew;
end

end

function L = lipschitz1(PSF,Fps,D,lambda)
% Estimate Lipschitz constant by power iteration
x = rand(size(PSF));
for k = 1:10
    x = (fftshift(ifft2(fft2(x).*Fps))+lambda*D*x)/norm(x,'fro');
end
    L = norm(x,'fro')^2;    % lambda(A'A) Assuming a symmetric matric A
end

function [f,g,r] = gradient(Fps,b,x,D,lambda)
    r = fftshift(ifft2(fft2(x).*Fps)) - b;
    f = 0.5*norm(r,'fro')^2;
    if (nargout > 1)
        g = fftshift(ifft2(fft2(r).*Fps))+lambda*(D.'*D)*x;
    end
end


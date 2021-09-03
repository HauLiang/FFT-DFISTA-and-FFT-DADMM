function x = FFT_DADMM(PSF, D, b, x0, lambda, tol)
% This code implements the FFT-DADMM 
% -- Fast Fourier Transform Differential Alternating Direction Multiplier Method 
%
%   minimize 1/2*|| F-1[F(x).*F(PSF)]-b]||_F^2 +lambda*||Dx||_1
%   s.t. x>=0
%
%   The solution is returned in the matrix x.
% 
% More information about ADMM can be found in the paper linked at:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%
%
% Inputs:
%    PSF:  corresponding point spread function 
%    D:  modified first-order difference matrix 
%    b:  beamforming map, obtained by DAS
%    x0: initialization matrix
%    lambda:  parameter controlling the sparsity
%    tol:     tolerance of convergence criterion
% Outputs:
%    x:  high-resolution beamforming map, obtained by FFT-DADMM
%
% Author: Hao Liang (haoliang@stu.xmu.edu.cn)
% Last modified by: 21/05/16
%

% Global constants and defaults
MAX_ITER = 100;

% ADMM solver
x = x0; z = x; u = x; rho = 1;
xDif = 1; k = 1; xold = 0;

while ( xDif > tol &&  k <= MAX_ITER ) 
    % x-update
    x = FISTA_solver(x,PSF,b,z,D,rho,u);

    % z-update 
    x_hat = D*x;
    z = shrinkage(D*x + u, lambda/rho);

    % u-update
    u = u + (x_hat - z);
    
    % convergence judgment
    xDif = (norm(x - xold)/norm(xold))^2;
    xold = x; k = k + 1;
end

end

function z = shrinkage(x, kappa)
% soft threshold function
temp = abs(x);
if temp < kappa
    z = zeros(size(x));
else
    z = (temp>kappa).*(temp-kappa).*(x./(temp+eps));
end
end



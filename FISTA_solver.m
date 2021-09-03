function x = FISTA_solver(x,PSF,b,z,D,rho,u)
%   minimize 1/2*|| F-1[F(x).*F(PSF)]-b]||_F^2 + rho/2*||z-Dx+u||_F^2 + I+(x)
%
%   The solution is returned in the matrix x.

% Compute Lipschitz constant
Fps = fft2(PSF);
L = calculate_lipschitz(PSF,Fps,D,rho);
y = x;
t = 1;

% For n = 1
[~,grady] = calculate_gradient(Fps,b,x,u,D,z,rho);

% Start iteration
n = 0;
while n < 50  
    
    n = n+1; xold = x;
    x = max(0,y - (1/L)*grady);
    
    tnew = (1+sqrt(1+4*t*t))/2;
    y = x + ((t-1)/tnew)*(x-xold);
    
    [~,grady] = calculate_gradient(Fps,b,y,u,D,z,rho);
end

end

function L = calculate_lipschitz(PSF,Fps,D, rho)
    % Estimate Lipschitz constant by power iteration
    xx = rand(size(PSF));
    for k = 1:10
        xx = (fftshift(ifft2(fft2(xx).*Fps))+rho*D*xx)/norm(xx,'fro');
    end
    L = norm(xx,'fro')^2;  % lambda(A'A) Assuming a symmetric matric A
end

function [f,g] = calculate_gradient(Fps,b,x,u,D,z,rho)

r1 = fftshift(ifft2(fft2(x).*Fps)) - b;
r2 = z-D*x+u;

f = 0.5*norm(r1,'fro')^2 + 0.5/rho*norm(r2,'fro');

% gradient
g = fftshift(ifft2(fft2(r1).*Fps))+rho*(z-D*x-u);

end
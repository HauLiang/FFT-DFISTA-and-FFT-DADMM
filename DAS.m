function [b,PSF,hn] = DAS(N,z0,f,phi,rn,source,SNR)
% This code implements the Delay-and-Sum (DAS) algorithm
%
% More information about DAS can be found in the paper:
%    L. de Santana, "Fundamentals of Acoustic Beamforming," 
%    Design and Operation of Aeroacoustic Wind Tunnel Tests for Group and Air Transport, 2017.
%
%
% Inputs:
%    N:  number of grid points in each dim 
%    z0: source distance
%    f:  sampling frequency
%    rn: coordinates of the microphone array
%    source: x,y position of sources
%    SNR:  signal-to-noise ratio
% Outputs:
%    b:  beamforming map, obtained by DAS, also termed the 'dirty map'
%    PSF:  point spread function (PSF)
%    ej:   steering vector
%
% Author: Hao Liang 
% Last modified by: 21/09/03
%

% Number of microphones in array
M = size(rn,1);

% Parameters
c = 343;         % Speed of sound
omega = 2*pi*f;  % Angular frequency
Np = round(N*1.3);     % PSF grid size

% Parameters initialization 
dj = zeros(Np,Np,M); PSF = zeros(Np,Np);
hn = zeros(Np,Np,M); gn = zeros(Np,Np,M); 
b = zeros(N,N);

% Scan plane
L = 2*z0*tand(phi);            
x = [-L/2 L/2];    
rx = linspace(x(1),x(2),N);
[X,Y] = meshgrid(rx);

% PSF grid
rx_psf = linspace(x(1),x(2),Np);
[Xp,Yp] = meshgrid(rx_psf);

% |d0|: Distance from (0,0,0) to all mesh points
d0 = sqrt(Xp.^2 + Yp.^2 + z0^2);    

% Distance dj from each microphone to each grid point
for n = 1:M
    dj(:,:,n) = sqrt((Xp-rn(n,1)).^2+(Yp-rn(n,2)).^2 + z0^2);
    hn(:,:,n) = (dj(:,:,n)./d0).*exp(1j*omega.*dj(:,:,n)./c);
    gn(:,:,n) = (d0./dj(:,:,n)).*exp(1j*omega.*dj(:,:,n)./c);
end

% Point spread function (PSF)
% Constructed by the DAS beamformer response of a unit strength point
% source in the middel of the grid
ind = round(Np/2);
e_unit = squeeze(gn(ind,ind,:));

for ii = 1:length(Xp)
    for jj = 1:length(Yp)
        PSF(ii,jj) = dot(squeeze(hn(ii,jj,:)),e_unit);
    end
end

PSF = rot90(abs(PSF).^2/M^2);      % Normalize PSF
if N ~= Np
    PSF = interp2(Xp,Yp,PSF,X,Y);
end

% Cross spectral matrix (CSM)
dj = zeros(N,N,M);
hn = zeros(N,N,M);
gn = zeros(N,N,M);

% Distance and steering vectors
d0 = sqrt(X.^2 + Y.^2 + z0^2);    % |d0|: Distance from (0,0,0) to all mesh points
for n = 1:M
    dj(:,:,n) = sqrt((X-rn(n,1)).^2+(Y-rn(n,2)).^2 + z0^2);
    hn(:,:,n) = (dj(:,:,n)./d0).*exp(1j*omega.*dj(:,:,n)./c);
    gn(:,:,n) = (d0./dj(:,:,n)).*exp(1j*omega.*dj(:,:,n)./c)+10^(-SNR/10)*(rand(N,N)+1j*rand(N,N));
end

% Create CSM
q = zeros(N,N);
CSM = zeros(M,M);

for k = 1:size(source,1)
    q(source(k,2),source(k,1)) = 1;
    CSM = CSM + squeeze(gn(source(k,2),source(k,1),:))*squeeze(gn(source(k,2),source(k,1),:))';
end

% Diagonal removal technique
CSM(logical(eye(size(CSM)))) = 0;

% DAS imaging
for ii = 1:length(X)
    for jj = 1:length(Y)
        e = squeeze(hn(ii,jj,:));
        b(ii,jj) = dot(e,CSM*e)/(M^2-M);
    end
end

end

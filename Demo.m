% ---------------- Demo for FFT-DFISTA and FFT-DADMM ----------------------
%
% This is a simple example to test the FFT-DFISTA and FFT-DADMM algorithm 
% -- Fast Fourier Transform Differential Fast Iterative Shrinkage Thresholding Algorithm (FFT-DFISTA)
% -- Fast Fourier Transform Differential Alternating Direction Multiplier Method (FFT-DADMM)
%
% Author: Hao Liang 
% Last modified by: 21/05/16
%

%% Parameter Setting
clc; clear; close all;
load('D.mat');   % load modified first-order difference matrix 
load('56_spiral_array.mat');   % load 56-channel microphone spatial location
rn = array; % coordinates of the microphone array
N = 50;     % number of grid points in each dim
z0 = 5;     % source distance 
phi = 15;   % off-axis angle 
f = 1500;   % sampling frequency 
SNR = 15;   % signal-to-noise ratio 
source = int64([N/2-N/4 N/2; N/2+N/4 N/2]);    % x,y position of sources

%% Traditional Beamforming (DAS)
[b,PSF] = DAS(N,z0,f,phi,rn,source,SNR);

subplot(221)
contourf(real(b)); hold on; colormap(hot)
plot(source(:,1),source(:,2),'k*')
axis equal; title('DAS')


%% FFT-DFISTA and FFT-DADMM Algorithms
N = size(b,1); b = real(zeropad(b));
PSF = zeropad(PSF); x0 = zeros(2*N);
lambda = 10; tol = 5e-5;

% FFT-DFISTA
fprintf('\t------------------------------------------\n');
fprintf('\tStart FFT-DFISTA algorithm...\n');
tic;
x_dfista = FFT_DFISTA(PSF, D, b, x0, lambda, tol);
x_dfista = x_dfista(int64(N/2)+1:int64(N/2 + N),int64(N/2)+1:int64(N/2 + N));  % remove zero-padding
time_DFISTA = toc;
fprintf('\tElapsed time is %.2f seconds...\n',time_DFISTA);

subplot(222)
contourf(x_dfista,'LineStyle','none'); hold on; colormap(hot)
plot(source(:,1),source(:,2),'k*')
axis equal; title('FFT-DFISTA')


% FFT-DADMM:
fprintf('\t------------------------------------------\n');
fprintf('\tStart FFT-DADMM algorithm...\n');
tic;
x_dadmm = FFT_DADMM(PSF, D, b, x0, lambda, tol);
x_dadmm = x_dadmm(int64(N/2)+1:int64(N/2 + N),int64(N/2)+1:int64(N/2 + N));   % remove zero-padding
time_DADMM = toc;
fprintf('\tElapsed time is %.2f seconds...\n',time_DADMM);

subplot(223)
contourf(x_dadmm,'LineStyle','none'); hold on; colormap(hot)
plot(source(:,1),source(:,2),'k*')
axis equal; title('FFT-DADMM')


% FFT-DADMM using FFT-DFISTA to initialize
fprintf('\t------------------------------------------\n');
fprintf('\tStart FFT-DADMM algorithm (with FFT-DFISTA initialization)...\n');
tic;
new_x0 = FFT_DFISTA(PSF, D, b, x0, lambda, tol);
x_dadmm = FFT_DADMM(PSF, D, b, new_x0, lambda, tol);
x_dadmm = x_dadmm(int64(N/2)+1:int64(N/2 + N),int64(N/2)+1:int64(N/2 + N));   % remove zero-padding
time_init_DADMM = toc;
fprintf('\tElapsed time is %.2f seconds...\n',time_init_DADMM);
fprintf('\t------------------------------------------\n');

subplot(224)
contourf(x_dadmm,'LineStyle','none'); hold on; colormap(hot)
plot(source(:,1),source(:,2),'k*')
axis equal; title('FFT-DADMM (Initialization)')

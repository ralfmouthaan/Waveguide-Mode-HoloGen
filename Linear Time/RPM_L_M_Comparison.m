% Ralf Mouthaan
% University of Cambridge
% May 2020
%
% This is a comparison of the L and M matrices for the linear-time
% calculation for fiber mode hologram generation.
%
% It is noted that M is quite small. Consequently, its Fourier transform L
% is (kind of) big. We thus have to apply a very high threshold value to
% ensure that LT is more efficient than DS, but this would give a lot of
% artefacts.

clc; clear variables; close all;
    
%% User-Entered Variables

identifier = '25um Step';

f1 = 15e-3;
f2 = 75;
f3 = 100;
f4 = 10e-3;
Nx = 2500;
lambda = 633e-9;
gamma = 1;

HoloDiameter = 1000;
structDF.dx = 3.74e-6;
SMF_NA = 0.14;

structTargetLP.n_co = 1.457; % Silica
structTargetLP.n_cl = 1.4403;
structTargetLP.l=3; 
structTargetLP.m=3;
structTargetLP.lambda = lambda;
structTargetLP.WaveguideDiameter = 65e-6;
    
%% Diffraction field

structDF.x = linspace(-Nx/2,Nx/2,Nx);
[structDF.x_mesh, structDF.y_mesh] = meshgrid(structDF.x, structDF.x.');
structDF.r_mesh = sqrt(structDF.x_mesh.^2 + structDF.y_mesh.^2);
structDF.Mask = ones(size(structDF.r_mesh));
structDF.Mask(structDF.r_mesh > HoloDiameter/2) = 0;
structDF.x = structDF.x*structDF.dx;
structDF.x_mesh = structDF.x_mesh*structDF.dx;
structDF.y_mesh = structDF.y_mesh*structDF.dx;
structDF.r_mesh = structDF.r_mesh*structDF.dx;
structDF.Illumination = structDF.Mask.*exp(-structDF.r_mesh.^2/(f1*SMF_NA)^2);

%% Replay field coordinates

structRF.dx = lambda*f4/(Nx*structDF.dx*f2/f3);
structRF.x = linspace(-Nx/2,Nx/2,Nx);
structRF.x = structRF.x*structRF.dx;
[structRF.x_mesh, structRF.y_mesh] = meshgrid(structRF.x, structRF.x.');
structRF.r_mesh = sqrt(structRF.x_mesh.^2 + structRF.y_mesh.^2);
structRF.Mask = ones(size(structRF.x_mesh));
structRF.Mask(structRF.r_mesh > structTargetLP.WaveguideDiameter) = 0;
structRF.Mask = logical(structRF.Mask);

%% M and L Study

M = structRF.Mask;
L = fftshift(ifft2(fftshift(M)));
L = abs(L);
L = L/max(max(L));

figure;
subplot(2,2,1);
imagesc(M);
title('M');
axis square;

subplot(2,2,2);
imagesc(L);
title('L');
axis square;

thresh = 0.01;
fprintf('Non-zero M: %d\n', sum(sum(M>thresh)));
fprintf('Non-zero L: %d\n', sum(sum(L>thresh)));

L(L<thresh) = 0;
subplot(2,2,3);
imagesc(L);
title('Thresholded L');
axis square;

M = fftshift(fft2(fftshift(L)));
M = abs(M);
subplot(2,2,4);
imagesc(M);
title('Thresholded M');
axis square;

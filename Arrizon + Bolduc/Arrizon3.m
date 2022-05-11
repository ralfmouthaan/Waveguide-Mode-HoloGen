% Ralf Mouthaan
% University of Cambridge
% November 2020
% 
% Coding up of Arrizon's first algorithm

clc; clear variables; close all;

Nx = 1000;
x = linspace(-20e-6, 20e-6, Nx);
global l m;
l = 2; m = 2;

%% Calculate Target

LP = LP_Mode(x);
LP.F = LP.F/max(max(abs(LP.F)));

figure;
imagesc(abs(LP.F));
axis square;
title('Target');

%% Add phase tilt to target

x = linspace(0, 1, Nx);
[x_mesh, y_mesh] = meshgrid(x, x.');
LP.F = LP.F.*exp(1i*x_mesh*250*2*pi + 1i*y_mesh*250*2*pi);

%% Calculate f function used for encoding hologram

syms symf;

arra = 0:0.01:1;
f = zeros(size(arra));
for i = 1:length(arra)
    
    a = arra(i);
    S = vpasolve(besselj(1, symf) == 0.5*a, symf, [0 1.85]);
    f(i) = double(S);
    
end

figure;
plot(arra, f);
xlabel('Amplitude a');
ylabel('f(a)');

%% Calculate hologram

Holo = interp1(arra, f, abs(LP.F));
Holo = exp(1i*Holo.*sin(angle(LP.F)));

figure('Position', [400 400 400 400]);
imagesc(angle(Holo));
axis square;
xticks(''); yticks('');

%% Apply Fourier transform, filter, Fourier transform

invHolo = fftshift(fft2(fftshift(Holo)))/Nx;
invHoloFiltered = zeros(Nx, Nx);
invHoloFiltered(Nx/2 - Nx/8:Nx/2 + Nx/8,Nx/2 - Nx/8:Nx/2 + Nx/8) = invHolo(3*Nx/4-Nx/8:3*Nx/4+Nx/8,3*Nx/4-Nx/8:3*Nx/4+Nx/8);
Replay = fftshift(fft2(fftshift(invHoloFiltered)))/Nx;

figure('Position', [400 400 400 400]);
ComplexPlot(Replay);
axis square;
xticks(''); yticks('');












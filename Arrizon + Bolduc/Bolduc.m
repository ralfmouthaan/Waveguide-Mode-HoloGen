% Ralf Mouthaan
% University of Cambridge
% November 2020
% 
% Coding up of Bolduc's algorithm

clc; clear variables; close all;

Nx = 1000;
x = linspace(-20e-6, 20e-6, Nx);
global l m;
l = 2; m = 2;

%% Calculate Target

LP = LP_Mode(x);
LP.F = LP.F/max(max(abs(LP.F)));

figure;
imagesc(LP.F);
axis square;
title('Target');

%% Add phase tilt to target

x = linspace(0, 1, Nx);
[x_mesh, y_mesh] = meshgrid(x, x.');
LP.F = LP.F;

%% Calculate f function used for encoding hologram

syms symf;

arra = 0:0.01:1;
f = zeros(size(arra));
for i = 1:length(arra)
    
    a = arra(i);
    S = vpasolve(sin(symf)./symf == a, symf, [-pi 0]);
    f(i) = double(S);

    
end

figure;
plot(arra, f);
xlabel('Amplitude a');
ylabel('f(a)');

%% Calculate hologram

M = interp1(arra, f, abs(LP.F));
M = 1 + M/pi;
F = angle(LP.F) - pi*M;

Holo = exp(1i*M.*mod(F + x_mesh*250*2*pi + y_mesh*250*2*pi, 2*pi));

%% Apply Fourier transform, filter, Fourier transform

invHolo = fftshift(fft2(fftshift(Holo)))/Nx;
figure; imagesc(abs(invHolo));
invHoloFiltered = zeros(Nx, Nx);
invHoloFiltered(Nx/2 - Nx/8:Nx/2 + Nx/8,Nx/2 - Nx/8:Nx/2 + Nx/8) = invHolo(3*Nx/4-Nx/8:3*Nx/4+Nx/8,3*Nx/4-Nx/8:3*Nx/4+Nx/8);

figure; imagesc(abs(invHoloFiltered));
Replay = fftshift(fft2(fftshift(invHoloFiltered)))/Nx;

figure;
imagesc(angle(Replay));
axis square;
title('Replay');












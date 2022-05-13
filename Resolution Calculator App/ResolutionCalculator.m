% Ralf Mouthaan
% University of Cambridge
% February 2021
%
% Calculator to determine size and resolution of replay field in fibre
% facet plane. Ideally the overall size of the replay field would cover the
% entire fibre facet, but not much more.
% Setup: SLM --> f1 --> f2 --> f3 --> replay field (fibre facet)

clc; clear variables; close all;

Nx = 1000;
dx = 3.74e-6;
NMacro = 32;
lambda = 633e-9;
f3 = 4e-3;
f1 = 100e-3;
f2 = 75e-3;
dxMacro  = 3.74e-6*Nx/NMacro;

du = lambda*f3/NMacro/dxMacro*f1/f2;

fprintf('Wavelength: %0.1fnm\n', lambda*1e9);
fprintf('Width: %d Minipixels\n', Nx);
fprintf('Minipixel pitch: %0.2fum\n', dx*1e6);
fprintf('Width: %d Macropixels\n', NMacro);
fprintf('Macropixel pitch: %0.1fum\n', dxMacro*1e6);
fprintf('Lens 1: %0.1fmm\n', f1*1e3);
fprintf('Lens 2: %0.1fmm\n', f2*1e3);
fprintf('Lens 3: %0.1fmm\n', f3*1e3);
fprintf('dx in fibre facet plane: %0.2fum\n', du*1e6);
fprintf('Width of fibre facet plane: %0.1fum\n', NMacro*du*1e6);
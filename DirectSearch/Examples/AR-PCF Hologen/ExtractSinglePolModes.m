% Ralf Mouthaan
% University of Cambridge
% November 2020
%
% Script to extract + save single-polarisation modes. These are the ones to
% subsequently be used for alignment.

clc; clear variables; close all;

load('FD Result - LP Modes.mat', 'RetVal');
ModeInds = [1 8 10 11 21 3 4];
alpha = -20*pi/180;
x1 = linspace(-1, 1, size(RetVal.Ex{1}, 1));
x2 = linspace(-0.70, 0.70, 100);

for i = 1:length(ModeInds)
    
    Ind = ModeInds(i);
    
    F = real(RetVal.Ex{Ind});
    F = interp2(x1, x1.', F, cos(alpha)*x2 - sin(alpha)*x2.', sin(alpha)*x2 + cos(alpha)*x2.');
    
    figure; imagesc(F); axis square;
    
    writematrix(F, ['Erlangen AR-PCF - Mode ' num2str(Ind) '.csv']);
    
end
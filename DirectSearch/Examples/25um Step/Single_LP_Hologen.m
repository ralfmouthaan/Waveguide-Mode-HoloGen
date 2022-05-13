% Ralf Mouthaan
% University of Cambridge
% November 2019
%
% Script to generate hologram corresponding for a given LP mode.

function Single_LP_Hologen

    clc; clear variables; close all;
    global l m;
    
    %% User-Entered
    
    l = 2;
    m = 1;
    lambda = 633e-9;
    f1 = 15e-3;
    f2 = 10e-3;
    Nx = 10000;
    Magnification = 100/150;
    HoloDiameter = 1250; % in Pixels
    FibreDiameter = 25e-6;
    dx = 3.74e-6;
    SMF_NA = 0.14;
    
    %% Illumination Mask
    
    x = linspace(-Nx/2,Nx/2,Nx);
    [x_mesh, y_mesh] = meshgrid(x, x.');
    r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);
    Illum = ones(size(r_mesh));
    Illum(r_mesh > HoloDiameter/2) = 0;
    
    w0 = f1*SMF_NA;
    x = linspace(0, Nx-1, Nx);
    x = x - (Nx-1)/2;
    x = x*dx;
    [x_mesh, y_mesh] = meshgrid(x, x.');
    r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);
    Illum = Illum.*exp(-r_mesh.^2/w0^2);
    
    figure;
    imagesc(x*1e3, x.'*1e3, Illum);
    axis square
    title('Illumination');
    xlabel('mm'); ylabel('mm');
    
    % Distances are normalised to 1. Delete to avoid confusion down the
    % line.
    clear x r_mesh x_mesh y_mesh;
    
    % Correct dx for magnification
    dx = dx*Magnification;
    
    %% Target Replay Field

    du = lambda*f2/(Nx*dx);
    u = linspace(0, Nx-1, Nx);
    u = u - (Nx-1)/2;
    u = u*du;
    Mode = LP_Mode(u);
    Target = Mode.F;
    fprintf('Fibre facet is %d pixels across\n', round(FibreDiameter/du));
    
    figure;
    subplot(1,2,1);
    imagesc(u*1e6, u.'*1e6, abs(Target));
    axis square
    title('|Target|');
    xlabel('\mum'); ylabel('\mum');
    xlim([-35 35]);
    ylim([-35 35]);
    subplot(1,2,2);
    imagesc(u*1e6, u.'*1e6, angle(Target));
    axis square
    title('\angleTarget');
    xlabel('\mum'); ylabel('\mum');
    xlim([-35 35]);
    ylim([-35 35]);
    
    %% Generate Mask
    
    [u_mesh, v_mesh] = meshgrid(u, u.');
    r_mesh = sqrt(u_mesh.^2 + v_mesh.^2);
    Mask = ones(size(Target));
    Mask(r_mesh > 20e-6) = 0;
    Mask = logical(Mask);
    
    %% Generate Hologram
    
    pause(0.1);
    Holo = DirectSearchSymmetryMulti(Illum, Target, Mask);
    
    %% Plot hologram
    figure;
    imagesc(angle(Holo));
    title('Hologram');
    axis square;
    
    %% Plot Replay
    
    Replay = fftshift(fft2(fftshift(Holo)));
    
    figure;
    subplot(1,2,1);
    imagesc(u*1e6, u.'*1e6, abs(Replay));
    axis square
    title('|Replay|');
    xlabel('\mum'); ylabel('\mum');
    xlim([-35 35]);
    ylim([-35 35]);
    subplot(1,2,2);
    imagesc(u*1e6, u.'*1e6, angle(Replay));
    axis square
    title('\angleReplay');
    xlabel('\mum'); ylabel('\mum');
    xlim([-35 35]);
    ylim([-35 35]);
    
    %% Save result
    
    Holo = Holo(size(Holo,1)/2 - HoloDiameter/2:size(Holo,1)/2 + HoloDiameter/2-1,...
        size(Holo,2)/2 - HoloDiameter/2:size(Holo,2)/2 + HoloDiameter/2-1);
    writematrix(angle(Holo), ['Holo - LP' num2str(l) num2str(m) '.csv']);

end
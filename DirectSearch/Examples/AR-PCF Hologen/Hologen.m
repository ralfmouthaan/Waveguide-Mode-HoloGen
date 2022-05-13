% Ralf Mouthaan
% University of Cambridge
% November 2019
%
% Script to generate all LP mode holograms for a set of modes calculated
% using the FDFD code.

function Hologen

    clc; clear variables; close all;
    
    %% Load target
    
    load('FD Result - LP Modes.mat', 'RetVal');    
    RetVal.x = (-size(RetVal.n, 1)/2+1/2:size(RetVal.n,1)/2-1/2)*RetVal.dx;
    
    %% User-Entered
    
    identifier = 'Erlangen Simplified PCF';
    lambda = 633e-9;
    f1 = 15e-3;
    f2 = 20e-3;
    Nx = 10000;
    HoloDiameter = 1000;
    FibreDiameter = 50e-6;
    Magnification = 75/100;
    dx = 3.74e-6*Magnification;
    SMF_NA = 0.14;
    
    %% Diffraction field coordinates
    
    x = linspace(-Nx/2,Nx/2,Nx);
    [x_mesh, y_mesh] = meshgrid(x, x.');
    r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);
    Illum = ones(size(r_mesh));
    Illum(r_mesh > HoloDiameter/2) = 0;
    
    %% Diffraction field illumination mask
    
    w0 = f1*SMF_NA*Magnification;
    x = linspace(0, Nx-1, Nx);
    x = x - (Nx-1)/2;
    x = x*dx;
    [x_mesh, y_mesh] = meshgrid(x, x.');
    r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);
    Illum = Illum.*exp(-r_mesh.^2/w0^2);
    
    % Distances are normalised to 1. Delete to avoid confusion down the
    % line.
    clear x r_mesh x_mesh y_mesh;

    %% Replay field coordinates

    du = lambda*f2/(Nx*dx);
    u = linspace(0, Nx-1, Nx);
    u = u - (Nx-1)/2;
    u = u*du;
    FibreRadiusPx = round(FibreDiameter/du/2);
    fprintf('Fibre facet is %d pixels in radius\n', FibreRadiusPx);
    
    %% Replay field ROI mask
    
    [u_mesh, v_mesh] = meshgrid(u, u.');
    r_mesh = sqrt(u_mesh.^2 + v_mesh.^2);
    Mask = ones(size(u_mesh));
    Mask(r_mesh > FibreDiameter) = 0;
    Mask = logical(Mask);
    
    %% The Loop
    
    for i = 1:length(RetVal.Ex)
        
        if RetVal.HorizontalPol(i) == 1
            
            %% Calculate Target
        
            F = RetVal.Ex{i};
            F = real(F);
            F = interp2(RetVal.x, RetVal.x.', F, u, u.', 'linear', 0);

            %% Generate hologram
            
            fprintf('%d x\n', i);
            HoloStruct = DirectSearchSymmetryBinary(Illum, F, Mask);

            Holo = HoloStruct.Holo(Nx/2 - HoloDiameter/2:Nx/2 + HoloDiameter/2 - 1, ...
                Nx/2 - HoloDiameter/2:Nx/2 + HoloDiameter/2 - 1);
            HoloFilename = [identifier ' - Holo - LP ' num2str(i) 'x.csv'];
            writematrix(Holo, HoloFilename);
        
        end
        
        if RetVal.VerticalPol(i) == 1
            
            %% Calculate Target
        
            F = RetVal.Ey{i};
            F = real(F);
            F = interp2(RetVal.x, RetVal.x.', F, u, u.', 'nearest', 0);

            %% Generate hologram

            fprintf('%d y\n', i);
            HoloStruct = DirectSearchRotSymmetryBinary(Illum, F, Mask);

            Holo = HoloStruct.Holo(Nx/2 - HoloDiameter/2:Nx/2 + HoloDiameter/2 - 1, ...
                Nx/2 - HoloDiameter/2:Nx/2 + HoloDiameter/2 - 1);
            HoloFilename = [identifier ' - Holo - LP ' num2str(i) 'y.csv'];
            writematrix(Holo, HoloFilename);
        
        end
        
    end

        
end
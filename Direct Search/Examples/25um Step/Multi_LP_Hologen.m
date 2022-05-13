% Ralf Mouthaan
% University of Cambridge
% November 2019
%
% Script to generate all LP mode holograms for a given step index fibre

function Multi_LP_Hologen

    clc; clear variables; close all;
    global l m;
    
    %% User-Entered
    
    identifier = '25um Step';
    lambda = 633e-9;
    f1 = 15e-3;
    f2 = 10e-3;
    HoloDiameter = 1000;
    FibreDiameter = 25e-6;
    Magnification = 75/100;
    dx = 3.74e-6*Magnification;
    Nx = 10000;
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
    Mask(r_mesh > 20e-6) = 0;
    Mask = logical(Mask);
    
    %% The Loop
    
    l = 0;
    m = 1;
    
    arrl = []; arrm = []; arrHoloFile = {}; arrModeFile = {}; arrBeta = []; arrEfficiency = [];
    
    while true

         %% Generate Mode

        ModeStruct = LP_Mode(u);
        
        %% Check if we need to continue
        if ModeStruct.bolGuided == false
            if m == 1
                break
            else
                m = 1;
                l = l + 1;
                continue
            end
        end
        
        disp(['LP ' num2str(l) ' ' num2str(m)]);
        pause(0.1);

        %% Mode, Replay field calculations - Non-rotated mode
        
        HoloStruct = DirectSearchSymmetryBinary(Illum, ModeStruct.F, Mask);
        
        Holo = HoloStruct.Holo(Nx/2 - HoloDiameter/2:Nx/2 + HoloDiameter/2 - 1, ...
            Nx/2 - HoloDiameter/2:Nx/2 + HoloDiameter/2 - 1);
        F = ModeStruct.F(Nx/2 - FibreRadiusPx:Nx/2+FibreRadiusPx - 1, ...
            Nx/2 - FibreRadiusPx:Nx/2+FibreRadiusPx - 1);

        %% Save hologram + mode - Non-rotated mode
        
        if l == 0
            HoloFilename = [identifier ' - Holo - LP ' num2str(l) num2str(m) '.csv'];
            ModeFilename = [identifier ' - Mode - LP ' num2str(l) num2str(m) '.csv'];
        else
            HoloFilename = [identifier ' - Holo - LP ' num2str(l) num2str(m) 'a.csv'];
            ModeFilename = [identifier ' - Mode - LP ' num2str(l) num2str(m) 'a.csv'];
        end
            
        writematrix(Holo, HoloFilename);
        writematrix(F, ModeFilename);

        arrl = [arrl l];
        arrm = [arrm m];
        arrHoloFile = [arrHoloFile HoloFilename];
        arrModeFile = [arrModeFile ModeFilename];
        arrBeta = [arrBeta ModeStruct.Beta];
        arrEfficiency = [arrEfficiency HoloStruct.globalc];
        
        %% Mode, Replay field calculations - rotated mode
        
        if l > 0
            
            disp(['LP ' num2str(l) ' ' num2str(m) 'b']);
            pause(0.1);
        
            HoloStruct = DirectSearchSymmetryBinary(Illum, ModeStruct.F_rotated, Mask);

            Holo = HoloStruct.Holo(Nx/2 - HoloDiameter/2:Nx/2 + HoloDiameter/2-1, ...
                Nx/2 - HoloDiameter/2:Nx/2 + HoloDiameter/2-1);
            F = ModeStruct.F_rotated(Nx/2 - FibreRadiusPx:Nx/2+FibreRadiusPx - 1, ...
                Nx/2 - FibreRadiusPx:Nx/2+FibreRadiusPx - 1);

            %% Save hologram + mode - Non-rotated mode

            HoloFilename = [identifier ' - Holo - LP ' num2str(l) num2str(m) 'b.csv'];
            ModeFilename = [identifier ' - Mode - LP ' num2str(l) num2str(m) 'b.csv'];

            writematrix(Holo, HoloFilename);
            writematrix(F, ModeFilename);

            arrl = [arrl l];
            arrm = [arrm m];
            arrHoloFile = [arrHoloFile HoloFilename];
            arrModeFile = [arrModeFile ModeFilename];
            arrBeta = [arrBeta ModeStruct.Beta];
            arrEfficiency = [arrEfficiency HoloStruct.globalc];
            
        end
        
        m = m + 1;

    end
    
    %% Write metadata
    
    % Sort
    [arrBeta, idxs] = sort(arrBeta, 'descend');
    arrl = arrl(idxs);
    arrm = arrm(idxs);
    arrHoloFile = arrHoloFile(idxs);
    arrModeFile = arrModeFile(idxs);
    
    % Print to text file
    fid = fopen([identifier ' - Metadata.txt'], 'w+');
    fprintf(fid, 'l\tm\tBeta\tHoloFile\tModeFile\tEfficiency\n');
    for i = 1:length(arrl)
        fprintf(fid, '%d\t%d\t%f\t%s\t%s\t%f\n', arrl(i), arrm(i), arrBeta(i), arrHoloFile{i}, arrModeFile{i}, arrEfficiency(i));
    end
    fclose(fid);
        
end
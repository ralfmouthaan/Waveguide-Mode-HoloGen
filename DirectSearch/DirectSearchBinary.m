% Ralf Mouthaan
% University of Cambridge
% November 2019
%
% Direct search algorithm based on Joel Carpenter's Optical Engineering
% paper on GPU-accelerated simulated annealing.
% This is specifically designed to optimise for fibre modes.

function [Holo] = DirectSearch(Holo, Target, Mask)
    % HOLOANNEAL Perform annealing
    % Author: RPM
    % An illumination field is modulated by a hologram to give a desired 
    % replay field response. Direct search is used for optimisation 
    % of the hologram. It is assumed the hologram is at the near focal 
    % point of a lens, and the image at the far focal point, with the two 
    % related by a simple Fourier transform.
    % Inputs:
    %   Illumination : SLM illumination field
    %   Target : Target field
    %   Mask : Binary mask with 1s corresponding to region of interest.
    % Outputs:
    %   OptimalHolo : Optimal hologram
    
    %% Error Checks
    
    if size(Holo, 1) ~= size(Holo, 2)
        error('HoloAnneal: hologram must be square');
    end
    
    if size(Holo) ~= size(Target)
        error('HoloAnneal: illumination and target must be same size')
    end
    
    if nargin > 3
        error('Too many input arguments')
    end
    
    %% Pre-Processing
    
    Nx = size(Holo, 1);
    
    % fftshifts
    Target = fftshift(Target);
    Mask = fftshift(Mask);
    Holo = fftshift(Holo);
    
    % Initialise a few parameters
    dispstr = '';
    Replay = fft2(Holo);
    x = linspace(0, Nx-1, Nx);
    [x_mesh, y_mesh] = meshgrid(x, x.');
    
    % Ensure power in target is half the power in the hologram. This
    % provides an extra degree of freedom.
    Target = Target ...
        / sqrt(0.25*sum(sum(abs(Target).^2))) * sqrt(sum(sum(abs(Holo).^2))) * ...
        sqrt(size(Target, 1)) * sqrt(size(Target, 2));
    
    % Apply mask to relevant arrays
    x_mesh = x_mesh(Mask);
    y_mesh = y_mesh(Mask);
    Replay = Replay(Mask);
    MaskedTarget = conj(Target(Mask));
    
    TargetPower = sum(sum(abs(Target).^2));
    MaskedTargetPower = sum(sum(abs(MaskedTarget).^2));
    
    %% Direct Search
    
    maxlocalc = 0;
    maxglobalc = 0;
    IterNo = 0;
    dispstr = '';
    gamma = 10;
    LastDisplayed = 0;
    
    while IterNo < 5e6 %|| (maxglobalc < 0.04)
        
        % Random pixel position + value
        m = randi(Nx);
        n = randi(Nx);
        OldPixel = Holo(m,n);
        if abs(OldPixel) < 0.01
            continue
        end
        NewPixel = -Holo(m,n);
        
        % Update replay
        NewReplay = Replay;
        DeltaReplay = (NewPixel - OldPixel)*exp(-2*pi*1i/Nx*((n-1)*x_mesh + (m-1)*y_mesh));
        NewReplay = NewReplay + DeltaReplay;
        
        % Calculate local overlap integral
        OverlapIntegral = abs(sum(sum(NewReplay.*MaskedTarget))).^2;
        MaskedReplayPower = sum(sum(abs(NewReplay).^2));
        localc = OverlapIntegral / MaskedReplayPower;
        localc = localc / MaskedTargetPower;
        globalc = MaskedReplayPower / TargetPower;
        
        % Determine if result is better
        if (localc^gamma*globalc > maxlocalc^gamma*maxglobalc)
                       
            maxlocalc = localc;
            maxglobalc = globalc;
            Holo(m,n) = NewPixel;
            Replay = NewReplay;
            
            if IterNo - LastDisplayed > 50
                fprintf(repmat('\b',1,length(dispstr)-1));
                dispstr = [num2str(IterNo) '; Local Overlap = ' num2str(maxlocalc) '; Global Overlap = ' num2str(maxglobalc) '\n'];
                fprintf(dispstr);
                LastDisplayed = IterNo;
            end
            
        end
        
        IterNo = IterNo + 1;
       
    end
    
    Holo = ifftshift(Holo);
    
end
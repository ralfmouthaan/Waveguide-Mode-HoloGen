% Ralf Mouthaan
% University of Cambridge
% December 2019
%
% Modified Gerchberg-Saxton algorithm. In his Nature paper, Cizmar
% specifies that he uses a modified Gerchberg-Saxton paper to generate his
% holograms. Specifically, in the replay field he constrains both amplitude
% and phase but only in a small region.

function [Holo] = ModifiedGerchbergSaxton(Holo, Target, Mask)

    % Holo is the entire hologram, and includes the illumination.
    % ReplayMask is the region of interest of the replay field.
    % Target is the target replay field, and is the same size as the
    % hologram. Really, only the bit in the ROI is relevant and the rest is
    % junk.
    
    bolDebug = true;

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
    
    %% Pre-processing
    
    Illum = abs(Holo);
    
    Nx = size(Holo, 1);
    
    % fftshifts
    Target = fftshift(Target);
    Mask = fftshift(Mask);
    Holo = fftshift(Holo);
    Illum = fftshift(Illum);
    
    % Scale target
    Target = Target ...
        / sqrt(0.1*sum(sum(abs(Target).^2))) * sqrt(sum(sum(abs(Holo).^2))) * ...
        sqrt(size(Target, 1)) * sqrt(size(Target, 2));
    Replay = Target;
    
    %% Main Loop
    
    if bolDebug == true
        figure;
    end
    
    while true
        
        Replay(Mask) = Target(Mask);
        
        if bolDebug == true
            subplot(2,2,3); imagesc(fftshift(abs(Replay))); title('|Replay|'); axis square;
            subplot(2,2,4); imagesc(fftshift(angle(Replay))); title('\angleReplay'); axis square;
            pause(0.1)
        end
        
        Holo = ifft2(Replay);
        
        if bolDebug == true
            subplot(2,2,1); imagesc(fftshift(abs(Holo))); title('|Holo|'); axis square;
            subplot(2,2,2); imagesc(fftshift(angle(Holo))); title('\angleHolo'); axis square;
            pause(0.1)
        end
        
        Holo = Illum.*exp(1i*angle(Holo));
        
        if bolDebug == true
            subplot(2,2,1); imagesc(fftshift(abs(Holo))); title('|Holo|'); axis square;
            subplot(2,2,2); imagesc(fftshift(angle(Holo))); title('\angleHolo'); axis square;
            pause(0.1)
        end
        
        Replay = fft2(Holo);
        
        if bolDebug == true
            subplot(2,2,3); imagesc(fftshift(abs(Replay))); title('|Replay|'); axis square;
            subplot(2,2,4); imagesc(fftshift(angle(Replay))); title('\angleReplay'); axis square;
            pause(0.1)
        end
                
    end

end
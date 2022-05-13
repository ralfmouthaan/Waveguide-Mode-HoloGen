% Ralf Mouthaan
% University of Cambridge
% November 2019
%
% Direct search algorithm based on Joel Carpenter's Optical Engineering
% paper on GPU-accelerated simulated annealing.
% This is specifically designed to optimise for fibre modes.

function [RetVal] = DirectSearchSymmetryBinary(Illum, Target, ReplayMask)
    % HOLOANNEAL Perform annealing
    % Author: RPM
    % An illumination field is modulated by a hologram to give a desired 
    % replay field response. Direct search is used for optimisation 
    % of the hologram. It is assumed the hologram is at the near focal 
    % point of a lens, and the image at the far focal point, with the two 
    % related by a simple Fourier transform.
    % Inputs:
    %   Holo : Hologram. Amplitudes must reflect illumination
    %   Target : Target field
    %   Mask : Binary mask with 1s corresponding to region of interest.
    % Outputs:
    %   OptimalHolo : Optimal hologram
    
    bolOutputToFile = false;
    
    %% Error Checks
    
    if size(Illum, 1) ~= size(Illum, 2)
        error('HoloAnneal: hologram must be square');
    end
    
    if size(Illum) ~= size(Target)
        error('HoloAnneal: illumination and target must be same size')
    end
    
    if size(Illum) ~= size(ReplayMask)
        error('HoloAnneal: illumination and replay mask must be same size')
    end
    
    if nargin > 3
        error('Too many input arguments')
    end
    
    %% Pre-Processing
    
    Nx = size(Illum, 1);
    
    % Determine UD axis of symmetry
    TargetTop = Target(1:Nx/2,:);
    TargetBottom = Target((Nx/2+1):Nx,:);
    TargetBottom = flipud(TargetBottom);
    if max(max(abs(TargetTop - TargetBottom))) < max(max(abs(TargetTop + TargetBottom)))
        fprintf('UD Symmetry Score: %0.5f\n', max(max(abs(TargetTop - TargetBottom)))/max(max(abs(TargetTop))))
        UpDownSymmetry = 1;
    else
        UpDownSymmetry = -1;
        fprintf('UD Symmetry Score: %0.5f\n', max(max(abs(TargetTop + TargetBottom)))/max(max(abs(TargetTop))))
    end
    clear TargetTop TargetBottom
    
    % Determine LR axis of symmetry
    TargetLeft = Target(:,1:Nx/2);
    TargetRight = Target(:,(1 + Nx/2):Nx);
    TargetRight = fliplr(TargetRight);
    if max(max(abs(TargetLeft - TargetRight))) < max(max(abs(TargetLeft + TargetRight)))
        fprintf('LR Symmetry Score: %0.5f\n', max(max(abs(TargetLeft - TargetRight)))/max(max(abs(TargetLeft))))
        LeftRightSymmetry = 1;
    else
        fprintf('LR Symmetry Score: %0.5f\n', max(max(abs(TargetLeft + TargetRight)))/max(max(abs(TargetLeft))))
        LeftRightSymmetry = -1;
    end
    
    % fftshifts
    Target = fftshift(Target);
    ReplayMask = fftshift(ReplayMask);
    Illum = fftshift(Illum);
    
    % Coordinate calculations
    x = linspace(0, Nx-1, Nx);
    [x_mesh, y_mesh] = meshgrid(x, x.');
    
    % Initial hologram value.
    % Have tried starting the hologram as the ifft of the replay. This
    % gives a large amount of power coupling in, but a very small overlap
    % integral that doesn't improve much. Seems to be some kind of local
    % minimum it gets trapped in.
    Holo = Illum.*(randi(2, size(Illum))*2 - 3);
    Replay = fft2(Holo);
    
    % Ensure power in target is half the power in the hologram. This
    % provides an extra degree of freedom.
    Target = Target ...
        / sqrt(sum(sum(abs(Target).^2))) * sqrt(sum(sum(abs(Holo).^2))) * ...
        sqrt(size(Target, 1)) * sqrt(size(Target, 2));
    
    % Take NW quadrant in all cases.
    Target = Target(1:Nx/2,1:Nx/2);
    ReplayMask = ReplayMask(1:Nx/2,1:Nx/2);
    Holo = Holo(1:Nx/2,1:Nx/2);
    x_mesh = x_mesh(1:Nx/2,1:Nx/2);
    y_mesh = y_mesh(1:Nx/2,1:Nx/2);
    
    % Apply mask to relevant arrays
    x_mesh = x_mesh(ReplayMask);
    y_mesh = y_mesh(ReplayMask);
    Replay = Replay(ReplayMask);
    MaskedTarget = conj(Target(ReplayMask));
    
    TargetPower = sum(sum(abs(Target).^2));
    MaskedTargetPower = sum(sum(abs(MaskedTarget).^2));
    
    %% Direct Search
    
    maxlocalc = 0;
    maxglobalc = 0;
    IterNo = 0;
    dispstr = '';
    gamma = 10;
    LastDisplayed = 0;
    IllumEdge = Nx/2;
    
    if bolOutputToFile == true
        fid_localc = fopen('Local c.txt', 'w+');
        fid_globalc = fopen('Global c.txt', 'w+');
    end
    
    while IterNo < 1.5e6 %|| (maxglobalc < 0.04)
        
        % Random pixel position + value
        %m = round((IllumEdge-1)*rand*cos(rand*pi/2))+1;
        %n = round((IllumEdge-1)*rand*cos(rand*pi/2))+1;
        m = randi(IllumEdge);
        n = randi(IllumEdge);
        OldPixel = Holo(m,n);
        if abs(OldPixel) < 0.001 && sqrt(m^2 + n^2) < IllumEdge
            IllumEdge = round(sqrt(m^2 + n^2));
            continue
        end
        NewPixel = -Holo(m,n);
        
        % Update replay
        NewReplay = Replay;
        DeltaReplay = (NewPixel - OldPixel)*...
            exp(-2*pi*1i/Nx*((n-1)*x_mesh + (m-1)*y_mesh));
        DeltaReplay = DeltaReplay + ...
            LeftRightSymmetry*(NewPixel - OldPixel)*...
            exp(-2*pi*1i/Nx*(-(n-1)*x_mesh + (m-1)*y_mesh));
        DeltaReplay = DeltaReplay + ...
            UpDownSymmetry*(NewPixel - OldPixel)*...
            exp(-2*pi*1i/Nx*((n-1)*x_mesh - (m-1)*y_mesh));
        DeltaReplay = DeltaReplay + ...
            LeftRightSymmetry*UpDownSymmetry*(NewPixel - OldPixel)*...
            exp(-2*pi*1i/Nx*(-(n-1)*x_mesh - (m-1)*y_mesh));
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
            
            if IterNo - LastDisplayed > 1000
                
                fprintf(repmat('\b',1,length(dispstr)-1));
                dispstr = [num2str(IterNo) '; Local Overlap = ' num2str(maxlocalc) '; Global Overlap = ' num2str(maxglobalc) '\n'];
                fprintf(dispstr);
                LastDisplayed = IterNo;
                
                if bolOutputToFile == true
                    fprintf(fid_localc, '%d\t%f\n', IterNo, maxlocalc);
                    fprintf(fid_globalc, '%d\t%f\n', IterNo, maxglobalc);
                end
                
            end
            
        end
        
        IterNo = IterNo + 1;
       
    end
    
    if bolOutputToFile == true
        fclose(fid_localc);
        fclose(fid_globalc);
    end
    
    Holo = [Holo LeftRightSymmetry*fliplr(Holo)];
    Holo = [Holo ; UpDownSymmetry*flipud(Holo)];
    Holo = fftshift(Holo);
    Holo(Holo>=0) = 1;
    Holo(Holo<0) = -1;
    
    RetVal.Holo = Holo;
    RetVal.globalc = maxglobalc;
    RetVal.localc = maxlocalc;
    
end
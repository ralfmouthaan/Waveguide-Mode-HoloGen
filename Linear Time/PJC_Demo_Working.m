% Ralf Mouthaan
% University of Cambridge
% May 2020
%
% Linear-time hologram generator for 25um fibre.

clc; clear variables; close all;
    
%% User-Entered Variables

identifier = '25um Step';

f1 = 15e-3;
f2 = 75;
f3 = 100;
f4 = 10e-3;
Nx = 2500;
lambda = 633e-9;
gamma = 1;

HoloDiameter = 1000;
structDF.dx = 3.74e-6;
SMF_NA = 0.14;

structTargetLP.n_co = 1.457; % Silica
structTargetLP.n_cl = 1.4403;
structTargetLP.l=3; 
structTargetLP.m=3;
structTargetLP.lambda = lambda;
structTargetLP.WaveguideDiameter = 25e-6;
    
%% Diffraction field

structDF.x = linspace(-Nx/2,Nx/2,Nx);
[structDF.x_mesh, structDF.y_mesh] = meshgrid(structDF.x, structDF.x.');
structDF.r_mesh = sqrt(structDF.x_mesh.^2 + structDF.y_mesh.^2);
structDF.Mask = ones(size(structDF.r_mesh));
structDF.Mask(structDF.r_mesh > HoloDiameter/2) = 0;
structDF.x = structDF.x*structDF.dx;
structDF.x_mesh = structDF.x_mesh*structDF.dx;
structDF.y_mesh = structDF.y_mesh*structDF.dx;
structDF.r_mesh = structDF.r_mesh*structDF.dx;
structDF.Illumination = structDF.Mask.*exp(-structDF.r_mesh.^2/(f1*SMF_NA)^2);

%% Replay field coordinates

structRF.dx = lambda*f4/(Nx*structDF.dx*f2/f3);
structRF.x = linspace(-Nx/2,Nx/2,Nx);
structRF.x = structRF.x*structRF.dx;
[structRF.x_mesh, structRF.y_mesh] = meshgrid(structRF.x, structRF.x.');
structRF.r_mesh = sqrt(structRF.x_mesh.^2 + structRF.y_mesh.^2);
structRF.Mask = ones(size(structRF.x_mesh));
structRF.Mask(structRF.r_mesh > structTargetLP.WaveguideDiameter) = 0;
structRF.Mask = logical(structRF.Mask);

%%

% Note, we use C1 here to represent the inside of the top term of the
% and C2 to represent the bottom
% coupling equation. i.e. C=|C1|^2/C2 and

structTargetLP.x = structRF.x;
structTargetLP = LP_Mode(structTargetLP);

targetT = structTargetLP.F;
targetG = ModeFT(structTargetLP);

H=exp(1i.*angle(targetG));
initialH=H;
initialR=ft2(initialH);
targetC1=sum(sum(targetG.*conj(H)));
targetC2=(sum(sum(abs(targetG).^2))*sum(sum(abs(H).^2)));
targetC=abs(targetC1)^2/targetC2;

M = structRF.Mask;
L = ift2(M);

otherC1s=zeros(4,4);
otherC2s=zeros(4,4);
otherGs=zeros(4,4,Nx,Nx);
for lOther=1:4
    for mOther=1:4
        if (lOther~=structTargetLP.l) || (mOther~=structTargetLP.m)
            structCurrLP = structTargetLP;
            structCurrLP.l = lOther;
            structCurrLP.m = mOther;
            structCurrLP = LP_Mode(structCurrLP);
            T = structCurrLP.F;
            G = ModeFT(structCurrLP);
            otherGs(lOther,mOther,:,:) = G;
            otherC1s(lOther,mOther)=sum(sum(G.*conj(H)));
            otherC2s(lOther,mOther)=(sum(sum(abs(G).^2))*sum(sum(abs(H).^2)));
        end
    end
end
otherCs=abs(otherC1s).^2./otherC2s;
otherCs(isnan(otherCs))=0;
otherC=sum(sum(abs(otherCs).^2));

% Note, this does not include symmetry in pixel changes
% Note, this does not include symmetry in mode checks
% Note, an alternative technique to the gamma function will make this
% linear and therefore guaranteeable by merely iterating over x and y
numItrs=10^7;
for itr=1:numItrs
    
    if mod(itr,1000) == 0
        fprintf("Iteration %8d of %8d; %5.2f%% complete; targetC = %12.10f; otherC = %12.10f; \n",itr,numItrs,itr/numItrs*100,targetC,otherC);
    end
    x = randi([0 Nx-1],1,1);
    y = randi([0 Nx-1],1,1);
    new = H(x+1,y+1).*exp(1i.*pi); % binary case
    
    newTargetC1=targetC1+targetG(x+1,y+1).*(conj(new)-conj(H(x+1,y+1)));
    newTargetC=abs(newTargetC1)^2/targetC2;
    newOtherC1s=otherC1s+otherGs(:,:,x+1,y+1).*(conj(new)-conj(H(x+1,y+1)));
    newOtherCs=abs(newOtherC1s).^2./otherC2s;
    newOtherCs(isnan(newOtherCs))=0;
    newOtherC=sum(sum(abs(newOtherCs).^2)); % newOtherC=sum(sum(abs(newOtherCs).^2)); may be better
    
    deltaTargetC=newTargetC-targetC;
    deltaOtherC=newOtherC-otherC;
    
    if (newTargetC/newOtherC^gamma > targetC/otherC^gamma)
        H(x+1,y+1) = new;
        targetC1=newTargetC1;
        targetC=newTargetC;
        otherC1s=newOtherC1s;
        otherCs=newOtherCs;
        otherC=newOtherC;
    end
end
R=ft2(H);

%% Overlap Integrals

fprintf('Local Overlap %2.2f\n', OverlapIntegral(M.*R, M.*targetT));
fprintf('Global Overlap Integral %2.2f\n', OverlapIntegral(R, targetT));

%% Plotting

subplot(2,3,1);
imagesc(structRF.x, structRF.x, abs(targetT));
axis square;
xlim([-50e-6 50e-6]);
ylim([-50e-6 50e-6]);
title('targetMode.T');

subplot(2,3,2);
imagesc(abs(targetG));
axis square;
title('targetMode.G');

subplot(2,3,3);
imagesc(angle(initialH));
axis square;
title('initialH');

subplot(2,3,4);
imagesc(structRF.x, structRF.x, abs(R));
axis square;
xlim([-50e-6 50e-6]);
ylim([-50e-6 50e-6]);
title('R');

subplot(2,3,5);
imagesc(structRF.x, structRF.x, abs(initialR));
axis square;
xlim([-50e-6 50e-6]);
ylim([-50e-6 50e-6]);
title('initialR');

subplot(2,3,6);
imagesc(angle(H));
axis square;
title('H');


function G = ModeFT(structMode)

    G = ift2(structMode.F);
    
    return;

    % Calculates the FT of the mode using an analytical approach:

    T = structMode.F;
    r = structMode.r_mesh;
    a = structMode.WaveguideRadius;
    U = structMode.U;
    W = structMode.W;
    l = structMode.l;
    
    G= 1i.^l.*a.^2.*(a.*r.*besselj(l,U).*besselj(l-1,r.*a)-U.*besselj(l-1,U).*besselj(l,r.*a))./(4*pi*(U.^2-r.^2.*a.^2)*besselj(l,U))+...
       1i.^l.*a.^2.*(a.*r.*besselk(l,W).*besselj(l+1,r.*a)-W.*besselk(l+1,W).*besselj(l,r.*a))./(4*pi*(W.^2+r.^2.*a.^2)*besselk(l,W));
    
    G=G * sqrt(sum(sum(abs(T).^2))/sum(sum(abs(G).^2)));

end



function B = ft2(A)
    N=sqrt(size(A,1)*size(A,2));
    B=fftshift(fft2(fftshift(A)))/N;
end

function B = ift2(A)
    N=sqrt(size(A,1)*size(A,2));
    B=fftshift(ifft2(fftshift(A)))*N;
end


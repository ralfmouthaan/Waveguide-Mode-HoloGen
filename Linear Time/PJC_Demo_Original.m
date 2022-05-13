
% close all;
clear;
clf;

global a lambda n_co n_cl scale;
diameter = 50e-6;
a = diameter/2;
lambda = 633e-9;
n_co = 1.457; % Silica
n_cl = 1.4403;
scale=15; % fudge factor for ease of manipulation
gamma=10;
Nx=512; Ny=512;
l=3; m=3;

% Note, we use C1 here to represent the inside of the top term of the
% and C2 to represent the bottom
% coupling equation. i.e. C=|C1|^2/C2 and

[targetT, targetG] = getMode2D(l,m,Nx,Ny,1);
H=1.*exp(1i.*angle(targetG));
% H=1.*exp(1i.*rand(Nx,Ny)*2*pi);
initialH=H;
initialR=ft2(initialH);
targetC1=sum(sum(targetG.*conj(H)));
targetC2=(sum(sum(abs(targetG).^2))*sum(sum(abs(H).^2)));
targetC=abs(targetC1)^2/targetC2;

[M,L] = getMask(Nx,Ny);

otherTsum=targetT*0;
otherGsum=targetG*0;
otherTs=zeros(4,4,Nx,Ny);
otherGs=zeros(4,4,Nx,Ny);
otherC1s=zeros(4,4);
otherC2s=zeros(4,4);
for lOther=1:4
    for mOther=1:4
        if (lOther~=l) || (mOther~=m)
            [T,G]=getMode2D(lOther,mOther,Nx,Ny,1);
            otherTs(lOther,mOther,:,:)=T;
            otherGs(lOther,mOther,:,:)=G;
            otherTsum=otherTsum+T;
            otherGsum=otherGsum+G;
            otherC1s(lOther,mOther)=sum(sum(G.*conj(H)));
            otherC2s(lOther,mOther)=(sum(sum(abs(G).^2))*sum(sum(abs(H).^2)));
        end
    end
end
otherCs=abs(otherC1s).^2./otherC2s;
otherCs(isnan(otherCs))=0;
otherC=sum(sum(abs(otherCs).^2)); % otherC=sum(sum(abs(otherCs).^2)); may be better

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
    y = randi([0 Ny-1],1,1);
    %         new = 1.*exp(1i.*2.*pi.*rand(1,1)); % non-binary case
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


% fig=figure('Color',[1 1 1]);
% set(gcf, 'Position',  [200, 200, 800, 800]);
hold on

subplot(2,3,1);
hold on; axis square;
imshow(DomainColoredImage2(targetT));
title('targetMode.T');

subplot(2,3,2);
hold on; axis square;
imshow(DomainColoredImage2(targetG));
title('targetMode.G');

subplot(2,3,3);
hold on; axis square;
imshow(DomainColoredImage2(initialH));
title('initialH');

subplot(2,3,4);
hold on; axis square;
imshow(DomainColoredImage2(R));
title('R');

subplot(2,3,5);
hold on; axis square;
imshow(DomainColoredImage2(initialR));
title('initialR');

subplot(2,3,6);
hold on; axis square;
imshow(DomainColoredImage2(H));
title('H');

% subplot(2,3,4);
% hold on; axis square;
% imshow(DomainColoredImage2(otherModes.T));
% title('otherModes.T');
% 
% subplot(2,3,5);
% hold on; axis square;
% imshow(DomainColoredImage2(otherModes.G));
% title('otherModes.G');

colormap(gca, 'hsv');
hL = colorbar('Ticks',[0.0,0.25,0.5,0.75,1.0],...
    'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'});
newPosition = [0.915,0.085,0.02,0.875];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);

function [T,G] = getMode2D(l,m,Nx,Ny,rot)

    global a lambda n_co n_cl scale;
    [U, V, W, beta] = RPM_Step_Index_Solver(a, lambda, n_co, n_cl, l, m);
    NA = sqrt(n_co^2 - n_cl^2);
    ratioX = max(Nx/Ny,1);
    ratioY = max(Ny/Nx,1);
    x = linspace(-scale*a*ratioX,scale*a*ratioX,Nx+1); 
    y = linspace(-scale*a*ratioY,scale*a*ratioY,Ny+1);
    x=x(1:Nx); y=y(1:Ny); % Note this addition
    [x_mesh, y_mesh] = meshgrid(y, x.');
    
    r = sqrt(x_mesh.^2 + y_mesh.^2)/a;
    theta = atan2(y_mesh, x_mesh);
    T = besselj(l, U.*r)./besselj(l, U);
    T(r >  1) = besselk(l, W.*r(r> 1))./besselk(l, W);
    if rot
        T = T.*cos(l*theta);
    end
    
    r=r/scale/scale*sqrt(Nx*Ny)/a*pi/2; % scale factors
    G= 1i.^l.*a.^2.*(a.*r.*besselj(l,U).*besselj(l-1,r.*a)-U.*besselj(l-1,U).*besselj(l,r.*a))./(4*pi*(U.^2-r.^2.*a.^2)*besselj(l,U))+...
       1i.^l.*a.^2.*(a.*r.*besselk(l,W).*besselj(l+1,r.*a)-W.*besselk(l+1,W).*besselj(l,r.*a))./(4*pi*(W.^2+r.^2.*a.^2)*besselk(l,W));
    G=G*sqrt(Nx*Ny)/a/a*pi/4; % scale factors
    if rot
        G=G.*cos(l*theta);
    end
    
    normaliseFactor = sqrt(mean(mean(abs(T).^2)));
    T=T/normaliseFactor;
    G=G/normaliseFactor;
    
end

function [t,g] = getMode1D(l,m,N)
    global a lambda n_co n_cl scale;
    [U, V, W, beta] = RPM_Step_Index_Solver(a, lambda, n_co, n_cl, l, m);
    NA = sqrt(n_co^2 - n_cl^2);
    
    r = linspace(0,scale*a,N/2)/a;
    t = besselj(l, U.*r)./besselj(l, U);
    t(r >  1) = besselk(l, W.*r(r> 1))./besselk(l, W);
    
    r=r/scale/scale*N/a*pi/2; % scale factors
    g= 1i.^l.*a.^2.*(a.*r.*besselj(l,U).*besselj(l-1,r.*a)-U.*besselj(l-1,U).*besselj(l,r.*a))./(4*pi*(U.^2-r.^2.*a.^2)*besselj(l,U))+...
       1i.^l.*a.^2.*(a.*r.*besselk(l,W).*besselj(l+1,r.*a)-W.*besselk(l+1,W).*besselj(l,r.*a))./(4*pi*(W.^2+r.^2.*a.^2)*besselk(l,W));
    g=g*N/a/a*pi/4; % scale factors
end

function [M,L] = getMask(Nx,Ny)
    global a scale;
    ratioX = max(Nx/Ny,1);
    ratioY = max(Ny/Nx,1);
    x = linspace(-scale*a*ratioX,scale*a*ratioX,Nx+1);
    y = linspace(-scale*a*ratioY,scale*a*ratioY,Ny+1);
    x=x(1:Nx); y=y(1:Ny); % Note this addition
    [x_mesh, y_mesh] = meshgrid(y, x.');

    r = sqrt(x_mesh.^2 + y_mesh.^2)/a;
    M=zeros(Nx,Ny);
    M(r<1.2)=1;
    M=normalise(M);
    L=ift2(M);
end

function C = coupling(A,B)
    C=abs(sum(sum(A.*conj(B))))^2/(sum(sum(abs(A).^2))*sum(sum(abs(B).^2)));
end

function B = ft2(A)
    N=sqrt(size(A,1)*size(A,2));
    B=fftshift(fft2(fftshift(A)))/N;
end

function B = ift2(A)
    N=sqrt(size(A,1)*size(A,2));
    B=fftshift(ifft2(fftshift(A)))*N;
end

function B = normalise(A)
    B=A/sqrt(mean(mean(abs(A).^2)));
end

function C = cconv2(A,B)
    C=ift2(ft2(A).*ft2(B));
end

function C = cconv2ALT(A,B)
    Nx=size(A,1);
    Ny=size(A,2);
    N=sqrt(size(A,1)*size(A,2));
    C=zeros(Nx,Ny);
    for xprime=0:Nx-1
        for yprime=0:Ny-1
            for a=0:Nx-1
                for b=0:Ny-1
                    x1=mod(xprime+Nx/2,Nx);
                    y1=mod(yprime+Ny/2,Ny);
                    C(x1+1,y1+1)=C(x1+1,y1+1)+A(a+1,b+1)*B(mod(xprime-a,Nx)+1,mod(yprime-b,Ny)+1)/N;
                end
            end
        end
    end
end

function B = dagger(A)
    Nx=size(A,1);
    Ny=size(A,2);
    B=A;
    for x=0:Nx-1
        for y=0:Ny-1
            B(x+1,y+1)=conj(A(mod(-x,Nx)+1,mod(-y,Ny)+1));
        end
    end
end

function B = altDagger(A)
    B=ift2(conj(ft2(A)));
end
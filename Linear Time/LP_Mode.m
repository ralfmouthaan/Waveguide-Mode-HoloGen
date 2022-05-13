% Ralf Mouthaan
% University of Cambridge
% August 2019
%
% Code to generate field distribution for a given LP mode. To subsequently
% be used in hologram generation algorithm.
%
% l and m parameters have to be declared globally. The spatial axis x is passed in.
% Other parameters, such as wavelength, diameter, refractive indices, are all set
% (i.e. hardcoded) in this function.
%
% Code is based on original analytical solver code, which in turn is based
% on
%   * Snyder & Love Ch. 13 + Ch. 14
%   * Gloge, Applied Optics, 1971.

function [Mode] = LP_Mode(Mode)

%% User-Defined Parameters + Derived Parameters

l = Mode.l;
m = Mode.m;
x = Mode.x;
diameter = Mode.WaveguideDiameter;
a = diameter/2;
lambda = Mode.lambda;
n_co = Mode.n_co;
n_cl = Mode.n_cl;
k = 2*pi/lambda;
NA = sqrt(n_co^2 - n_cl^2); 
V = a*k*NA;
Delta = (1 - n_cl^2/n_co^2)/2;
[x_mesh, y_mesh] = meshgrid(x, x.');
normradius = sqrt(x_mesh.^2 + y_mesh.^2)/a;
r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);
theta = atan2(y_mesh, x_mesh);

Mode.r_mesh = r_mesh;
Mode.WaveguideDiameter = diameter;
Mode.WaveguideRadius = Mode.WaveguideDiameter/2;
Mode.Wavelength = lambda;
Mode.n_co = n_co;
Mode.n_cl = n_cl;
Mode.k = k;
Mode.WaveguideNA = NA;
Mode.IndexContrast = Delta;
Mode.l = l;
Mode.m = m;
Mode.x = x;
Mode.r = x(x>=0);
Mode.V = V;

%% Error Checks

if n_co < n_cl
    disp('WARNING: n_co must be > n_cl');
    return
end

if l < 0
    disp('WARNING: l must be >= 0');
    return
end

if m < 1
    disp('WARNING: m must be >= 1');
    return
end

if Delta > 1
    disp('WARNING: Waveguide is not weakly guiding');
    return
end

%% Eigenvalue equation

% The eigenvalue equation corresponds to the boundary matching conditions
% between the core and the cladding. The values U that satisfy the
% eigenvalue equation correspond to propagating modes.

% The eigenvalue equation f is evaluated at a range of U-values between 0
% and V.
arrU = 0:0.001:V;
funEigEq = @(x) EigenvalueEquation(x);
funEigEq(Mode);

% The findpeaks() function is used on arrEigEq to find the approximate
% position of the roots
[~,idxs] = findpeaks(-abs(EigenvalueEquation(arrU)));

if m > length(idxs)
    Guided = false;
    Mode.bolGuided = Guided;
else

    % fzero() is used to get a better approximation for the root of interest
    %U = fzero(funEigEq, arrU(idxs(m)));
    U = arrU(idxs(m));

    % Associated values are then calculated.
    W = sqrt(V^2 - U^2);
    beta = V/a/(2*Delta)^(1/2)*(1-2*Delta*U^2/V^2)^(1/2);

    if beta < k*n_cl || beta > k*n_co
        Guided = false;
    else
        Guided = true;
    end
    
    % Set values in returned structure
    Mode.U = U;
    Mode.W = W;
    Mode.bolGuided = Guided;
    Mode.Beta = beta;

end

%% Generate Fields

if Guided == true

    F = zeros(size(normradius));
    F(normradius <= 1) = besselj(l, U.*normradius(normradius<=1))/besselj(l, U);
    F(normradius >  1) = besselk(l, W.*normradius(normradius> 1))/besselk(l, W);

    F = F.*cos(l*theta);
    Mode.F = F;
    
    if l ~= 0
        
        % Rotated mode
        theta = theta + pi/2/l;
        
        F = zeros(size(normradius));
        F(normradius <= 1) = besselj(l, U.*normradius(normradius<=1))/besselj(l, U);
        F(normradius >  1) = besselk(l, W.*normradius(normradius> 1))/besselk(l, W);

        F = F.*cos(l*theta);
        Mode.F_rotated = F;
        
    end
    
end

end
function RetVal = EigenvalueEquation(U)

    % Eigenvalue equation as given in Snyder & Love. Is equivalent to the
    % Eigenvalue equation given in Gloge, who just used a different
    % recursive bessel function relationship for the derivation.
    %
    % Can pass in a fibre mode structure, in which case l and V are
    % updated. Else, we pass in U in which case the result of the
    % Eigenvalue equation is returned. This is done in this way so that we
    % can update l and V and still use this function in a optimiser without
    % having to resort to global variables.
    
    persistent l V;
    
    if isstruct(U)
        l = U.l;
        V = U.V;
        return;
    end

    W = sqrt(V^2 - U.^2);
    RetVal = U.*besselj(l+1,U)./besselj(l,U) - ...
        W.*besselk(l+1,W)./besselk(l,W);

end
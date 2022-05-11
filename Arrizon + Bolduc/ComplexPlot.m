% Ralf Mouthaan
% University of Cambridge
% November 2018
% September 2019: Turned into a function
%
% Script to plot complex data on a single graph

function ComplexPlot(Raw)

    maxmag = max(max(abs(Raw)));

    Real = real(Raw);
    Imag = imag(Raw);

    R = max(0, Real)/maxmag;
    G = (-min(0, Real) + max(0, Imag))/sqrt(2)/maxmag;
    B = (-min(0, Real) - min(0, Imag))/sqrt(2)/maxmag;

    M = cat(3, R, G, B);
    image(M);
    axis square;
    xticks('');
    yticks('');
    colormap jet

end
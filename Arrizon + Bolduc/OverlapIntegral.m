% Ralf Mouthaan
% University of Cambridge
% July 2019
% Calculation of overlap integral

function [RetVal] = OverlapIntegral(M1, M2)

    RetVal = abs(sum(sum(M1.*conj(M2))))^2;
    RetVal = RetVal/sum(sum(abs(M1.*conj(M1))));
    RetVal = RetVal/sum(sum(abs(M2.*conj(M2))));
    RetVal = sqrt(RetVal);

end
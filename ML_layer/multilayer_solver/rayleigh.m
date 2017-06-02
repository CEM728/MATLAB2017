%
%  D   : U_RB or D_RB fictitous interface for enforcing radiation condition
%  d   : period
%  M   : number of points in D
%  BR  : truncation of Rayleigh-Bloch expansion (from -BR to BR)
%  k1x : k1*cos(ki)
%  k   : wave number
%
function [W Wn] = rayleigh(D, k1x, k, d, M, BR)

W  = zeros(M, BR+BR+1);
Wn = zeros(M, BR+BR+1);
for I=1:M
    for J=-BR:BR
        kappa_n = k1x+2*pi*J/(d);
        k_n     = sqrt(k^2-kappa_n^2);
        W(I,J+BR+1) =   exp(1i*kappa_n*real(D.x(I)));
        Wn(I,J+BR+1)=   imag(D.nx(I))*1i*k_n*W(I,J+BR+1);  % normal deriv (downwards sense)
     end
end

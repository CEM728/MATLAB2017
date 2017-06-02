%
%Compute field using Rayleigh-Bloch expasnion 
%using the coefficient au and ad above and below the radiation line
%
function [field_up, field_down, gx_up, gy_up, gx_down, gy_down] = compute_field_n_layer(au, ad, k1x, k1y, k1, k2, d, up_height, down_height, BR)
nx = 150; 
gx_up =  -0.5:1/(nx-1):0.5;
gy_up =  up_height:1.0/(nx-1):up_height+1.0;
[xx, yy] =  meshgrid(gx_up,gy_up); 
zz1      = (xx+1i*yy); 
clear xx yy

test_point1.x  = zz1(:); 
test_point1.nx = [];

gx_down =  -0.5:1/(nx-1):0.5;
gy_down = down_height-1.0:1.0/(nx-1):down_height;
[xx, yy] =  meshgrid(gx_down,gy_down); 
zz2      = (xx+1i*yy); 
clear xx yy

test_point2.x  = zz2(:); 
test_point2.nx = [];

field_up = zeros(nx*nx,1);
field_down = zeros(nx*nx,1);
for J=-BR:BR
    kappa_n = k1x+2*pi*J/(d);
    ku_n     = sqrt(k1^2-kappa_n^2);
    kd_n     = sqrt(k2^2-kappa_n^2);
    field_up   = field_up+au(J+BR+1)*exp(1i*kappa_n.*real(test_point1.x)).*exp(1i*ku_n.*(imag(test_point1.x)-up_height));
    field_down = field_down+ad(J+BR+1)*exp(1i*kappa_n.*real(test_point2.x)).*exp(-1i*kd_n.*(imag(test_point2.x)-down_height));
end

field_up = field_up+exp(1i*k1x*real(test_point1.x)+1i*k1y*imag(test_point1.x));


field_up = reshape(field_up, [nx,nx]);
field_down = reshape(field_down, [nx,nx]);
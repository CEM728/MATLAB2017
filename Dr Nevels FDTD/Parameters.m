function Parameters()
clear all
close all
c = 2.99792458e8; %% Speed of light
eps0 = 8.845e-12; %% Free space permittivity
epsr = 4; %% Relative Permitivitty of medium
xmu0 = 1.25663706e-6; %% Free space permeability
a = 10; b = 10; %% Space Dimensions
nx = 14; ny = 14; nx2 = nx/2; ny2 = ny/2; nt = 100; nskip = 5;
% % % % % % % % % % 
% % % % % % % % % %  Initializing
Ez = zeros( nx, ny ); %%  Z E-field initialize to zero
Hx = zeros( nx, ny ); %%  X H-field initialize to zero
Hy = zeros( nx, ny ); %%  Y H-field initialize to zero

dx = a/(nx -1); %% 
dy = dx;
ds = dx;
dt = ds/sqrt(2);
nsnap = nt/nskip;
% % % Field Coefficients
dte = dt/(ds*eps0);
dtm = dt/(ds*xmu0);

save fdtdcom.mat
end

    

function [Ez, Hx, Hy, dx, dy, ds, dt, nsnap] = initialize()
load fdtdcom.mat

% % % % Fields
Ez = zeros( nx, ny ); %%  Z E-field initialize to zero
Hx = zeros( nx, ny ); %%  X H-field initialize to zero
Hy = zeros( nx, ny ); %%  Y H-field initialize to zero

% % % % Spatial Increments
dx = a/(nx -1); %% cell width in x-direction
dy = dx;        %% cell width in y-direction
ds = dx;        %% standard cell width

% % % % Time increments
dt = ds/sqrt(2);%% stability condition

% % % % Number of plots on-screen
nsnap = nt/nskip;
end
clear all

% Physical constants
eps0 = 8.8541878e-12;         % Permittivity of vacuum
mu0  = 4e-7 * pi;             % Permeability of vacuum
c0   = 299792458;             % Speed of light in vacuum

% Parameter initiation
Lx = .05; Ly = .04; Lz = .03; % Cavity dimensions in meters
Nx =  40; Ny =  30; Nz =  20; % Number of cells along each axis
Cx = Nx / Lx;                 % Inverse cell dimensions
Cy = Ny / Ly;
Cz = Nz / Lz;
Nt = (2^1)*8192;              % Number of time steps
Dt = 1/(c0*norm([Cx Cy Cz])); % Time step

% Allocate field matrices
Ex = zeros(Nx  , Ny+1, Nz+1);
Ey = zeros(Nx+1, Ny  , Nz+1);
Ez = zeros(Nx+1, Ny+1, Nz  );
Hx = zeros(Nx+1, Ny  , Nz  );
Hy = zeros(Nx  , Ny+1, Nz  );
Hz = zeros(Nx  , Ny  , Nz+1);

% Allocate time signals
Et = zeros(Nt,3);

% Initiate fields with noise (except on the boundary)
Ex( :  , 2:Ny, 2:Nz) = rand(Nx  , Ny-1, Nz-1) - 0.5;
Ey(2:Nx,  :  , 2:Nz) = rand(Nx-1, Ny  , Nz-1) - 0.5;
Ez(2:Nx, 2:Ny,  :  ) = rand(Nx-1, Ny-1, Nz  ) - 0.5;

% Time stepping
for n = 1:Nt
  
  if (mod(n,1000) == 0)
      str = sprintf('n = %i going to %i ', n,Nt);
      disp(str), drawnow
  end
    
  % Update H everywhere
  Hx = Hx + (Dt/mu0)*(diff(Ey,1,3)*Cz - diff(Ez,1,2)*Cy);
  Hy = Hy + (Dt/mu0)*(diff(Ez,1,1)*Cx - diff(Ex,1,3)*Cz);
  Hz = Hz + (Dt/mu0)*(diff(Ex,1,2)*Cy - diff(Ey,1,1)*Cx);

  % Update E everywhere except on boundary
  Ex(:,2:Ny,2:Nz) = Ex(:,2:Ny,2:Nz) + (Dt /eps0) * ...
      (diff(Hz(:,:,2:Nz),1,2)*Cy - diff(Hy(:,2:Ny,:),1,3)*Cz);
  Ey(2:Nx,:,2:Nz) = Ey(2:Nx,:,2:Nz) + (Dt /eps0) * ...
      (diff(Hx(2:Nx,:,:),1,3)*Cz - diff(Hz(:,:,2:Nz),1,1)*Cx);
  Ez(2:Nx,2:Ny,:) = Ez(2:Nx,2:Ny,:) + (Dt /eps0) * ...
      (diff(Hy(:,2:Ny,:),1,1)*Cx - diff(Hx(2:Nx,:,:),1,2)*Cy);

  % Sample the electric field at chosen points
  Et(n,:) = [Ex(4,4,4) Ey(4,4,4) Ez(4,4,4)];
end



% --- post-processing of numerical spectrum ---
Df       = 1/(Nt*Dt);     % Spacing between frequency points
freq_vtr = Df*(0:Nt-1);   % Frequency sampling points

fet      = abs(fft(Et(:,1) + Et(:,2) + Et(:,3) ) );

idx      = find(freq_vtr < 1e10);    % Get values in the frequency interval of interest
idx      = setdiff(idx,1);           % Remove the static solution


% --- analytical eigenfrequencies ---
fv = [];
nm = [];
for nx = 0:10
    for ny = 0:10
        for nz = 0:10
            if (((nz >= 1) & ~(nx == 0 & ny == 0)) | ...
                ((nx >= 1) & (ny >= 1)))
                fv = [fv; ...
                    c0/(2*pi)*sqrt( ...
                    (pi*nx/Lx)^2 ...
                    + (pi*ny/Ly)^2 ...
                    + (pi*nz/Lz)^2)];
                nm = [nm; nx ny nz];
            end
        end
    end
end

[fv, id] = sort(fv);
nm = nm(id,:);


% --- plotting ---
plot(freq_vtr/1e9, fet, 'k-'), hold on

ymax     = ceil(max(fet(idx)));
for idx = 1:100
    plot(fv(idx)*[1 1]/1e9, ymax*[0 1], 'k--')
end

xlabel('f [GHz]')
ylabel('E [V/m]')
axis([0 10 0 ymax])


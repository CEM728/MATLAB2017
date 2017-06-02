clear all
close all
figure('Name', 'sdf','Position',[100 100 850 600])
% % % % % % % % % % PARAMETERS

c = 2.99792458e8; %% Speed of light
eps0 = 8.854187817e-12; %% Free space permittivity
epsr = 4; %% Relative Permitivitty of medium
xmu0 = 4*pi*1e-7; %% Free space permeability
a = 10; b = 10; %% Space Dimensions
nx = 150; ny = 150; 
nx2 = ceil(nx/2); ny2 = ceil(ny/2); 
nt = 100; nskip = 5;
% % % % % % % % % % 

% % % % % % % % % %  INITIALIZE
Ez = zeros( nx, ny ); %%  Z E-field initialize to zero
Hx = zeros( nx, ny ); %%  X H-field initialize to zero
Hy = zeros( nx, ny ); %%  Y H-field initialize to zero
Ezsz = zeros(1,nt);   %%  Shape of the Source function with time
dx = a/(nx - 1); %% 
dy = dx;
ds = dx;
dt = 1.0005*ds/(c*sqrt(2)); %% Stability Condition
nsnap = ceil(nt/nskip);

% % % % % % % % % % MEDIUM DEFINITION
% % % This function sets up the field coefficients subject to the medium
% parameters
% % % Field Coefficients
dte = ones(nx,ny)*dt/(ds*eps0);
dtm = ones(nx,ny)*dt/(ds*xmu0);
Da = ones(nx,ny);
Db = dtm;
Ca = ones(nx,ny);
Cb = dte;
% % % % % % % % % % 

% % % % This part is used to create a PEC box in computational space
dte(100:110, 100:110) = 0;
dte(100:110, 100:110) = 0;
Da(100:110, 100:110) = 0;
Db = dtm;
Ca(100:110, 100:110) = 0;
Cb = dte;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         


% % % % % % % % % MAIN PROGRAM
for n = 1:nt
    for i = 2 : nx - 1
        for j = 2 : ny - 1
             if (i == nx2 && j == ny2)
                Ez(nx2,ny2) = Ez_inc(n);
                Ezsz(n) = Ez(i,j);
%              elseif (i == 2 || j == 2 || i == nx-1 || j == ny-1 )   %%% PEC Boundary Conditions
%                  Ez(i,j) = 0;
                 
             else
                 Ez(i,j) = Ez(i,j)*Ca(i,j) + Cb(i,j)*(Hy(i,j) - Hy(i-1,j)...
                - (Hx(i,j) - Hx(i,j-1)));
                 
             end
        end
    end
    
    for i = 1 : nx 
        for j = 1 : ny - 1
                      
             Hx(i,j) = Hx(i,j)*Da(i,j) - Db(i,j)*(Ez(i,j+1) - Ez(i,j));
        end
    end
    for i = 1 : nx - 1
        for j = 1 : ny 

             Hy(i,j) = Hy(i,j)*Da(i,j) + Db(i,j)*(Ez(i+1,j) - Ez(i,j));
        end
    end

%     V = interp2(abs(Ez),1);
     hSurf = surf((Ez),'EdgeColor','interp','LineStyle','none','BackFaceLighting', 'reverselit',...
         'FaceLighting','phong','FaceColor','interp','SpecularColorReflectance',.1);
     camlight ('right', 'infinite');
%      shading faceted
%      light
%      lighting flat
     
% set(gca, 'CameraPosition', [45 35 9.8])
%     hSurf = surf(V);
     set(gcf,'Color','white');
     set (gca,'FontName','times new roman')
%     set(hSurf,'Color','black','LineWidth',1.4)
     xlabel(' X Cells','FontSize',11,'Interpreter','latex','FontName','times new roman')
     ylabel(' Y Cells','FontSize',11,'Interpreter','latex','FontName','times new roman')
     title(['E-field at Time Step ',int2str(n)],'FontSize',11,'Interpreter','latex','FontName','times new roman')
%    axis([400 780 -2 2.2])
     axis tight%     colorbar
    grid off
    drawnow
%     axesLabelsAlign3D
    box off
%     if n == 20
%         export_fig E_field_20 '-pdf'  -nocrop -r300 -transparent  -native -painters -q101
%         saveas(hSurf,'imag1.eps')
%     elseif n == 60
%         export_fig E_field_60 '-pdf'  -nocrop -r300 -transparent  -native -painters -q101
%         saveas(hSurf,'imag2.eps')
%     elseif n == 80
%         export_fig E_field_80 '-pdf'  -nocrop -r300 -transparent  -native -painters -q101
%         saveas(hSurf,'imag3.eps')
%     elseif n == 100
%         export_fig E_field_100 '-pdf'  -nocrop -r300 -transparent  -native -painters -q101
%         saveas(hSurf,'imag4.eps')
%     end
        
    
%     view([0,0,1])
     
%     colormap('hot')
    
    M(:,n) = getframe ;
    
end
movie(M,1,30);




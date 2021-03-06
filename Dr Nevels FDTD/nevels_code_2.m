clear all
close all
% figure('Name', 'Field Plot','Position',[0 0 1366 768])
% % % % % % % Parameters
c = 2.99792458e8;
xmu = 4*pi*1e-7;
eps0 = 8.854187817e-12;
asize = 5;
nx = 200;      %% Number of cells in x-direction
ny = 200;     %% Number of cells in y-direction
nt = 400;     %% Number of time steps
mxst = 85;    %% Start of PEC section in x-direction
mxnd = 117;    %% End of PEC section in x-direction
myst = 85;    %% Start of PEC section in y-direction
mynd = 117;    %% End of PEC section in y-direction
ky = 49;

% % % % % % % Initialize
Ez = zeros(nx,ny); %% z-component of E-field
Hx = zeros(nx,ny); %% x-component of H-field
Hy = zeros(nx,ny); %% y-component of H-field

mediaEz = ones(nx,ny);
mediaHx = ones(nx,ny);
mediaHy = ones(nx,ny);

Ez1 = zeros(1,nx);
Ez2 = zeros(1,nx);
Ez3 = zeros(1,nx);
Ez4 = zeros(1,nx);
Ez5 = zeros(1,nx);
Ez6 = zeros(1,nx);
Ez7 = zeros(1,nx);

Ca = zeros(2,1);
Cb = zeros(2,1);
Da = zeros(2,1);
Db = zeros(2,1);

ds = asize/(mxnd - mxst - 1); %% Length Increment
dt = ds/(2*c);%% Time increment for 2-D

iflaga = 2;   %% 1 if Free space; 2 if PEC

% % % % % % % % Structure Definition
if (iflaga == 2)
    
    for  i = 1:nx
        for j = 1:ny
            if (i >= mxst && i <= mxnd)
                if ( j >= myst && j <= mynd)
                    mediaEz(i,j) = 2;
                end
            end
        end
    end

    for  i = 1:nx
        for j = 1:ny
            if (i >= mxst && i <= mxnd)
                if ( j >= myst && j <= mynd-1)
                    mediaHx(i,j) = 2;
                end
            end
        end
    end

    for  i = 1:nx
        for j = 1:ny
            if (i >= mxst && i <= mxnd-1)
                if ( j >= myst && j <= mynd)
                    mediaHy(i,j) = 2;
                end
            end
        end
    end
end

% % % % % % % % Field Coefficients
dte = dt/(ds*eps0);
dtm = dt/(ds*xmu);
Da(1) = 1;
Db(1) = dtm;
Ca(1) = 1;
Cb(1) = dte;
Da(2) = 0;
Db(2) = 0;
Ca(2) = 0;
Cb(2) = 0;

% % % % % % % % Main Program
for n = 1:nt
    
% % %     Compute z-component of E-field
    
    for i = 1 : nx      
        for j = 1 : ny 
            m  = mediaEz(i,j);
             if (i == 1)   %% Incident Field Source Excitation
                Ez(1,j) = Ez_inc(n);
             elseif (i >= 2 && j >=2)
                Ez(i,j) = Ez(i,j)*Ca(m) + Cb(m)*(Hy(i,j) - Hy(i-1,j)...
                - (Hx(i,j) - Hx(i,j-1)));
             elseif (j == 1) %% Field at the bottom edge of the boundary
                     Ez(i,j) = Ez(i,j)*Ca(m) + Cb(m)*(Hy(i,j) - Hy(i-1,j)...
                         - Hx(i,j));                      
             end
        end
    end
    
% % %     Compute x-component of H-field
    
    for i = 1 : nx 
        for j = 1 : ny - 1
             m = mediaHx(i,j);         
             Hx(i,j) = Hx(i,j)*Da(m) - Db(m)*(Ez(i,j+1) - Ez(i,j));
        end
    end
    
% % %     Compute y-component of H-field
    
    for i = 1 : nx - 1
        for j = 1 : ny 
             m = mediaHy(i,j);
             Hy(i,j) = Hy(i,j)*Da(m) + Db(m)*(Ez(i+1,j) - Ez(i,j));
        end
    end
 
% % % % % % %  PLOT Routine
 
%   hold on
 hSurf = surf((Ez),'EdgeColor','none','LineStyle','none','BackFaceLighting', 'reverselit',...
     'FaceLighting','phong','FaceColor','interp','SpecularColorReflectance',.1);
 caxis manual
 rectangle('Position',[mxst,myst,32,32],'FaceColor','r','LineWidth',2)
view([ 0 0 1])
%  brighten(1)
%  camlight ('headlight', 'infinite');
%    hold off
%      shading faceted
%      light
%      lighting flat

% set(gca, 'CameraPosition', [45 35 9.8])
%     hSurf = surf(V);
axis equal
 set(gcf,'Color','white');
 set (gca,'FontName','times new roman')
%     set(hSurf,'Color','black','LineWidth',1.4)
 xlabel(' X Cells','FontSize',11,'Interpreter','latex','FontName','times new roman')
 ylabel(' Y Cells','FontSize',11,'Interpreter','latex','FontName','times new roman')
 title(['E-field at Time Step ',int2str(n)],'FontSize',11,'Interpreter','latex','FontName','times new roman')
%    axis([400 780 -2 2.2])
% axis ([ 0 nx 0 ny -1 1]) %     colorbar
grid off
% view([10,010,5])

%     drawnow
%     axesLabelsAlign3D
box off

%     if n == 20
%         export_fig E1_field_20 '-pdf'  -nocrop -r300 -transparent  -native -painters -q101
%         saveas(hSurf,'imag11.eps')
%     elseif n == 60
%         export_fig E1_field_60 '-pdf'  -nocrop -r300 -transparent  -native -painters -q101
%         saveas(hSurf,'imag22.eps')
%     elseif n == 120
%         export_fig E1_field_80 '-pdf'  -nocrop -r300 -transparent  -native -painters -q101
%         saveas(hSurf,'imag33.eps')
%     elseif n == 160
%         export_fig E1_field_100 '-pdf'  -nocrop -r300 -transparent  -native -painters -q101
%         saveas(hSurf,'imag42.eps')
%     end
if n == 10
    Ez1 = Ez(:,30);
elseif n == 30
    Ez2 = Ez(:,30);
elseif n == 60
    Ez3 = Ez(:,30);
elseif n == 90
    Ez4 = Ez(:,30);
elseif n == 120
    Ez5 = Ez(:,30);
elseif n == 150
    Ez6 = Ez(:,30);
elseif n == 180
    Ez7 = Ez(:,30);
end



%     colormap('hot')

M(:,n) = getframe ;
        
end
figure('Name', 'Time Steps of Electric Fields','Position',[100 100 850 600])
set(gcf,'Color','white');
title(['E-field at Time Step ',int2str(n)],'FontSize',11,'Interpreter','latex','FontName','times new roman')
subplot(7,1,1)
plot(Ez1,'Color','black','LineWidth',1.4)
axis([0 100 -1 1])
set(gca,'XTickLabel',{})
set (gca,'FontName','times new roman')
subplot(7,1,2)
plot(Ez2,'Color','black','LineWidth',1.4)
axis([0 100 -.5 .5])
set(gca,'XTickLabel',{})
set (gca,'FontName','times new roman')
subplot(7,1,3)
plot(Ez3,'Color','black','LineWidth',1.4)
axis([0 100 -.5 .5])
set(gca,'XTickLabel',{})
set (gca,'FontName','times new roman')
subplot(7,1,4)
plot(Ez4,'Color','black','LineWidth',1.4)
axis([0 100 -.5 .5])
set(gca,'XTickLabel',{})
set (gca,'FontName','times new roman')
subplot(7,1,5)
plot(Ez5,'Color','black','LineWidth',1.4)
axis([0 100 -.5 .5])
set(gca,'XTickLabel',{})
set (gca,'FontName','times new roman')
subplot(7,1,6)
plot(Ez6,'Color','black','LineWidth',1.4)
axis([0 100 -.5 .5])
set(gca,'XTickLabel',{})
set (gca,'FontName','times new roman')
subplot(7,1,7)
plot(Ez7,'Color','black','LineWidth',1.4)
axis([0 100 -.5 .5])
set(gcf,'Color','white');
set (gca,'FontName','times new roman')
axis([0 100 -.5 .5])
xlabel(' Cells in x-direction','FontSize',11,'Interpreter','latex','FontName','times new roman')
export_fig E1_fields '-pdf'  -nocrop -r300 -transparent  -native -painters -q101
 
%    axis([400 780 -2 2.2])
axis tight%     colorbar
grid off
figure('Name', 'Movie','Position',[10 10 850 600])
set(gcf,'Color','white');
set (gca,'FontName','times new roman')
axis ([ 0 100 0 80 -1 1])
movie(M,1,30);
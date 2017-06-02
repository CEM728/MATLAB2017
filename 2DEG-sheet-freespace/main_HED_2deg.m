% This program computes the Sommerfeld integral for a Horizontal Electric
% Dipole
clear all; 
close all
tic
tol = 1e-15; % tolerance of the routine
num = 100; %Size of the arrays 200 is a good number
%% Global Parameters
global i % index number of the distance array
global p % distance
global a % Breakpoint location
global nu % Switch for TE/TM case (alpha = 0 -> TE, else -> TM)
global h
global f omega 
global ep1 ep2 
global k1 k2 
global ep0 mu0

%% GaN
% f = .79e12;
% f = .88e12; % good conductivity
% f = .81e12; % good conductivity
% f = .875e12; % good conductivity
% f = .583e12; % good conductivity
% f = .2663e12; % good conductivity

%% GaN
f = 5.6e12;


omega = 2*pi*f;
ep1 = 1;
% 3 K
ep2 = -13.46 - 1i*0.11;
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

k1 = omega*sqrt(mu0*ep0*ep1);
k2 = omega*sqrt(mu0*ep0*ep2);

load besselzeros.mat First_15_zeros_J0 First_15_zeros_J1

a = 2*k1; % Set breakpoint
% p = linspace(1e-3/k1,1e2/k1, num); % Define distance array
p = 1/k1*logspace(-3,2, num); % Define distance array

% TE case
nu = 0;


% % TM case
% nu = 1;
% 

% Define bessel functions
S_0 = @(kp) besselj(0, kp * p(i));
S_1 = @(kp) besselj(1, kp * p(i));



for i = 1 : length(p)
    if nu == 0
        q = First_15_zeros_J0/p(i);
    else
        q = First_15_zeros_J1/p(i);
    end
    % Avoid branch points
    if nu == 0 % TE case
        
        h = 1;
        val_1(i) = TanhSinhQuad_2deg(0, k1 + .001i, tol); % Integrate upto k through DE
        h = 1;
        val_2(i) = TanhSinhQuad_2deg(k1 + .001i, a, tol); % Integrate k upto a through DE
        h = 1;
        val_3(i) = PE_Levin_2deg(a, tol, q); % Tail through PE Levin with Lucas
    else
        h = 1;
        val_1(i) = TanhSinhQuad_2deg(.001i, k1, tol); % Integrate upto k through DE
        h = 1;
        val_2(i) = TanhSinhQuad_2deg(k1, a, tol); % Integrate k upto a through DE
        h = 1;
        val_3(i) = PE_Levin_2deg(a, tol, q); % Tail through PE Levin with Lucas
    end
    
    val(i) = val_1(i) + val_2(i) + val_3(i);
end

clf
% Individual Contribution
figure (1)

N = 3; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = loglog(p*k1, abs(val_1/k1), 'linewidth',1.3);
hold on
h2 = loglog(p*k1, abs(val_2/k1), 'linewidth',1.3);
h3 = loglog(p*k1, abs(val_3/k1), 'linewidth',1.3);

loglog(p*k1, abs(val_1/k1), 's', 'markersize',4);
loglog(p*k1, abs(val_2/k1), 's', 'markersize',4);
loglog(p*k1, abs(val_3/k1), 's', 'markersize',4);
xlabel('$k_1\rho$','interpreter','latex')
ylabel('$I(z, \rho, \tau)$','interpreter','latex')
legend([h1 h2 h3],{'DE Rule 0 to k', 'DE Rule k to a', 'PE Rule'},...
     'location','northwest');
if nu == 0
    title('TE case');
else
    title('TM case');
end
box on
set(gcf,'color','white');
set(gca,'FontName','times new roman','FontSize',15);
hold off

% Tail Contribution
figure (2)
N = 2; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h4 = loglog(p*k1, abs(val_3/k1), 'linewidth',1.4);
hold on
% loglog(p*k1, abs(val_3/k1), 's', 'markersize',4);
xlabel('$k_1\rho$','interpreter','latex')
set(gca,'FontName','times new roman','FontSize',15);


box on
set(gcf,'color','white');
hold off
% cleanfigure();
% if nu == 0
% 
%     ylabel('$G_{xx}^{\mathrm{A}}$','interpreter','latex')
%     matlab2tikz('filename',sprintf('figures/Gxx_20.tikz'),'showInfo', false)
%     
% %         ylabel('$G_{zz}^{\mathrm{A}}$','interpreter','latex')
% %         matlab2tikz('filename',sprintf('figures/Gzz_20.tikz'),'showInfo', false)
% else
%     
%     ylabel('$G_{zx}^{\mathrm{A}}$','interpreter','latex')
%     matlab2tikz('filename',sprintf('figures/Gzx_20.tikz'),'showInfo', false)
% end



% Overall integrals
figure(3)

N = 2; % Number of colors to be used
% Use Brewer-map color scheme 'Set1'
axes('ColorOrder',brewermap(N,'set1'),'NextPlot','replacechildren')

h5 = loglog(p*k1, abs(val)/k1, 'linewidth',1.4);
hold on
h6 = loglog(p*k1, abs(exp(-1i*k1*p)./(4*pi*p))./k1, 'linewidth',1.);
% loglog(p*k1, abs(val)/k1, 's', 'markersize',4);
xlabel('$k_1\rho$','interpreter','latex')
% ylabel('$I(z, \rho, \tau)$','interpreter','latex')

% if nu == 0
%     title('TE case');
% else
%     title('TM case');
% end
box on
xlim([1e-3 1e2])
set(gcf,'color','white');
set(gca,'FontName','times new roman','FontSize',15);
legend([h5 h6],{'\mathcal{S}_n\{ \cdot \}', 'Freespace GF'},...
     'location','northeast');
hold off
if nu == 0

%     ylabel('$G_{xx}^{\mathrm{A}}$','interpreter','latex')
%     matlab2tikz('filename',sprintf('figures/Gxx_gan_sheet_56_rt.tikz'),'showInfo', false)
   
        ylabel('$G_{zz}^{\mathrm{A}}$','interpreter','latex')
        matlab2tikz('filename',sprintf('figures/Gzz_gan_sheet_56_rt.tikz'),'showInfo', false)
else
    
    ylabel('$G_{zx}^{\mathrm{A}}$','interpreter','latex')
    matlab2tikz('filename',sprintf('figures/Gzx_gan_sheet_56_rt.tikz'),'showInfo', false)
end
toc
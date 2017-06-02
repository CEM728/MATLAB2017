close all
clf;clear all
% Hasan Tahir Abbas
% 10/14/15
%
%
% Based on the paper:
% Yoon et. al., "Plasmonics with
% two-dimensional conductors", 
% Phil. Trans. R. Soc. A 372: 20130104.
% http://dx.doi.org/10.1098/rsta.2013.0104
%
%
%% Initial Parameters
eps0 = 8.85418782e-12; % Free-space permitivitty
c = 2.99792458e8; % speed of light
e = 1.60217662e-19; % electron charge
m_e = 9.10938356e-31;  % electron mass
m_star = .067*m_e; % Effective mass of GaAs/AlGaAs
h = 6.62607004e-34; % Planck's constant
h_bar = h/(2*pi);
% Assumptions
n_2d = 1e16; % Sheet Charge Density in /m^2
kf = sqrt ( 2*pi*n_2d ); % Fermi wavenumber
vf = h_bar*kf/m_star; % Fermi velocity
Ef = h_bar^2 * kf^2 /(2*m_star); % Fermi Energy
K = 4; % effective dielectrinc constant

kp = 0 : 10e4 : 2*kf;
omega = sqrt( e^2*Ef /(2*eps0*K*pi*h_bar) * kp ...
    + 3/4 *vf^2*kp.^2);
% omega = sqrt( g*e^2*vf*kf*kp /( 8*pi*K*eps0* h_bar));
plot(kp,omega,'k','LineWidth',1.2); % Disperion Curve
hold on
box on
 plot(kp, (1e6*kp),'--k','LineWidth',1.2); % Light Line
hold off
axis tight
set (gca,'FontName','times new roman'...
    ,'FontWeight','bold','FontSize',11); % Set axes fonts to
set(gcf,'Color','white');
h = gca;
h.LineWidth = 1.2;
% h.XTick ={};
% h.YTick ={};


% Create ylabel
xlabel('Wavenumber $k_p $',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel

ylabel('Frequency $\omega $',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

plotTickLatex2D

% Create Legend
legend({'2DEG', 'Light Line'}...
    ,'Location','northwest','FontWeight','normal')
% Create title
title(['Dispersion Plot of GaAs/AlGaAs 2DEG'],...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

matlab2tikz('filename',sprintf('2DEG_dispersion.tex'));

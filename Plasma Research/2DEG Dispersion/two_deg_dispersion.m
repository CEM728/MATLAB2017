% This program arbitrarily computes and plots the dispersion curves of 3D
% and 2D plasmonic structures
% 10-05-2014
%
%
% Based on the paper:
% B. M. Santoyo, "Plasmons in one two and three dimensions", Revista
% Mexicana de Fisica pg 640-652, 1993.
%
%
clf;clear all;close all
%%
n0 = 3e12*1e-4;
m0 = 9.1e-31;
m  = .042*m0;
d = 10e-9;
e = 1.602e-19;
eps_r = 8.845e-12*13.22;
q = 0 : .01: 1;
% e = 1*ones(1,length(q));
% w_3d = sqrt(4*pi*n0*e.^2/m);
w_ungated = sqrt(e^2*n0/(m*eps_r).* q);
w_gated =sqrt(e^2*n0*d/(m*eps_r)) .*q;
w_ll = q;
str = 'Dispersion curves';
%% Plot Formatting
set (gca,'FontName','times new roman'...
    ,'FontWeight','bold','FontSize',11); % Set axes fonts to
set(gcf,'Color','white');
h = gca;
h.LineWidth = 1.2;
% figure('Name',str)
% h = plot(q, w_3d,'LineStyle','--','LineWidth',1.5);
plot(q, .4*w_ungated./max(w_ungated),'LineStyle','-','LineWidth',1.3,'Color','black'); % ungated
hold on
plot(q, .2*w_gated./max(w_gated),'LineStyle','-','LineWidth',1.3,'Color','red'); % gated 
plot(q, .5*ones(length(q)),'LineStyle','-','LineWidth',1.2,'Color','green'); % 3d
plot(q, .02+1.5*w_ll./max(w_ll),'LineStyle','--','LineWidth',1.,'Color','black'); % light line
ylim ([0 1])
set(gca,'FontSize',10, 'FontName' , 'times new roman')
% Create ylabel
xlabel('Wavenumber, $k$ [arb. units]',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',8,...
    'Interpreter','latex');

% Create xlabel

ylabel('Frequency [arb. units]',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',8,...
    'Interpreter','latex');

plotTickLatex2D
cleanfigure
% Create Legend
legend({'Ungated','Gated','3D','Light Line'}...
    ,'Location','northeast','FontWeight','bold',...
    'FontSize',10,...
    'Interpreter','latex')
% Create title
% title(['Dispersion Plot of GaAs/AlGaAs 2DEG'],...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',8,...
%     'Interpreter','tex');
box on
set(gca,'FontSize',010, 'FontName' , 'times new roman')
matlab2tikz('filename',sprintf('2DEG_disper.tex'),'interpretTickLabelsAsTex',true);
% export_fig 2DEG_disper '-eps'  -nocrop -r300 -transparent  -native -painters -q101
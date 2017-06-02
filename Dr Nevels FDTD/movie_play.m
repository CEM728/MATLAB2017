clear all
close all
load fdtd_sine_wave.mat
figure('Name', 'E-field','Position',[0 0 1366 768])
 set(gcf,'Color','white');
 set (gca,'FontName','times new roman')
 set(gca,'XTickLabel',{})
 set(gca,'YTickLabel',{})
 axis off
 movie(M,1,30)
 axis equal
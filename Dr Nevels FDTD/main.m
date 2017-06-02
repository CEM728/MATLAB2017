% % % % Here is the main program where all the functions are invoked
clear all
close all
load fdtdcom.mat
% % % % % % % 
Parameters();
initialize();
medium();
% % % % % % % 
for n = 1:nt
    for i = 2 : nx-1
        for j = 2 : ny-1
    adv_efield(i,j,n);
    adv_hfield(i,j,n);
        end
    end
end



        

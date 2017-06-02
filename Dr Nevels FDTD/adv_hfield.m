function [Hx, Hy] = adv_hfield (i,j,n)
% % % % Computes time stepped Magnetic Fields, Hx and Hy
load fdtdcom.mat
[ Da, Db, Ca, Cb] = medium();
% % % Evaluation

       Hx(i,j) = Hx(i,j)*Da(i,j) - Db(i,j)*(Ez(i,j+1) - Ez(i,j));

       Hy(i,j) = Hy(i,j)*Da(i,j) + Db(i,j)*(Ez(i+1,j) - Ez(i,j));
end
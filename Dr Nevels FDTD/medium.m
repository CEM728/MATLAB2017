function [ dte, dtm, Da, Db, Ca, Cb] = medium()
% % % This function sets up the field coefficients subject to the medium
% parameters
% % % Field Coefficients
load fdtdcom.mat
dte = ones(nx,ny)*dt/(ds*eps0);
dtm = ones(nx,ny)*dt/(ds*xmu0);
Da = ones(nx,ny);
Db = dtm;
Ca = ones(nx,ny);
Cb = dte;
end
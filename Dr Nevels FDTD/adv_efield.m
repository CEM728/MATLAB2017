function Ez = adv_efield (i,j,n)
% % % % Computes time stepped Electric Field, Ez
load fdtdcom.mat
[ Da, Db, Ca, Cb] = medium();
% % % Evaluation
        if (i == nx2 && j == ny2)
            Ez(i,j) = Ez_inc(n);
        else
            Ez(i,j) = Ez(i,j)*Ca(i,j) + Cb(i,j)*(Hy(i,j) - Hy(i-1,j)...
                - (Hx(i,j) - Hx(i,j-1)));
        end
end
function [B2, B2n] = construct_B2(i, k, Layer_part, MFS, nop)
B2=[];B2n=[];
for j=1:nop(i)
    [B_temp, B_temp_n] = double_layer_seg(k(i), Layer_part{i,j}, MFS(i), 0, 0);
    [BS_temp, BS_temp_n] = single_layer_seg(k(i), Layer_part{i,j}, MFS(i), 0, 0);
    B2 = [B2;B_temp+1i*k(i)*BS_temp];
    B2n = [B2n;B_temp_n+1i*k(i)*BS_temp_n];
end
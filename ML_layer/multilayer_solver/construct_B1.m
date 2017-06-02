function [B1, B1n] = construct_B1(i, k, Layer_part, MFS, nop)
B1=[];B1n=[];
for j=1:nop(i-1)
    [B_temp, B_temp_n] = double_layer_seg(k(i), Layer_part{i-1,j}, MFS(i), 0, 0);
    [BS_temp, BS_temp_n] = single_layer_seg(k(i), Layer_part{i-1,j}, MFS(i), 0, 0);
    B1 = [B1;B_temp+1i*k(i)*BS_temp];   
    B1n = [B1n;B_temp_n+1i*k(i)*BS_temp_n];
end
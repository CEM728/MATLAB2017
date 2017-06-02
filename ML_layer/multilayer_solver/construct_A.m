function [AD11, AS11, AD11n, AS11n, AD12, AS12, AD12n, AS12n, AD21, AS21, AD21n, AS21n] = construct_A(i, k, Layer_part, nei, alpha, d, MI, MI_total, nop)

%(i,i) Diagonal Block
AD11  = zeros(MI_total(i),MI_total(i));
AS11  = zeros(MI_total(i),MI_total(i));
AD11n = zeros(MI_total(i),MI_total(i));
AS11n = zeros(MI_total(i),MI_total(i));

%(i,i+1) Off-diagonal Block
AD12  = zeros(MI_total(i),MI_total(i+1));
AS12  = zeros(MI_total(i),MI_total(i+1));
AD12n = zeros(MI_total(i),MI_total(i+1));
AS12n = zeros(MI_total(i),MI_total(i+1));

%(i-1,i) Off-diagonal Block
AD21  = zeros(MI_total(i),MI_total(i-1));
AS21  = zeros(MI_total(i),MI_total(i-1));
AD21n = zeros(MI_total(i),MI_total(i-1));
AS21n = zeros(MI_total(i),MI_total(i-1));

M11_index = [0,cumsum(MI(i,:))];
M12_index = [0,cumsum(MI(i+1,:))];
M21_index  = [0,cumsum(MI(i-1,:))];

for j=-nei:nei
    a = j*d;
    %(i,i)-Diagonal Block
    for ii=1:nop(i)
        for jj=1:nop(i)
            AD_temp = zeros(MI(i,ii), MI(i,jj));
            AD_temp_n = zeros(MI(i,ii), MI(i,jj));
            AS_temp = zeros(MI(i,ii), MI(i,jj));
            AS_temp_n = zeros(MI(i,ii), MI(i,jj));
            
            [AD_temp, AD_temp_n] = double_layer_seg(k(i), Layer_part{i,ii}, Layer_part{i,jj}, 1, a);
            [AS_temp, AS_temp_n] = single_layer_seg(k(i), Layer_part{i,ii}, Layer_part{i,jj}, 1, a);
   
            
            AD1( M11_index(ii)+1:M11_index(ii+1),M11_index(jj)+1:M11_index(jj+1)) = AD_temp;
            AD1n( M11_index(ii)+1:M11_index(ii+1),M11_index(jj)+1:M11_index(jj+1)) = AD_temp_n;
            AS1( M11_index(ii)+1:M11_index(ii+1),M11_index(jj)+1:M11_index(jj+1)) = AS_temp;
            AS1n( M11_index(ii)+1:M11_index(ii+1), M11_index(jj)+1:M11_index(jj+1)) = AS_temp_n;
            
            
            AD_temp = zeros(MI(i,ii), MI(i,jj));
            AD_temp_n = zeros(MI(i,ii), MI(i,jj));
            AS_temp = zeros(MI(i,ii), MI(i,jj));
            AS_temp_n = zeros(MI(i,ii), MI(i,jj));
                        
            [AD_temp, AD_temp_n] = double_layer_seg(k(i+1), Layer_part{i,ii}, Layer_part{i,jj}, 1, a);
            [AS_temp, AS_temp_n] = single_layer_seg(k(i+1), Layer_part{i,ii}, Layer_part{i,jj}, 1, a);
            
            AD2(M11_index(ii)+1:M11_index(ii+1),M11_index(jj)+1:M11_index(jj+1)) = AD_temp;
            AD2n(M11_index(ii)+1:M11_index(ii+1),M11_index(jj)+1:M11_index(jj+1)) = AD_temp_n;
            AS2(M11_index(ii)+1:M11_index(ii+1),M11_index(jj)+1:M11_index(jj+1)) = AS_temp;
            AS2n(M11_index(ii)+1:M11_index(ii+1),M11_index(jj)+1:M11_index(jj+1)) = AS_temp_n;
            
            
        end
    end
    AD11  = AD11+(alpha^j)*(AD1-AD2);
    AS11  = AS11+(alpha^j)*(AS1-AS2);   
    AD11n = AD11n+(alpha^j)*(AD1n-AD2n);
    AS11n = AS11n+(alpha^j)*(AS1n-AS2n);
    
    for ii=1:nop(i)
        for jj=1:nop(i+1)
            
            AD_temp = zeros(MI(i,ii), MI(i+1,jj));
            AD_temp_n = zeros(MI(i,ii), MI(i+1,jj));
            AS_temp = zeros(MI(i,ii), MI(i+1,jj));
            AS_temp_n = zeros(MI(i,ii), MI(i+1,jj));
            
            [AD_temp, AD_temp_n] = double_layer_seg(k(i+1), Layer_part{i,ii},  Layer_part{i+1,jj},   1, a);
            [AS_temp, AS_temp_n] = single_layer_seg(k(i+1), Layer_part{i,ii},  Layer_part{i+1,jj},   1, a);
            
            AD2C(M11_index(ii)+1:M11_index(ii+1), M12_index(jj)+1:M12_index(jj+1)) = AD_temp;
            AD2Cn(M11_index(ii)+1:M11_index(ii+1), M12_index(jj)+1:M12_index(jj+1)) = AD_temp_n;
            AS2C(M11_index(ii)+1:M11_index(ii+1), M12_index(jj)+1:M12_index(jj+1)) = AS_temp;
            AS2Cn(M11_index(ii)+1:M11_index(ii+1), M12_index(jj)+1:M12_index(jj+1)) = AS_temp_n;

            
            
        end
    end
    AD12  = AD12+(alpha^j)*AD2C;
    AS12  = AS12+(alpha^j)*AS2C;
    AD12n = AD12n+(alpha^j)*AD2Cn;
    AS12n = AS12n+(alpha^j)*AS2Cn;
        
    for ii=1:nop(i)
        for jj=1:nop(i-1)
            
            AD_temp = zeros(MI(i,ii), MI(i-1,jj));
            AD_temp_n = zeros(MI(i,ii), MI(i-1,jj));
            AS_temp = zeros(MI(i,ii), MI(i-1,jj));
            AS_temp_n = zeros(MI(i,ii), MI(i-1,jj));
            
            [AD_temp, AD_temp_n] = double_layer_seg(k(i), Layer_part{i,ii},  Layer_part{i-1,jj},   1, a);
            [AS_temp, AS_temp_n] = single_layer_seg(k(i), Layer_part{i,ii},  Layer_part{i-1,jj},   1, a);
            
            
            AD3C(M11_index(ii)+1:M11_index(ii+1), M21_index(jj)+1:M21_index(jj+1)) = AD_temp;
            AD3Cn(M11_index(ii)+1:M11_index(ii+1), M21_index(jj)+1:M21_index(jj+1)) = AD_temp_n;
            AS3C(M11_index(ii)+1:M11_index(ii+1), M21_index(jj)+1:M21_index(jj+1)) = AS_temp;
            AS3Cn(M11_index(ii)+1:M11_index(ii+1), M21_index(jj)+1:M21_index(jj+1)) = AS_temp_n;
        end
    end
    AD21  = AD21+(alpha^j)*AD3C;
    AS21  = AS21+(alpha^j)*AS3C;
    AD21n = AD21n+(alpha^j)*AD3Cn;
    AS21n = AS21n+(alpha^j)*AS3Cn;
end
    
    
    

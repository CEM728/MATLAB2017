function C = construct_C(i, k, Layer_part, L, R, nop, MI, MI_total, ML, nei, d, alpha)

    CDL1  = zeros(ML, MI_total(i-1));
    CDL1n = zeros(ML, MI_total(i-1));
    CSL1  = zeros(ML, MI_total(i-1));
    CSL1n = zeros(ML, MI_total(i-1));
    CDR1  = zeros(ML, MI_total(i-1));
    CDR1n = zeros(ML, MI_total(i-1));
    CSR1  = zeros(ML, MI_total(i-1));
    CSR1n = zeros(ML, MI_total(i-1));

    M1_index = [0,cumsum(MI(i-1,:))];
    for ii=1:nop(i-1)
        [CDL_temp1, CDL_temp1_n] = double_layer_seg_side(k(i), L(i-1), Layer_part{i-1,ii}, 0,  -nei*d);
        CDL1(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDL_temp1;
        CDL1n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDL_temp1_n;
    
        [CSL_temp1, CSL_temp1_n] = single_layer_seg_side(k(i), L(i-1), Layer_part{i-1,ii}, 0,  -nei*d);
        CSL1(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSL_temp1;
        CSL1n(1:ML,M1_index(ii)+1:M1_index(ii+1)) = CSL_temp1_n;
    
        [CDR_temp1, CDR_temp1_n] = double_layer_seg_side(k(i), R(i-1), Layer_part{i-1,ii}, 0,   nei*d);
        CDR1(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDR_temp1;
        CDR1n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDR_temp1_n;
    
        [CSR_temp1, CSR_temp1_n] = single_layer_seg_side(k(i), R(i-1), Layer_part{i-1,ii}, 0,   nei*d);
        CSR1(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSR_temp1;
        CSR1n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSR_temp1_n;
    end

    CDL2  = zeros(ML, MI_total(i));
    CDL2n = zeros(ML, MI_total(i));
    CSL2  = zeros(ML, MI_total(i));
    CSL2n = zeros(ML, MI_total(i));
    CDR2  = zeros(ML, MI_total(i));
    CDR2n = zeros(ML, MI_total(i));
    CSR2  = zeros(ML, MI_total(i));
    CSR2n = zeros(ML, MI_total(i));

    M1_index = [0,cumsum(MI(i,:))];
    for ii=1:nop(i)
        [CDL_temp1, CDL_temp1_n] = double_layer_seg_side(k(i), L(i-1), Layer_part{i,ii},   0,  -nei*d);
        CDL2(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDL_temp1;
        CDL2n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDL_temp1_n;
    
        [CSL_temp1, CSL_temp1_n] = single_layer_seg_side(k(i), L(i-1), Layer_part{i,ii},   0,  -nei*d);
        CSL2(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSL_temp1;
        CSL2n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSL_temp1_n;
    
        [CDR_temp1, CDR_temp1_n] = double_layer_seg_side(k(i), R(i-1), Layer_part{i,ii},   0,   nei*d);
        CDR2(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDR_temp1;
        CDR2n(1:ML,M1_index(ii)+1:M1_index(ii+1)) = CDR_temp1_n;
    
        [CSR_temp1, CSR_temp1_n] = single_layer_seg_side(k(i), R(i-1), Layer_part{i,ii},   0,   nei*d);
        CSR2(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSR_temp1;
        CSR2n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSR_temp1_n;
    end
    
    C = [alpha^(-nei-1)*CDR1-alpha^(nei)*CDL1     alpha^(-nei-1)*CSR1-alpha^(nei)*CSL1     alpha^(-nei-1)*CDR2-alpha^(nei)*CDL2     alpha^(-nei-1)*CSR2-alpha^(nei)*CSL2;
         alpha^(-nei-1)*CDR1n-alpha^(nei)*CDL1n   alpha^(-nei-1)*CSR1n-alpha^(nei)*CSL1n   alpha^(-nei-1)*CDR2n-alpha^(nei)*CDL2n   alpha^(-nei-1)*CSR2n-alpha^(nei)*CSL2n];

    
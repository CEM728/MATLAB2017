% Rokhlin-Muller + quasi periodizing scheme + Proxy points for multilayered media
% 
% by Min Hyung Cho and Alex Barnett Oct. 2014
% (mhcho@math.dartmouth.edu, ahb@math.dartmouth.edu)
%
%      incidence wave (incident angle = theta_i)
%            \
%            _\|
%
%               :---------:
%               :  k(1)   :
%   ____________:_________:____________   Layer(1)
%               :  k(2)   :
%   ____________:_________:____________   Layer(2)
%               :   :     :
%               :   :     :
%   ____________:_________:____________   Layer(noi-1)
%               : k(noi)  :
%   ____________:_________:____________   Layer(noi)
%               : k(noi+1):
%               :---------:
%
%
% k(i) = omega*sqrt(epsilon(i)*mu(i))
%
%
%
%
% Note : 
%
%  The code is set to mimic figure 4 in Optics Express paper at lower resolution 
%  (due to randomness, it will look different from actual figure 4 in the paper).
%
%  The variable "nx" can be changed to larger number for high resolution figure (1000 is used for figure 4)
%
%  The exact geometry can be obtained from the Author.
%
%  1. The code requires MPSpack from https://code.google.com/p/mpspack/
%  2. parallel toolbox is used but not required to run the code
%  3. This is not yet user friendly. It will be updated as we make more easy to use.
%
%
%% Initialize parameters
clear; 
warning('off','all');

noi         = 100;                                 %number of interfaces
epsilon     = load('epsilon_random_1_and_2.m');    %load epsilon distribution
epsilon(1)  = 1;                                   %Force the first layer to be vaccum
omega       = 9*pi;                                %omega
theta_i     = -acos(1-2*pi/omega);                 %incident angle at Wood anomaly
figure_on   = 1;                                   %0 for no figure and 1 for figure
nx          = 200;                                  %number of grid points for field plot
field_layer = 10;                                  %plot the field in the top  10 layers
field_layer_bottom = noi-10;                       %plot the field in the last 10 layers
mu = zeros(noi+1,1);
k  = zeros(noi+1,1);
for i=1:noi+1
    mu(i) = 1;                                     %mu = 1 for all layers
    k(i)  = omega*sqrt(epsilon(i)*mu(i));          %k = omega*sqrt(epsilon*mu)
end

k1x = k(1)*cos(theta_i);
k1y = k(1)*sin(theta_i);

%Number of quadrature points on each interface
%Each interface is consist of 'nop' pieces and each piece has MI quadrature points
%For sine interfaces
for i=1:noi
    nop(i)  = 1;
    MI(i,1) = 260;
end
%For triangle interfaces
for i=3:8:noi
     nop(i)  = 2;
     MI(i,1) = 100;
     MI(i,2) = 100;
end
%For rectangle interfaces
for i=4:15:noi
    nop(i)  = 5;
    MI(i,1) = 90;
    MI(i,2) = 90;
    MI(i,3) = 90;
    MI(i,4) = 90;
    MI(i,5) = 90;
end
%Add all number of quadrature points in each piece to make total number of
%quadrature points
MI_total = zeros(noi,1);
for i=1:noi
    for j=1:nop(i)
        MI_total(i) = MI_total(i)+MI(i,j);
    end
end

% Number of points for all other domains
ML          = 120;            %number of points on the left walls
MR          = ML;             %number of points on the right walls
MU          = 60;             %number of points on for Rayleigh-Bloch expansion
MD          = MU;             %number of points on for Rayleigh-Bloch expansion
MP          = 120;            %number of proxy points
BR          = 10;             %number of Rayleigh-Bloch expansion mode
nei         = 1;              %number of neighbor
d           = 1.0;            %periodicity
alpha       = exp(1i*k1x*d);  %quasi-periodic parameter

%% Construct multilayered media
t1 = tic;
%To avoid collision between interfaces
min_thickness = 1;

%Sine interfaces
for i=1:noi
    if i > 1
        height = abs(imag(Layer_part{i-1}.Z(0)))+0.1*rand(1);
        h = -0.3+(0.3+0.3)*rand(1);   
        phase = ((1-0)*rand(1));
        thickness = 0.5*rand(1);
        freq = randi(2,1);

        Layer_part{i,1} = segment(MI(i,1), {@(t) -0.5+t + 1i*h*sin(2*pi*freq*(t-phase))-(height+min_thickness)*1i, @(t) 1+1i*h*2*pi*freq*cos(2*pi*freq*(t-phase)) ,@(t) -1i*h*(2*pi*freq)^2*sin(2*pi*freq*(t-phase))  },'p');
        Layer_part{i,1}.qpblocha = 1/alpha;
            
    else
        h = -0.3+(0.3+0.3)*rand(1);   
        phase = ((1-0)*rand(1));
        thickness = 0.5*rand(1);
        freq = randi(2,1);
        
        Layer_part{i,1} = segment(MI(i,1), {@(t) -0.5+t + 1i*h*sin(2*pi*freq*(t-phase)), @(t) 1+1i*h*2*pi*freq*cos(2*pi*freq*(t-phase)) ,@(t) -1i*h*(2*pi*freq)^2*sin(2*pi*freq*(t-phase))  },'p');
        Layer_part{i,1}.qpblocha = 1/alpha;
    end
end

%Triangle interfaces (starting from 3rd interface and every 8th interface)
for i=3:8:noi
     x1 = -0.25+(0.25+0.25)*rand(1);
     height = abs(imag(Layer_part{i-1}.Z(0)))+0.1*rand(1);
     
     o = []; 
     o.kressq =6; %It needs to be adjusted for high freqeuncy and it will determine convergence
     
     Layer_part{i,1} = segment(MI(i,1), [-0.5-(height+min_thickness)*1i       x1+0.4i-(height+min_thickness)*1i], 'pc', o);
     Layer_part{i,2} = segment(MI(i,2), [x1+0.4i-(height+min_thickness)*1i    0.5-(height+min_thickness)*1i], 'pc', o);
end

%Rectangle interfaces (starting from 4th interface and every 15th interface)
for i=4:15:noi
    x1 = -0.25+(-0.125+0.25)*rand(1);
    x2 =  0.125+(0.25-0.125)*rand(1);
    height = abs(imag(Layer_part{i-1}.Z(0)))+0.1*rand(1);
    h = (0.3+0.1*rand(1))*1i;
    
    o = []; 
    o.kressq =5; %It needs to be adjusted for high freqeuncy and it will determine convergence
    
    Layer_part{i,1} = segment(MI(i,1), [-0.5-(height+min_thickness)*1i       x1-(height+min_thickness)*1i], 'pc',o);
    Layer_part{i,2} = segment(MI(i,2), [x1-(height+min_thickness)*1i         x1+h-(height+min_thickness)*1i], 'pc',o);
    Layer_part{i,3} = segment(MI(i,3), [x1+h-(height+min_thickness)*1i       x2+h-(height+min_thickness)*1i], 'pc',o);
    Layer_part{i,4} = segment(MI(i,4), [x2+h-(height+min_thickness)*1i       x2-(height+min_thickness)*1i], 'pc',o);
    Layer_part{i,5} = segment(MI(i,5), [x2-(height+min_thickness)*1i         0.5-(height+min_thickness)*1i], 'pc', o);
end

%Up and Down Layers to enforece radiation condition via Rayleigh-Bloch expansion
RB_height       = min_thickness;       %height of the layer for Rayleigh-Bloch expansion
RB_height_down  = min_thickness;       %height of the layer for Rayleigh-Bloch expansion

height1 = imag(Layer_part{1,1}.Z(0));
height2 = imag(Layer_part{noi,1}.Z(0));
U_RB = segment(MU, {@(t) -0.5+(1-t) + height1*1i+ RB_height*1i, @(t) -1 ,@(t) 0 },'p');
U_RB.nx = 1i*ones(MU,1);

D_RB = segment(MD, {@(t) -0.5+(1-t) + height2*1i- RB_height_down*1i, @(t) -1 ,@(t) 0 },'p');
D_RB.nx = 1i*ones(MD,1);


%Left and Right boundaries of the unit cell
height1 = imag(Layer_part{1,1}.Z(0));
height2 = imag(Layer_part{noi,1}.Z(0));

L_Top_RB = segment(ML, {@(t) -0.5 + height1*1i+1i*RB_height*t, @(t) 1i*RB_height ,@(t) 0 },'g');
L_Top_RB.nx = ones(ML,1);

R_Top_RB = segment(ML, {@(t) -0.5+d + height1*1i+1i*RB_height*t, @(t) 1i*RB_height ,@(t) 0 },'g');
R_Top_RB.nx = ones(ML,1);

L_Down_RB = segment(ML, {@(t) -0.5-1i*1*(-height2+RB_height_down-RB_height_down*t), @(t) 1i*RB_height_down ,@(t) 0 },'g');
L_Down_RB.nx = ones(ML,1);

R_Down_RB = segment(ML, {@(t) -0.5+d-1i*1*(-height2+RB_height_down-RB_height_down*t), @(t) 1i*RB_height_down ,@(t) 0 },'g');
R_Down_RB.nx = ones(ML,1);

for i=1:noi-1
    height1 = imag(Layer_part{i,1}.Z(0));
    height2 = imag(Layer_part{i+1,1}.Z(0));
    
    L(i) = segment(ML, {@(t) -0.5+1i*(height1*(t)+height2*(1-t)), @(t) 1i*(height1-height2),@(t) 0 },'g');
    L(i).nx = ones(ML,1);

    R(i) = segment(ML, {@(t) -0.5+d+1i*(height1*(t)+height2*(1-t)), @(t) 1i*(height1-height2) ,@(t) 0 },'g');
    R(i).nx = ones(ML,1);
end

%Proxy circles
MFS_radius = 2; %radius of proxy circle
center = imag(Layer_part{1,1}.Z(0))*1i;
MFS(1) = segment(MP, {@(t) -0.5+d/2+MFS_radius*exp(2i*pi*t)+center @(t) 2i*pi*MFS_radius*exp(2i*pi*t)   @(t) 2i*pi*2i*pi*MFS_radius*exp(2i*pi*t)}, 'p' );
for i=2:noi
    center = (imag(Layer_part{i-1,1}.Z(0))+imag(Layer_part{i,1}.Z(0)))*0.5i;
    MFS(i) = segment(MP, {@(t) -0.5+d/2+MFS_radius*exp(2i*pi*t)+center @(t) 2i*pi*MFS_radius*exp(2i*pi*t)   @(t) 2i*pi*2i*pi*MFS_radius*exp(2i*pi*t)}, 'p' );
end
center = imag(Layer_part{noi,1}.Z(0))*1i;
MFS(noi+1) = segment(MP, {@(t) -0.5+d/2+MFS_radius*exp(2i*pi*t)+center @(t) 2i*pi*MFS_radius*exp(2i*pi*t)   @(t) 2i*pi*2i*pi*MFS_radius*exp(2i*pi*t)}, 'p' );

t2=toc(t1);
fprintf('Constructing domain done\n');


%% Evaluate matrix components
%% Construct A matrix for n interfaces using Rokhlin-Muller formulation
t3 = tic;
%For the 1st layer
AD1  = zeros(MI_total(1), MI_total(1));
AD1n = zeros(MI_total(1), MI_total(1));
AS1  = zeros(MI_total(1), MI_total(1));
AS1n = zeros(MI_total(1), MI_total(1));

AD2  = zeros(MI_total(1), MI_total(1));
AD2n = zeros(MI_total(1), MI_total(1));
AS2  = zeros(MI_total(1), MI_total(1));
AS2n = zeros(MI_total(1), MI_total(1));

AD12  = zeros(MI_total(1),MI_total(2));
AS12  = zeros(MI_total(1),MI_total(2));
AD12n = zeros(MI_total(1),MI_total(2));
AS12n = zeros(MI_total(1),MI_total(2));

AD11  = zeros(MI_total(1),MI_total(1));
AS11  = zeros(MI_total(1),MI_total(1));
AD11n = zeros(MI_total(1),MI_total(1));
AS11n = zeros(MI_total(1),MI_total(1));

AD_off   = zeros(MI_total(1), MI_total(2));
AD_off_n = zeros(MI_total(1), MI_total(2));
AS_off   = zeros(MI_total(1), MI_total(2));
AS_off_n = zeros(MI_total(1), MI_total(2));

M1_index = [0,cumsum(MI(1,:))];
M2_index = [0,cumsum(MI(2,:))];
for j=-nei:nei
    a = j*d;
    %diagonal block term
    for ii=1:nop(1)
        for jj=1:nop(1)
            [AD_temp, AD_temp_n] = double_layer_seg(k(1), Layer_part{1,ii}, Layer_part{1,jj}, 1, a);
            [AS_temp, AS_temp_n] = single_layer_seg(k(1), Layer_part{1,ii}, Layer_part{1,jj}, 1, a);
                
            AD1(M1_index(ii)+1:M1_index(ii+1),   M1_index(jj)+1:M1_index(jj+1))=AD_temp;
            AD1n(M1_index(ii)+1:M1_index(ii+1),  M1_index(jj)+1:M1_index(jj+1))=AD_temp_n;
            AS1(M1_index(ii)+1:M1_index(ii+1),   M1_index(jj)+1:M1_index(jj+1))=AS_temp;
            AS1n(M1_index(ii)+1:M1_index(ii+1),  M1_index(jj)+1:M1_index(jj+1))=AS_temp_n;
        end
    end

    for ii=1:nop(1)
        for jj=1:nop(1)
            [AD_temp, AD_temp_n] = double_layer_seg(k(2), Layer_part{1,ii}, Layer_part{1,jj}, 1, a);
            [AS_temp, AS_temp_n] = single_layer_seg(k(2), Layer_part{1,ii}, Layer_part{1,jj}, 1, a);
    
            AD2(M1_index(ii)+1:M1_index(ii+1),   M1_index(jj)+1:M1_index(jj+1))=AD_temp;
            AD2n(M1_index(ii)+1:M1_index(ii+1),  M1_index(jj)+1:M1_index(jj+1))=AD_temp_n;
            AS2(M1_index(ii)+1:M1_index(ii+1),   M1_index(jj)+1:M1_index(jj+1))=AS_temp;
            AS2n(M1_index(ii)+1:M1_index(ii+1),  M1_index(jj)+1:M1_index(jj+1))=AS_temp_n;
        end
    end

    AD11  = AD11 +(alpha^j)*(AD1-AD2);
    AS11  = AS11 +(alpha^j)*(AS1-AS2);   
    AD11n = AD11n+(alpha^j)*(AD1n-AD2n);
    AS11n = AS11n+(alpha^j)*(AS1n-AS2n);

    %Off-diagonal blocks
    for ii=1:nop(1)
        for jj=1:nop(2)
            [AD_temp, AD_temp_n] = double_layer_seg(k(2), Layer_part{1,ii}, Layer_part{2,jj}, 1, a);
            [AS_temp, AS_temp_n] = single_layer_seg(k(2), Layer_part{1,ii}, Layer_part{2,jj}, 1, a);
            
            AD_off(M1_index(ii)+1:M1_index(ii+1),   M2_index(jj)+1:M2_index(jj+1))=AD_temp;
            AD_off_n(M1_index(ii)+1:M1_index(ii+1), M2_index(jj)+1:M2_index(jj+1))=AD_temp_n;
            AS_off(M1_index(ii)+1:M1_index(ii+1),   M2_index(jj)+1:M2_index(jj+1))=AS_temp;
            AS_off_n(M1_index(ii)+1:M1_index(ii+1), M2_index(jj)+1:M2_index(jj+1))=AS_temp_n;
            
        end
    end
    AD12  = AD12 +(alpha^j)*AD_off;
    AS12  = AS12 +(alpha^j)*AS_off;
    AD12n = AD12n+(alpha^j)*AD_off_n;
    AS12n = AS12n+(alpha^j)*AS_off_n;
end

A{1,1} = [-eye(MI_total(1),MI_total(1))+AD11   AS11;
                                        AD11n  eye(MI_total(1),MI_total(1))+AS11n];
A{1,2} = [ -AD12  -AS12 ;
           -AD12n -AS12n];


       
%For the last layer       
AD11  = zeros(MI_total(noi),MI_total(noi));
AS11  = zeros(MI_total(noi),MI_total(noi));
AD11n = zeros(MI_total(noi),MI_total(noi));
AS11n = zeros(MI_total(noi),MI_total(noi));

AD1 = zeros(MI_total(noi), MI_total(noi));
AD1n = zeros(MI_total(noi), MI_total(noi));
AS1 = zeros(MI_total(noi), MI_total(noi));
AS1n = zeros(MI_total(noi), MI_total(noi));

AD2 = zeros(MI_total(noi), MI_total(noi));
AD2n = zeros(MI_total(noi), MI_total(noi));
AS2 = zeros(MI_total(noi), MI_total(noi));
AS2n = zeros(MI_total(noi), MI_total(noi));

AD_off   = zeros(MI_total(noi), MI_total(noi-1));
AD_off_n = zeros(MI_total(noi), MI_total(noi-1));
AS_off   = zeros(MI_total(noi), MI_total(noi-1));
AS_off_n = zeros(MI_total(noi), MI_total(noi-1));

AD21  = zeros(MI_total(noi),MI_total(noi-1));
AS21  = zeros(MI_total(noi),MI_total(noi-1));
AD21n = zeros(MI_total(noi),MI_total(noi-1));
AS21n = zeros(MI_total(noi),MI_total(noi-1));


M1_index = [0,cumsum(MI(noi,:))];
M2_index = [0,cumsum(MI(noi-1,:))];
for j=-nei:nei
    a = j*d;
    %Diagonal block
    for ii=1:nop(noi)
        for jj=1:nop(noi)
            [AD_temp, AD_temp_n] = double_layer_seg(k(noi), Layer_part{noi,ii}, Layer_part{noi,jj}, 1, a);
            [AS_temp, AS_temp_n] = single_layer_seg(k(noi), Layer_part{noi,ii}, Layer_part{noi,jj}, 1, a);
            
            AD1(M1_index(ii)+1:M1_index(ii+1),   M1_index(jj)+1:M1_index(jj+1))=AD_temp;
            AD1n(M1_index(ii)+1:M1_index(ii+1),  M1_index(jj)+1:M1_index(jj+1))=AD_temp_n;
            AS1(M1_index(ii)+1:M1_index(ii+1),   M1_index(jj)+1:M1_index(jj+1))=AS_temp;
            AS1n(M1_index(ii)+1:M1_index(ii+1),  M1_index(jj)+1:M1_index(jj+1))=AS_temp_n;
            
            [AD_temp, AD_temp_n] = double_layer_seg(k(noi+1), Layer_part{noi,ii}, Layer_part{noi,jj}, 1, a);
            [AS_temp, AS_temp_n] = single_layer_seg(k(noi+1), Layer_part{noi,ii}, Layer_part{noi,jj}, 1, a);
            
            AD2(M1_index(ii)+1:M1_index(ii+1),   M1_index(jj)+1:M1_index(jj+1))=AD_temp;
            AD2n(M1_index(ii)+1:M1_index(ii+1),  M1_index(jj)+1:M1_index(jj+1))=AD_temp_n;
            AS2(M1_index(ii)+1:M1_index(ii+1),   M1_index(jj)+1:M1_index(jj+1))=AS_temp;
            AS2n(M1_index(ii)+1:M1_index(ii+1),  M1_index(jj)+1:M1_index(jj+1))=AS_temp_n;
        end
    end

    AD11  = AD11 +(alpha^j)*(AD1-AD2);
    AS11  = AS11 +(alpha^j)*(AS1-AS2);   
    AD11n = AD11n+(alpha^j)*(AD1n-AD2n);
    AS11n = AS11n+(alpha^j)*(AS1n-AS2n);
    
    %Off-diagonal block
    for ii=1:nop(noi)
        for jj=1:nop(noi-1)
            [AD_temp, AD_temp_n] = double_layer_seg(k(noi), Layer_part{noi,ii}, Layer_part{noi-1,jj}, 1, a);
            [AS_temp, AS_temp_n] = single_layer_seg(k(noi), Layer_part{noi,ii}, Layer_part{noi-1,jj}, 1, a);
            
            AD_off(M1_index(ii)+1:M1_index(ii+1),   M2_index(jj)+1:M2_index(jj+1))=AD_temp;
            AD_off_n(M1_index(ii)+1:M1_index(ii+1), M2_index(jj)+1:M2_index(jj+1))=AD_temp_n;
            AS_off(M1_index(ii)+1:M1_index(ii+1),   M2_index(jj)+1:M2_index(jj+1))=AS_temp;
            AS_off_n(M1_index(ii)+1:M1_index(ii+1), M2_index(jj)+1:M2_index(jj+1))=AS_temp_n;
        end
    end
    AD21  = AD21 +(alpha^j)*AD_off;
    AS21  = AS21 +(alpha^j)*AS_off;
    AD21n = AD21n+(alpha^j)*AD_off_n;
    AS21n = AS21n+(alpha^j)*AS_off_n;
end
A{noi, noi}   =  [-eye(MI_total(noi),MI_total(noi))+AD11                                      AS11;... 
                                                    AD11n    eye(MI_total(noi),MI_total(noi))+AS11n];       
A{noi, noi-1} =  [AD21     AS21 ;... 
                  AD21n    AS21n];         

              
%For the layers between 2nd and (noi-1)th layers (parallel for)
parfor i=2:noi-1
    [AD11, AS11, AD11n, AS11n, AD12, AS12, AD12n, AS12n, AD21, AS21, AD21n, AS21n] = construct_A(i, k, Layer_part, nei, alpha, d, MI, MI_total, nop);
     Atemp{i} = [AD21     AS21  -eye(MI_total(i),MI_total(i))+AD11                                 AS11      -AD12   -AS12;... 
                 AD21n    AS21n                               AD11n   eye(MI_total(i),MI_total(i))+AS11n     -AD12n  -AS12n];         
end         

for i=2:noi-1
    A{i,i-1} = Atemp{i}(1:2*MI_total(i),1:2*MI_total(i-1));
    A{i,i}   = Atemp{i}(1:2*MI_total(i),2*MI_total(i-1)+1:2*MI_total(i-1)+2*MI_total(i));
    A{i,i+1} = Atemp{i}(1:2*MI_total(i),2*MI_total(i-1)+2*MI_total(i)+1:2*MI_total(i-1)+2*MI_total(i)+2*MI_total(i+1));
end

t4 = toc(t3);
fprintf('Constructing A done\n');

%% Construct B matrix (layer due to proxy circles)
t5 = tic;
B = cell(noi+1);

%For the 1st layer
B1  = []; 
B1n = [];
for ii=1:nop(1)
    [B_temp, B_temp_n]   = double_layer_seg(k(1),  Layer_part{1,ii},   MFS(1), 0, 0);
    [BS_temp, BS_temp_n] = single_layer_seg(k(1),  Layer_part{1,ii},   MFS(1), 0, 0);
    B1  = [B1;B_temp+1i*k(1)*BS_temp];
    B1n = [B1n;B_temp_n+1i*k(1)*BS_temp_n];
end
B{1} = [B1; B1n];

%For the last layer
B1  = []; 
B1n = [];
for ii=1:nop(noi)
    [B_temp, B_temp_n]   = double_layer_seg(k(noi+1),  Layer_part{noi,ii},   MFS(noi+1), 0, 0);
    [BS_temp, BS_temp_n] = single_layer_seg(k(noi+1),  Layer_part{noi,ii},   MFS(noi+1), 0, 0);
    B1  = [B1;B_temp+1i*k(noi+1)*BS_temp];
    B1n = [B1n;B_temp_n+1i*k(noi+1)*BS_temp_n];
end
B{noi+1}  = [-B1;-B1n];

for i=2:noi
    [B1, B1n] = construct_B1(i, k, Layer_part, MFS, nop);
    [B2, B2n] = construct_B2(i, k, Layer_part, MFS, nop);
    B{i}      = [-B1;-B1n;B2;B2n];
end

t6 = toc(t5);
fprintf('Constructing B done\n');


%% Construct C matrix (left-right matching quasi-periodic condition)
t7 = tic;
C = cell(noi+1);

%For the 1st layer
CDL1  = zeros(ML, MI_total(1));
CDL1n = zeros(ML, MI_total(1));
CSL1  = zeros(ML, MI_total(1));
CSL1n = zeros(ML, MI_total(1));
CDR1  = zeros(ML, MI_total(1));
CDR1n = zeros(ML, MI_total(1));
CSR1  = zeros(ML, MI_total(1));
CSR1n = zeros(ML, MI_total(1));  

M1_index = [0,cumsum(MI(1,:))];
for ii=1:nop(1)
    [CDL_temp1, CDL_temp1_n] = double_layer_seg_side(k(1), L_Top_RB, Layer_part{1,ii}, 0, -nei*d);
    CDL1(1:ML, M1_index(ii)+1:M1_index(ii+1))  = CDL_temp1;
    CDL1n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDL_temp1_n;
    
    [CSL_temp1, CSL_temp1_n] = single_layer_seg_side(k(1), L_Top_RB, Layer_part{1,ii}, 0, -nei*d);
    CSL1(1:ML, M1_index(ii)+1:M1_index(ii+1))  = CSL_temp1;
    CSL1n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSL_temp1_n;

    [CDR_temp1, CDR_temp1_n] = double_layer_seg_side(k(1), R_Top_RB, Layer_part{1,ii}, 0, nei*d);
    CDR1(1:ML, M1_index(ii)+1:M1_index(ii+1))  = CDR_temp1;
    CDR1n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDR_temp1_n;
    
    [CSR_temp1, CSR_temp1_n] = single_layer_seg_side(k(1), R_Top_RB, Layer_part{1,ii}, 0, nei*d);
    CSR1(1:ML, M1_index(ii)+1:M1_index(ii+1))  = CSR_temp1;
    CSR1n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSR_temp1_n;
end

C{1} = [alpha^(-nei-1)*CDR1-alpha^(nei)*CDL1     alpha^(-nei-1)*CSR1-alpha^(nei)*CSL1   ;
        alpha^(-nei-1)*CDR1n-alpha^(nei)*CDL1n   alpha^(-nei-1)*CSR1n-alpha^(nei)*CSL1n];


%For the last layer
CDL1  = zeros(ML, MI_total(noi));
CDL1n = zeros(ML, MI_total(noi));
CSL1  = zeros(ML, MI_total(noi));
CSL1n = zeros(ML, MI_total(noi));
CDR1  = zeros(ML, MI_total(noi));
CDR1n = zeros(ML, MI_total(noi));
CSR1  = zeros(ML, MI_total(noi));
CSR1n = zeros(ML, MI_total(noi));  

M1_index = [0,cumsum(MI(noi,:))];

for ii=1:nop(noi)    
    [CDL_temp1, CDL_temp1_n] = double_layer_seg_side(k(noi+1), L_Down_RB, Layer_part{noi,ii}, 0,  -nei*d);
    CDL1(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDL_temp1;
    CDL1n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDL_temp1_n;
    
    [CSL_temp1, CSL_temp1_n] = single_layer_seg_side(k(noi+1), L_Down_RB, Layer_part{noi,ii}, 0,  -nei*d);
    CSL1(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSL_temp1;
    CSL1n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSL_temp1_n;
    
    [CDR_temp1, CDR_temp1_n] = double_layer_seg_side(k(noi+1), R_Down_RB, Layer_part{noi,ii}, 0,   nei*d);
    CDR1(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDR_temp1;
    CDR1n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CDR_temp1_n;
    
    [CSR_temp1, CSR_temp1_n] = single_layer_seg_side(k(noi+1), R_Down_RB, Layer_part{noi,ii}, 0,   nei*d);
    CSR1(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSR_temp1;
    CSR1n(1:ML, M1_index(ii)+1:M1_index(ii+1)) = CSR_temp1_n;
end
C{noi+1} = [alpha^(-nei-1)*CDR1-alpha^(nei)*CDL1     alpha^(-nei-1)*CSR1-alpha^(nei)*CSL1;
            alpha^(-nei-1)*CDR1n-alpha^(nei)*CDL1n   alpha^(-nei-1)*CSR1n-alpha^(nei)*CSL1n];

%For the 2nd to noi layers (parfor is used)
parfor i=2:noi
    C{i} = construct_C(i, k, Layer_part, L, R, nop, MI, MI_total, ML, nei, d, alpha);
end

t8 = toc(t7);
fprintf('Constructing C done\n');

     
%% Construct Q matrix (left and right due to proxy)
t9 = tic;
Q = cell(noi+1);

%For the 1st layer
[Q1L1, Q1L1n]   = double_layer_seg(k(1), L_Top_RB, MFS(1), 0, 0);
[QS1L1, QS1L1n] = single_layer_seg(k(1), L_Top_RB, MFS(1), 0, 0);

[Q1R1, Q1R1n]   = double_layer_seg(k(1), R_Top_RB, MFS(1), 0, 0);
[QS1R1, QS1R1n] = single_layer_seg(k(1), R_Top_RB, MFS(1), 0, 0);

Q{1}     = [(Q1R1+1i*k(1)*QS1R1)/alpha-(Q1L1+1i*k(1)*QS1L1);
            (Q1R1n+1i*k(1)*QS1R1n)/alpha-(Q1L1n+1i*k(1)*QS1L1n)];

%For the last layer
[Q1L1, Q1L1n]   = double_layer_seg(k(noi+1), L_Down_RB, MFS(noi+1), 0, 0);
[QS1L1, QS1L1n] = single_layer_seg(k(noi+1), L_Down_RB, MFS(noi+1), 0, 0);

[Q1R1, Q1R1n]   = double_layer_seg(k(noi+1), R_Down_RB, MFS(noi+1), 0, 0);
[QS1R1, QS1R1n] = single_layer_seg(k(noi+1), R_Down_RB, MFS(noi+1), 0, 0);

Q{noi+1} = [(Q1R1+1i*k(noi+1)*QS1R1)/alpha-(Q1L1+1i*k(noi+1)*QS1L1);
            (Q1R1n+1i*k(noi+1)*QS1R1n)/alpha-(Q1L1n+1i*k(noi+1)*QS1L1n)];
        
%For all other layers
for i=2:noi 
    [Q1L1, Q1L1n]   = double_layer_seg(k(i), L(i-1), MFS(i), 0, 0);
    [QS1L1, QS1L1n] = single_layer_seg(k(i), L(i-1), MFS(i), 0, 0);
    
    [Q1R1, Q1R1n]   = double_layer_seg(k(i), R(i-1), MFS(i), 0, 0);
    [QS1R1, QS1R1n] = single_layer_seg(k(i), R(i-1), MFS(i), 0, 0);
    
    Q{i}=[(Q1R1+1i*k(i)*QS1R1)/alpha-(Q1L1+1i*k(i)*QS1L1);
          (Q1R1n+1i*k(i)*QS1R1n)/alpha-(Q1L1n+1i*k(i)*QS1L1n)];
end

t10 = toc(t9);
fprintf('Constructing Q done\n');



%% Construct Z matrix (Up and Down Rayleigh-Bloch due to Top and bottom interfaces)
t11 = tic;
ZU11 = zeros(MU,MI_total(1));
ZU12 = zeros(MU,MI_total(1));
ZU21 = zeros(MU,MI_total(1));
ZU22 = zeros(MU,MI_total(1));

ZD11 = zeros(MD,MI_total(noi));
ZD12 = zeros(MD,MI_total(noi));
ZD21 = zeros(MD,MI_total(noi));
ZD22 = zeros(MD,MI_total(noi));

M1_index = [0,cumsum(MI(1,:))];
M2_index = [0,cumsum(MI(noi,:))];


for i=-nei:nei
    a = i*d;
    for ii=1:nop(1)
        [Z_temp1, Z_temp1_n] = double_layer_seg(k(1),  U_RB, Layer_part{1,ii},   0, a);
        ZU_D1(1:MU,  M1_index(ii)+1:M1_index(ii+1)) = Z_temp1;
        ZU_D1n(1:MU, M1_index(ii)+1:M1_index(ii+1)) = Z_temp1_n;
        
        [Z_temp1, Z_temp1_n] = single_layer_seg(k(1),  U_RB, Layer_part{1,ii},   0, a);
        ZU_S1(1:MU,  M1_index(ii)+1:M1_index(ii+1)) = Z_temp1;
        ZU_S1n(1:MU, M1_index(ii)+1:M1_index(ii+1)) = Z_temp1_n;
    end
    
    for ii=1:nop(noi)
        [Z_temp1, Z_temp1_n] = double_layer_seg(k(noi+1), D_RB, Layer_part{noi,ii}, 0, a);
        ZD_D2(1:MD,  M2_index(ii)+1:M2_index(ii+1)) = Z_temp1;
        ZD_D2n(1:MD, M2_index(ii)+1:M2_index(ii+1)) = Z_temp1_n;
        
        [Z_temp1, Z_temp1_n] = single_layer_seg(k(noi+1), D_RB, Layer_part{noi,ii}, 0, a);
        ZD_S2(1:MU,  M2_index(ii)+1:M2_index(ii+1)) = Z_temp1;
        ZD_S2n(1:MU, M2_index(ii)+1:M2_index(ii+1)) = Z_temp1_n;
    end
    
    ZU11 = ZU11+(alpha^i)*ZU_D1;
    ZU21 = ZU21+(alpha^i)*ZU_D1n;
    ZU12 = ZU12+(alpha^i)*ZU_S1;   
    ZU22 = ZU22+(alpha^i)*ZU_S1n;
    
    ZD11 = ZD11+(alpha^i)*ZD_D2;
    ZD21 = ZD21+(alpha^i)*ZD_D2n;
    ZD12 = ZD12+(alpha^i)*ZD_S2;
    ZD22 = ZD22+(alpha^i)*ZD_S2n;
end
ZU = [ZU11  ZU12 ;
      ZU21  ZU22];

ZD = [ZD11  ZD12 ;
      ZD21  ZD22];
  
t12 = toc(t11);
fprintf('Constructing Z done\n');



%% Construct V matrix (Up and Down Rayleigh-Bloch due to proxy)      
t13 = tic;
[VU1, VU1n] = double_layer_seg(k(1),     U_RB, MFS(1),     0, 0);
[VSU1, VSU1n] = single_layer_seg(k(1),     U_RB, MFS(1),     0, 0);

[VD2, VD2n] = double_layer_seg(k(noi+1), D_RB, MFS(noi+1), 0, 0);
[VSD2, VSD2n] = single_layer_seg(k(noi+1), D_RB, MFS(noi+1), 0, 0);

VU = [VU1+1i*k(1)*VSU1;VU1n+1i*k(1)*VSU1n];
VD = [VD2+1i*k(noi+1)*VSD2;VD2n+1i*k(noi+1)*VSD2n];  

t14 = toc(t13);
fprintf('Constructing V done\n');

%% Construct W matrix (Rayleigh-Bloch coefficients at top and bottom interfaces)
t15 = tic;
[WU1, WU1n] = rayleigh(U_RB, k1x, k(1),     d, MU, BR);
[WD1, WD1n] = rayleigh(D_RB, k1x, k(noi+1), d, MD, BR);

WU = [-WU1;-WU1n];
WD = [-WD1; WD1n];

t16 = toc(t15);
fprintf('Constructing W done\n');

%Rearrange matrix for Schur Complement
B{1}     = [B{1}      zeros(2*MI_total(1),   (2*BR+1))];
B{noi+1} = [B{noi+1}  zeros(2*MI_total(noi), (2*BR+1))];

C{1}     = [C{1};    ZU];
C{noi+1} = [C{noi+1};ZD];

Q{1}     = [Q{1}      zeros(2*ML, 2*BR+1);
            VU        WU];
Q{noi+1} = [Q{noi+1}  zeros(2*ML, 2*BR+1);
            VD        WD];
        

t17=tic;
%Construct right-hand side
f=[];fn=[];
for ii=1:nop
    f_temp = -exp(1i*(k1x*real(Layer_part{1,ii}.x)+k1y*imag(Layer_part{1,ii}.x)));  % -u_inc  on top interfaces
    f_temp_n = 1i*f_temp.*(k1x*real(Layer_part{1,ii}.nx)+k1y*imag(Layer_part{1,ii}.nx));   % -un_inc on top interfaces
    f=[f;f_temp];
    fn = [fn;f_temp_n];
end

rhs_schur = cell(noi);
rhs_schur{1} = [f;fn];

t18 = toc(t17);
t_matrix_filling = toc(t3);

fprintf('Constructing RHS done\n');

%% Solve the linear system with Schur Complement
%Schur complement
t19 = tic;
%Evalute Q^{-dag}*C for all Q
QdagC = cell(noi+1);
for i=1:noi+1
    QdagC{i} = Q{i}\C{i};
end

%Evalute schur complement A-B*Q^{-dag}*C and overwrite A
A{1,1} = A{1,1}-B{1}*QdagC{1};
for i = 1:noi-1
     A_temp = [A{i,i}  A{i,i+1}; A{i+1,i}  A{i+1,i+1}]-B{i+1}*QdagC{i+1};
     A{i,i}     = A_temp(1:2*MI_total(i),                                 1:2*MI_total(i));
     A{i,i+1}   = A_temp(1:2*MI_total(i),                                 2*MI_total(i)+1:2*MI_total(i)+2*MI_total(i+1));
     A{i+1,i}   = A_temp(2*MI_total(i)+1:2*MI_total(i)+2*MI_total(i+1),   1:2*MI_total(i));
     A{i+1,i+1} = A_temp(2*MI_total(i)+1:2*MI_total(i)+2*MI_total(i+1),   2*MI_total(i)+1:2*MI_total(i)+2*MI_total(i+1));
end
A{noi,noi} = A{noi,noi}-B{noi+1}*QdagC{noi+1};

t_schur_complement = toc(t19);
fprintf('Schur Complement done\n');


%Solve tri-diagonal block matrix
%Eliminate (i+1,1) block using (i,i) block and construct upper triangular block matrix
t21 = tic;
for i=1:noi-2
    A_inv = inv(A{i,i});
    A{i+1,i+1} = A{i+1,i}*(A_inv*A{i,i+1})-A{i+1,i+1};
    A{i+1,i+2} = -A{i+1, i+2};
    rhs_schur{i+1} = A{i+1,i}*(A_inv*rhs_schur{i});
end

A_inv=inv(A{noi-1,noi-1});
A{noi,noi}    = A{noi,noi-1}*(A_inv*A{noi-1,noi})-A{noi,noi};
rhs_schur{noi} = A{noi,noi-1}*(A_inv*rhs_schur{noi-1});

%Solve upper triangular block matrix with back substitution to find eta
eta{noi} = A{noi,noi}\rhs_schur{noi};
for i=noi-1:-1:1
    eta{i} = A{i,i}\(rhs_schur{i}-A{i,i+1}*eta{i+1});
end

%Seperate eta into tau and sigma
tau_schur   = cell(noi, MI_total(i));
sigma_schur = cell(noi, MI_total(i));
for i=1:noi
    tau_schur{i} = eta{i}(1:MI_total(i));
    sigma_schur{i} = eta{i}(MI_total(i)+1:2*MI_total(i));
end

%Use eta to find Rayleigh Block coefficients (au_schur, ad_schur) and Proxy solution (b_schur) saved in x
x = cell(noi+1);
x{1} = -QdagC{1}*[tau_schur{1};sigma_schur{1}];
for i=2:noi
    x{i} = -QdagC{i}*[tau_schur{i-1};sigma_schur{i-1};tau_schur{i};sigma_schur{i};];
end
x{noi+1} = -QdagC{noi+1}*[tau_schur{noi};sigma_schur{noi}];

t_block_solve = toc(t21);
t_matrix_solve = toc(t19);

fprintf('Block Matrix Solve done\n');

%Seprate b_schur, au_schur, ad_schur from x
b_schur    = cell(noi+1);
b_schur{1} = x{1}(1:MP);
au_schur = x{1}(MP+1:MP+2*BR+1);
for i=2:noi
    b_schur{i} = x{i};
end
b_schur{noi+1} = x{noi+1}(1:MP);
ad_schur = x{noi+1}(MP+1:MP+2*BR+1);

%Evalute flux
[flux_up, flux_down] = compute_flux(au_schur, ad_schur, k1x, k(1), k(noi+1), d, BR);
total_flux_up        = sum(real(flux_up));
total_flux_down      = sum(real(flux_down));
flux_error           = abs((total_flux_up+total_flux_down-abs(k1y))/abs(k1y));
fprintf('k1y            = %15.12f\n',   abs(k1y));
fprintf('Total flux     = %15.12f\n\n', total_flux_up+total_flux_down);
fprintf('Flux error     = %g\n\n',      flux_error);
fprintf('Reflection     = %15.12f\n',   (total_flux_up/abs(k1y)));
fprintf('Transmission   = %15.12f\n\n', (total_flux_down/abs(k1y)));

fprintf('Order                Up                                     Down\n' );
for i=-5:5
    fprintf('%3d       %e+1i*(%e)       %e+1i*(%e) \n', i, real(flux_up(i+BR+1)), imag(flux_up(i+BR+1)),   real(flux_down(i+BR+1)), imag(flux_down(i+BR+1)) );
end

%Compute total memory used
ram = 0;
total_variables = whos;
for i=1:length(total_variables)
    ram = ram+total_variables(i).bytes;
end
ram = (ram/1024)/1024;


%% Compute solution at a test point
%Evaluate the solution at testpoint in the 1st layer
testpoint.x = 0.4+1i*0.6;
testpoint.nx = [];
utest = (double_layer_seg(k(1), testpoint, MFS(1), 0, 0)...
        +1i*k(1)*single_layer_seg(k(1), testpoint, MFS(1), 0, 0))*b_schur{1};


M1_index = [0,cumsum(MI(1,:))];
for i=-nei:nei 
    a = i*d;
    for jj=1:nop(1)
          utest = utest + (alpha^i)*double_layer_seg(k(1), testpoint, Layer_part{1,jj}, 1, a)*tau_schur{1}(M1_index(jj)+1:M1_index(jj+1))...
                        + (alpha^i)*single_layer_seg(k(1), testpoint, Layer_part{1,jj}, 1, a)*sigma_schur{1}(M1_index(jj)+1:M1_index(jj+1)); 
    end
end

fprintf('\n\nSolution at a test point\n');
fprintf('u(%e, %e) = (%15.12f)+i(%15.12f)\n\n', real(testpoint.x), imag(testpoint.x), real(utest), imag(utest));

fprintf('\nConstructing domain  : %f\n\n', t2);
fprintf('Constructing A       : %f\n', t4);
fprintf('Constructing B       : %f\n', t6);
fprintf('Constructing C       : %f\n', t8);
fprintf('Constructing Q       : %f\n', t10);
fprintf('Constructing Z       : %f\n', t12);
fprintf('Constructing V       : %f\n', t14);
fprintf('Constructing W       : %f\n', t16);
fprintf('Constructing rhs     : %f\n\n', t18);
fprintf('Total matrix filling : %f\n\n', t_matrix_filling);
fprintf('Schur Complement     : %f\n', t_schur_complement);
fprintf('Matrix Solve         : %f\n\n', t_block_solve);
fprintf('Total matrix Solve   : %f\n\n', t_matrix_solve);

fprintf('\n\nTotal memory = %15.12f\n\n',ram);


%% Compute solution for graph
if figure_on    
    t=tic;
    %plotting grid on top layers
    gx1      =  -0.5:1/(nx-1):0.5;
    gy1      =  imag(Layer_part{field_layer,1}.Z(0)) : (imag(U_RB.x(1))-imag(Layer_part{field_layer,1}.Z(0)))/(nx-1):imag(U_RB.x(1));
    [xx, yy]  =  meshgrid(gx1,gy1); 
    zz1      = (xx+1i*yy); 
    clear xx yy
    total_grid.x  = zz1(:); 
    total_grid.nx = [];
    
    %plotting grid on bottom layers
    gx2      =  -0.5:1/(nx-1):0.5;
    gy2      =  imag(D_RB.x(1)) : (imag(Layer_part{field_layer_bottom,1}.Z(0))-imag(D_RB.x(1)))/(nx-1):imag(Layer_part{field_layer_bottom,1}.Z(0));
    [xx, yy]  =  meshgrid(gx2,gy2); 
    zz2      = (xx+1i*yy); 
    clear xx yy
    total_grid_bottom.x  = zz2(:); 
    total_grid_bottom.nx = [];
    
    
    %Determine where the grid point is in which layer
    %Find grid points in the 1st layer
    xv=[];
    yv=[];
    for i=1:nop(1)
        xv = [xv real(Layer_part{1,i}.x)'];
        yv = [yv imag(Layer_part{1,i}.x)'];
    end
    xv = [real(Layer_part{1,1}.Z(0))    xv   real(Layer_part{1,nop(1)}.Z(1))  real(Layer_part{1,nop(1)}.Z(1))                    real(Layer_part{1,1}.Z(0))                        real(Layer_part{1, 1}.Z(0))];
    yv = [imag(Layer_part{1,1}.Z(0))    yv   imag(Layer_part{1,1}.Z(0))  imag(Layer_part{1,1}.Z(0))+RB_height                    imag(Layer_part{1,1}.Z(0))+RB_height         imag(Layer_part{1, 1}.Z(0))];
    
    in{1} = inpolygon(real(total_grid.x),imag(total_grid.x),xv,yv);
    
    %Find grid points located between 2nd and field_layer
    for j=2:field_layer+1
        clear xv yv
        xv1=[];
        yv1=[];
        xv2=[];
        yv2=[];
        for i=1:nop(j)
            xv1 = [xv1 real(Layer_part{j,i}.x)'];
            yv1 = [yv1 imag(Layer_part{j,i}.x)'];
        end
        
        for i=1:nop(j-1)
            xv2 = [xv2 real(Layer_part{j-1,i}.x)'];
            yv2 = [yv2 imag(Layer_part{j-1,i}.x)'];
        end
        xv2 = [xv2 real(Layer_part{j-1,nop(j-1)}.Z(1))];
        yv2 = [yv2 imag(Layer_part{j-1,nop(j-1)}.Z(1))];
        
        xv2 = fliplr(xv2);
        yv2 = fliplr(yv2);

        xv = [real(Layer_part{j,1}.Z(0)) xv1 real(Layer_part{j,nop(j)}.Z(1)) xv2  real(Layer_part{j-1,1}.Z(0))  real(Layer_part{j,1}.Z(0))];
        yv = [imag(Layer_part{j,1}.Z(0)) yv1 imag(Layer_part{j,nop(j)}.Z(1)) yv2  imag(Layer_part{j-1,1}.Z(0))  imag(Layer_part{j,1}.Z(0))];

        in{j} = inpolygon(real(total_grid.x),imag(total_grid.x),xv,yv);
    end
    
    
    for j=field_layer_bottom:noi
        clear xv yv
        xv1=[];
        yv1=[];
        xv2=[];
        yv2=[];
        for i=1:nop(j)
            xv1 = [xv1 real(Layer_part{j,i}.x)'];
            yv1 = [yv1 imag(Layer_part{j,i}.x)'];
        end
        
        for i=1:nop(j-1)
            xv2 = [xv2 real(Layer_part{j-1,i}.x)'];
            yv2 = [yv2 imag(Layer_part{j-1,i}.x)'];
        end
        xv2 = [xv2 real(Layer_part{j-1,nop(j-1)}.Z(1))];
        yv2 = [yv2 imag(Layer_part{j-1,nop(j-1)}.Z(1))];
        
        xv2 = fliplr(xv2);
        yv2 = fliplr(yv2);

        xv = [real(Layer_part{j,1}.Z(0)) xv1 real(Layer_part{j,nop(j)}.Z(1)) xv2  real(Layer_part{j-1,1}.Z(0))  real(Layer_part{j,1}.Z(0))];
        yv = [imag(Layer_part{j,1}.Z(0)) yv1 imag(Layer_part{j,nop(j)}.Z(1)) yv2  imag(Layer_part{j-1,1}.Z(0))  imag(Layer_part{j,1}.Z(0))];

        in{j} = inpolygon(real(total_grid_bottom.x),imag(total_grid_bottom.x),xv,yv);
        
    end
    
    
    %Find grid points located in the last layer
    xv=[];
    yv=[];
    xv2=[];
    yv2=[];
    for i=1:nop(noi)
        xv = [xv real(Layer_part{noi,i}.x)'];
        yv = [yv imag(Layer_part{noi,i}.x)'];
    end
    xv = [real(Layer_part{noi,1}.Z(0)) xv real(Layer_part{noi,nop(noi)}.Z(1))];
    yv = [imag(Layer_part{noi,1}.Z(0)) yv imag(Layer_part{noi,nop(noi)}.Z(1))];
    
    xv = fliplr(xv);
    yv = fliplr(yv);

    xv2 = [real(Layer_part{noi,1}.Z(0))                 real(Layer_part{noi,nop(noi)}.Z(1))               xv       real(Layer_part{noi,1}.Z(0)) ];
    yv2 = [imag(Layer_part{noi,1}.Z(0))-RB_height_down  imag(Layer_part{noi,1}.Z(0))-RB_height_down       yv       imag(Layer_part{noi,1}.Z(0))-RB_height_down  ];

    in{noi+1} = inpolygon(real(total_grid_bottom.x), imag(total_grid_bottom.x),xv2,yv2);
    
    %Compute MFS
    ug = zeros(nx*nx,1);
    for j=1:field_layer+1
        indices       = find(in{j});
        grid_point.x  = total_grid.x(indices);
        grid_point.nx = [];
        ug(indices)   = (double_layer_seg(k(j), grid_point, MFS(j), 0, 0)...
                        +1i*k(j)*single_layer_seg(k(j), grid_point, MFS(j), 0, 0))*b_schur{j};
    end
    
    ug_bottom = zeros(nx*nx,1);
    for j=field_layer_bottom:noi+1
       indices       = find(in{j});
       grid_point.x  = total_grid_bottom.x(indices);
       grid_point.nx = [];
       ug_bottom(indices) = (double_layer_seg(k(j), grid_point, MFS(j), 0, 0)...
                            +1i*k(j)*single_layer_seg(k(j), grid_point, MFS(j), 0, 0))*b_schur{j};
    end

    %First layer
    indices        = find(in{1});
    grid_point.x   = total_grid.x(indices);
    grid_point.nx  = [];
    M1_index = [0,cumsum(MI(1,:))];
    for i=-nei:nei
        a = i*d;
        for j=1:nop(1)
              ug(indices) = ug(indices) + (alpha^i)*double_layer_seg(k(1), grid_point, Layer_part{1,j}, 1, a)*tau_schur{1}(M1_index(j)+1:M1_index(j+1))...
                                        + (alpha^i)*single_layer_seg(k(1), grid_point, Layer_part{1,j}, 1, a)*sigma_schur{1}(M1_index(j)+1:M1_index(j+1)); 
        end
    end
    ug(indices)=ug(indices)+exp(1i*(k1x.*real(grid_point.x)+k1y.*imag(grid_point.x)));
   
    %second to field_layer
    for j=2:field_layer+1
        indices       = find(in{j});
        grid_point.x  = total_grid.x(indices);
        grid_point.nx = [];
        M1_index = [0,cumsum(MI(j-1,:))];
        M2_index = [0,cumsum(MI(j,:))];
        for i=-nei:nei
            a = i*d;
            for jj=1:nop(j-1)
                ug(indices) = ug(indices) + (alpha^i)*double_layer_seg(k(j), grid_point, Layer_part{j-1,jj}, 1, a)*tau_schur{j-1}(M1_index(jj)+1:M1_index(jj+1))...
                                          + (alpha^i)*single_layer_seg(k(j), grid_point, Layer_part{j-1,jj}, 1, a)*sigma_schur{j-1}(M1_index(jj)+1:M1_index(jj+1)); 
            end
            for jj=1:nop(j)
                ug(indices) = ug(indices) + (alpha^i)*double_layer_seg(k(j), grid_point, Layer_part{j,jj}, 1, a)*tau_schur{j}(M2_index(jj)+1:M2_index(jj+1))...
                                          + (alpha^i)*single_layer_seg(k(j), grid_point, Layer_part{j,jj}, 1, a)*sigma_schur{j}(M2_index(jj)+1:M2_index(jj+1)); 
            end
        end
    end
    ug_total = reshape(ug, size(zz1));
    
    for j=field_layer_bottom:noi
        indices       = find(in{j});
        grid_point.x  = total_grid_bottom.x(indices);
        grid_point.nx = [];
        M1_index = [0,cumsum(MI(j-1,:))];
        M2_index = [0,cumsum(MI(j,:))];
        for i=-nei:nei
            a = i*d;
            for jj=1:nop(j-1)
                ug_bottom(indices) = ug_bottom(indices) + (alpha^i)*double_layer_seg(k(j), grid_point, Layer_part{j-1,jj}, 1, a)*tau_schur{j-1}(M1_index(jj)+1:M1_index(jj+1))...
                                          + (alpha^i)*single_layer_seg(k(j), grid_point, Layer_part{j-1,jj}, 1, a)*sigma_schur{j-1}(M1_index(jj)+1:M1_index(jj+1)); 
            end
            for jj=1:nop(j)
                ug_bottom(indices) = ug_bottom(indices) + (alpha^i)*double_layer_seg(k(j), grid_point, Layer_part{j,jj}, 1, a)*tau_schur{j}(M2_index(jj)+1:M2_index(jj+1))...
                                          + (alpha^i)*single_layer_seg(k(j), grid_point, Layer_part{j,jj}, 1, a)*sigma_schur{j}(M2_index(jj)+1:M2_index(jj+1)); 
            end
        end
    end
    
    
    %Last layer
    indices       = find(in{noi+1});
    grid_point.x  = total_grid_bottom.x(indices);
    grid_point.nx = [];
    M1_index = [0,cumsum(MI(noi,:))];
    for ii=-nei:nei
        a = ii*d;
        for jj=1:nop(noi)
            ug_bottom(indices) = ug_bottom(indices) + (alpha^ii)*double_layer_seg(k(noi+1), grid_point, Layer_part{noi,jj}, 1, a)*tau_schur{noi}(M1_index(jj)+1:M1_index(jj+1))...
                                                    + (alpha^ii)*single_layer_seg(k(noi+1), grid_point, Layer_part{noi,jj}, 1, a)*sigma_schur{noi}(M1_index(jj)+1:M1_index(jj+1)); 
        end
    end
    
    ug_total_bottom = reshape(ug_bottom, size(zz2));


    %Evalute total field above and below the artificial layer
    [field_up, field_down, gx_up, gy_up, gx_down, gy_down] = compute_field_n_layer(au_schur, ad_schur, k1x, k1y, k(1), k(noi+1), d, imag(U_RB.x(1)), imag(D_RB.x(1)), BR);

    str1 = num2str(noi);
    str2 = 'data';
    str3 = [str2 str1];
    save(str3)
    
    t_compute_field=toc(t);
    fprintf('Computing Fields   : %f\n\n', t_compute_field);

    
    %% Plot domain, total field, tau, and sigma
    %Plot Domain
    p_copy = 2;   %number of periodic copy to plot
    figure;
    subplot(1,3,1)
    hold on
    for i=1:10
        for ii=-p_copy:p_copy
            for j=1:nop(i)
                plot(real(Layer_part{i,j}.x)+ii, imag(Layer_part{i,j}.x));
            end
        end
    end
    title('First 10 layers');
    axis equal;
    axis tight;
    axis([-p_copy-0.5 p_copy+0.5 imag(Layer_part{field_layer,1}.x(1)) imag(U_RB.x(1))+1]);
   
    %plot the total field in top layers
    %figure;
    subplot(1,3,2)
    hold on
    for i=-p_copy:p_copy
        imagesc(gx1+i,gy1,real(ug_total*(alpha^i)))
        imagesc(gx_up+i,gy_up, real(field_up*(alpha^i)));
    end

    for i=1:field_layer
        for j=-p_copy:p_copy
            for ii=1:nop(i)
                plot(real(Layer_part{i,ii}.x)+j, imag(Layer_part{i,ii}.x),'k','LineWidth',2);
            end
        end
    end
    caxis([-1.7,1.7]);
    axis equal
    axis tight;
    axis([-p_copy-0.5 p_copy+0.5 imag(Layer_part{field_layer,1}.x(1)) imag(U_RB.x(1))+1]);
    colorbar
    title('Re(Total field) in Top');

    %plot the total field in bottom layers
    subplot(1,3,3)
    hold on
    for i=-p_copy:p_copy
        imagesc(gx2+i,gy2,real(ug_total_bottom*(alpha^i)))
        imagesc(gx_down+i,gy_down, real(field_down*(alpha^i)));
    end

    hold on
    for i=field_layer_bottom:noi
        for j=-p_copy:p_copy
            for ii=1:nop(i)
                plot(real(Layer_part{i,ii}.x)+j, imag(Layer_part{i,ii}.x),'k','LineWidth',2);
            end
        end
    end
    caxis([-1.0,1.0]);
    axis equal
    axis tight;
    axis([-p_copy-0.5 p_copy+0.5 imag(D_RB.x(1))+0.3 imag(Layer_part{field_layer_bottom,1}.x(1))]);
    colorbar
    title('Re(Total field) in Bottom');
end

function [A, An] = single_layer_seg(k, t, s, isinterface, a) 
if nargin < 5, 
    a = 0; 
end
wantder = nargout>1;

s_copy = copy(s);
s_copy.x = s_copy.x+a;

if isinterface==1
    opt.quad = 'a'; 
    opt.ord  = 16;      
else
    opt.quad = 'm'; 
end


if t.x(1)==s_copy.x(1)   %Self interaction
    %for fast hankel function (hank106)
%     N = numel(s_copy.x); 
%     M = numel(s_copy.x);                                               % # target pts
%     opt.displ = repmat(s_copy.x, [1 N]) - repmat(s_copy.x.', [M 1]);   % C-# displ mat
%     opt.rdist = abs(opt.displ);
%     [opt.Sker opt.Dker_noang] = utils.greengardrokhlinhank106(k*opt.rdist);
%     opt.Sker = (1i/4) * opt.Sker;
%     opt.Dker_noang = (1i*k/4) * opt.Dker_noang;
    
    % quadrature must be Alpert
    B  = layerpot.S(k, s_copy, [], opt);  
    if wantder
        opt.derivSLP = 1;  % option for D^T
        Bn = layerpot.D(k, s_copy, [], opt);  
    end
else
    %for fast hankel function (hank106)
%     N = numel(s_copy.x); 
%     M = numel(t.x);                                               % # target pts
%     opt.displ = repmat(t.x, [1 N]) - repmat(s_copy.x.', [M 1]);   % C-# displ mat
%     opt.rdist = abs(opt.displ);
%     [opt.Sker opt.Dker_noang] = utils.greengardrokhlinhank106(k*opt.rdist);
%     opt.Sker = (1i/4) * opt.Sker;
%     opt.Dker_noang = (1i*k/4) * opt.Dker_noang;

    B  = layerpot.S(k, s_copy, t, opt);  
    if wantder, 
        opt.derivSLP = 1;  % option for D^T
        Bn = layerpot.D(k, s_copy, t, opt);   
    end
end

A=B;
if wantder               
 An = Bn;
end

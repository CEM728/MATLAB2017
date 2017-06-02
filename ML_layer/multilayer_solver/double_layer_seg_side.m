function [A, An] = double_layer_seg_side(k, t, s, isinterface, a) 
if nargin < 5, 
    a = 0; 
end
wantder = nargout>1;

t_copy = copy(t);
t_copy.x = t_copy.x+a;

if isinterface==1
    opt.quad = 'm'; 
else
    opt.quad = 'm'; 
end


if t_copy.x(1)==s.x(1) %self interaction
    %for fast hankel function (hank106)
%     N = numel(t_copy.x); 
%     M = numel(t_copy.x);                                             % # target pts
%     opt.displ = repmat(t_copy.x, [1 N]) - repmat(t_copy.x.', [M 1]); % C-# displ mat
%     opt.rdist = abs(opt.displ);
%     [opt.Sker opt.Dker_noang] = utils.greengardrokhlinhank106(k*opt.rdist);
%     opt.Sker = (1i/4) * opt.Sker;
%     opt.Dker_noang = (1i*k/4) * opt.Dker_noang;
%     
    B  = layerpot.D(k, s, [],opt);  
    if wantder
        Bn = layerpot.T(k, s, [],opt);  
    end
else
    %for fast hankel function (hank106)
%     N = numel(s.x); 
%     M = numel(t_copy.x);                                             % # target pts
%     opt.displ = repmat(t_copy.x, [1 N]) - repmat(s.x.', [M 1]);      % C-# displ mat
%     opt.rdist = abs(opt.displ);
%     [opt.Sker opt.Dker_noang] = utils.greengardrokhlinhank106(k*opt.rdist);
%     opt.Sker = (1i/4) * opt.Sker;
%     opt.Dker_noang = (1i*k/4) * opt.Dker_noang;
    
    B  = layerpot.D(k, s, t_copy,opt);  
    if wantder, 
        Bn = layerpot.T(k, s, t_copy,opt);  
    end
end

A=B;
if wantder                                           
 An = Bn;
end



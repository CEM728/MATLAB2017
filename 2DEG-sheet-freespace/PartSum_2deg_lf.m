function s = PartSum_2deg_lf(a,b,eh, e2h, n, eta)
% This function performs the summation of the terms generated by the
% product of weights and nodes

ekh = eh;
s = Term_2deg_lf(a,b,ekh, eta);
for k = 2 : n
    ekh = ekh * e2h;
    s = s + Term_2deg_lf(a,b,ekh,eta);
end
end

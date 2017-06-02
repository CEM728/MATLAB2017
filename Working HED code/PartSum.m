function s = PartSum(a,b,eh, e2h, n, eta)
% This function performs the summation of the terms generated by the
% product of weights and nodes

ekh = eh;
s = Term(a,b,ekh, eta);
for k = 2 : n
    ekh = ekh * e2h;
    s = s + Term(a,b,ekh,eta);
end
end

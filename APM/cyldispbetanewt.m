function newt = cyldispbetanewt(betasqr, cyl)

% evaluates ratio of dispersion relation and its derivative 
% for beta squared as variable

% set beta if not provided
if nargin == 1
  betasqr = cyl.beta^2;
end

% interior exterior arguments
kperpi = sqrt(cyl.k.^2.*cyl.epi.*cyl.mu - betasqr);
% adjust branch cut and choose branch
kperp = -sqrtbr(cyl.k.^2.*cyl.ep.*cyl.mu - betasqr, -pi/2);

ka = kperp.*cyl.a;
kia = kperpi.*cyl.a;

% scaled bessel functions
bessj = besselj(cyl.orders, kia, 1);
bessjm = besselj(cyl.orders-1, kia, 1);
bessjp = besselj(cyl.orders+1, kia, 1);

bessh = besselk(cyl.orders, -1i.*ka, 1)./i.^(cyl.orders+1);
besshm = besselk(cyl.orders-1, -1i.*ka, 1)./i.^(cyl.orders);
besshp = besselk(cyl.orders+1, -1i.*ka, 1)./i.^(cyl.orders+2);

% equation ratios
jrat = 0.5.*(bessjm-bessjp)./bessj./kia;
hrat = 0.5.*(besshm-besshp)./bessh./ka;

% derivatives
djratdk = (bessjp.*bessjm - bessj.*(bessjm-bessjp)./kia - bessj.^2)./kperpi./bessj.^2;
dhratdk = (besshp.*besshm - bessh.*(besshm-besshp)./ka - bessh.^2)./kperp./bessh.^2;
dkdbeta = -0.5./kperp;
dkidbeta = -0.5./kperpi;

% determinant derivative
detd = cyl.k.^2.*(djratdk.*dkidbeta-dhratdk.*dkdbeta).*(cyl.epi.*jrat-cyl.ep.*hrat) + cyl.k.^2.*(jrat-hrat).*(cyl.epi.*djratdk.*dkidbeta-cyl.ep.*dhratdk.*dkdbeta) - cyl.orders.^2.*(kia.^-2-ka.^-2).^2.*(1+2.*cyl.a.^2.*betasqr.*(kia.^-2+ka.^-2));

% determinant
detm = cyl.k.^2.*(jrat-hrat).*(cyl.epi.*jrat-cyl.ep.*hrat) - cyl.orders.^2*betasqr.*(kia.^-2-ka.^-2).^2;

newt = detm./detd;

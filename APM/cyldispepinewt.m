function newt = cyldispepinewt(epi, cyl)

% evaluates ratio of dispersion relation and its derivative 
% for cylinder epsilon as variable

% set epi if not provided
if nargin == 1
  epi = cyl.epi;
end

% interior exterior arguments
% branch cut irrelevant
kperp = sqrt(cyl.k.^2.*cyl.ep.*cyl.mu - cyl.beta.^2);
kperpi = sqrt(cyl.k.^2.*epi.*cyl.mu - cyl.beta.^2);

ka = kperp.*cyl.a;
kia = kperpi.*cyl.a;

% scaled bessel functions
bessj = besselj(cyl.orders, kia, 1);
bessjm = besselj(cyl.orders-1, kia, 1);
bessjp = besselj(cyl.orders+1, kia, 1);

bessh = besselh(cyl.orders, 1, ka, 1);
besshm = besselh(cyl.orders-1, 1, ka, 1);
besshp = besselh(cyl.orders+1, 1, ka, 1);

% equation ratios
jrat = 0.5.*(bessjm-bessjp)./bessj./kia;
hrat = 0.5.*(besshm-besshp)./bessh./ka;

% derivatives
djratdk = (bessjp.*bessjm - bessj.*(bessjm-bessjp)./kia - bessj.^2)./kperpi./bessj.^2;
dkdepi = cyl.k.^2./2./kperpi;

% determinant derivative
detd = cyl.k.^2.*djratdk.*dkdepi.*(epi.*jrat-cyl.ep.*hrat) + cyl.k.^2.*(jrat-hrat).*(jrat + epi.*djratdk.*dkdepi) + 4.*cyl.orders.^2*cyl.beta.^2./kperpi./kia.^2.*(kia.^-2-ka.^-2).*dkdepi;

% determinant
detm = cyl.k.^2.*(jrat-hrat).*(epi.*jrat-cyl.ep.*hrat) - cyl.orders.^2*cyl.beta.^2.*(kia.^-2-ka.^-2).^2;

newt = detm./detd;

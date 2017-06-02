function [roots, epis, cen, rad] = disprootsepi(cyl, nroots)

% Finds roots of the cylinder dispersion relation using argument principle method 
%
% Usage
% roots = disprootsepi(cyl, nroots)
% [roots, sings, cen, rad] = disprootsepi(cyl, nroots)
%
% Inputs
% cyl - structure containing cylinder parameters
% nroots - number of roots desired, beginning from fundamental
%
% Outputs
% roots - locations of roots
% sings - (optional) locations of singularities of dispersion relation
% cen - (optional) centers of contours
% rad - (optianal) radii of contours

% number of contours
% each singularity known to have 2 roots nearby
ncntrs = ceil(nroots/2);

% get zero of J, corresponds to kperpi*a
jzero = besselzero(cyl.orders, max(4, ncntrs+1));
% calculate corresponding epsilon
epis = ((jzero./cyl.a).^2 + cyl.beta.^2)./cyl.k.^2./cyl.mu;

% find plasmonic root
trustr = realmax;
kperp = sqrt(cyl.k.^2.*cyl.ep.*cyl.mu - cyl.beta.^2);
if real(kperp) < eps(1)
  proot = newton(@cyldispepinewt, -1, trustr, cyl);
elseif cyl.orders ~= 0
  proot = newton(@cyldispepinewt, -1-0.1i, trustr, cyl);
else
  proot = [];
end

% contour method for positive plasmonic root
if real(proot) > 0
  cen = 0; rad = epis(1)/2;
  obl = 1; phi = 0;
  sings = cyl.beta^2/cyl.k^2; proots = [];
  roots = apm(@cyldispepinewt, cen, rad, obl, phi, sings, proots, 1e-5, cyl);
  proot = newton(@cyldispepinewt, min(roots), trustr, cyl);
end

% integration contours: avoid plasmonic root if necessary
[rad, cen] = contours(epis, proot);

% recalculate number of contours necessary
ncntrs = ceil((nroots-length(proot))/2);
droots = zeros([1 2*ncntrs]);

% root search by contours centered on singularities
for k = 1:ncntrs
  droots([2*k-1 2*k]) = search(cyl, epis, [proot droots(1:2*k-2)], cen(k), rad(k));
end

% check for reality
roots = [proot droots];
if all(abs(imag(roots)./real(roots)) < eps)
  roots = real(roots);
end

% sort roots
[~, order] = sort(real(roots));
roots = roots(order);

% exclude unwanted roots
roots = roots(1:nroots);

function [rad, cen] = contours(sings, proot)
% lower boundary is 3/4 distance to next lower singularity 
dz = diff(sings);
llim = [0 sings(1:end-2)+dz(1:end-1)./4];

% upper boundary is 3/4 distance to next higher singularity
dz = diff(sings);
ulim = sings(1:end-1)+3/4.*dz;

% lowest boundary is extrapolation of lowest singularities
%ising = interp1(1:4, sings(1:4), 0, 'pchip');
ising = interp1(1:2, sings(1:2), 0, 'pchip');

% extend boundary if lowest boundary if close to unity
if ising > 50
  llim(1) = ising;
else
  if isempty(proot) || abs(proot) > 20;
    % contour excludes plasmonic root
    llim(1) = -10;
  else
    llim(1) = -abs(proot)-10;
  end
end

% calculate center and radius of contours
cen = (ulim+llim)./2;
rad = (ulim-llim)./2;

function roots = search(cyl, sing, proot, cen, rad)
% course contour parameters
crstol = 1e-5;
conbuffer = 2;

% fine contour parameters
fintol = 1e-10;
sdistlim = 1e3*crstol;
fconsize = 1.5;

% newton-raphson parameters
trustr = 5;

% detect if known roots are close to contour
% use only after first contour
if length(proot) > 1
  nearc = abs(abs(cen-proot)-rad) < conbuffer;

  % adjust contour
  if any(nearc)
    % existing contour
    llim = cen-rad; ulim = cen+rad;

    % new contour bisects root and next lowest contour
    buffroot = min(proot(nearc));
    lroot = max(sing(sing < buffroot));
    llim = (buffroot+lroot)/2;
    cen = (ulim+llim)/2; rad = (ulim-llim)/2;
  end
end

% detect singularities in contour
incs = abs(cen - sing) < rad;
% bessel singularities known to be double poles
sings = [sing(incs) sing(incs)];

% include zero alpha singularity if it exists in large contour
if cyl.orders ~= 0 && abs(cen - cyl.beta^2/cyl.k^2) < rad
  sings = [cyl.beta^2/cyl.k^2 sings];
end

% detect previous roots in large contour
incp = abs(cen - proot) < rad;
proots = proot(incp);

% large course contour to include at least two roots
phi = 0; obl = 1;
roots = apm(@cyldispepinewt, cen, rad, obl, phi, sings, proots, crstol, cyl);
roots = sort(roots);

% test if roots are close to central singularity
sdist = abs(sings(end)-roots)/abs(sings(end));
csing = sdist < sdistlim;

% small fine contour for roots close to singularity
if any(csing)
  % center: average of points, radius: larger than furtherst point
  fpts = [roots(csing); sings(end)];
  fcen = mean(fpts);
  frad = max([fconsize.*abs(fcen-fpts); 1e-3*rad]);
  fphi = 1; fobl = 1;
  fsings = [sings(end) sings(end)]; fzeros = [];
  roots(csing) = apm(@cyldispepinewt, fcen, frad, fobl, fphi, fsings, fzeros, fintol, cyl);
end

% polish roots using Newton-Raphson
% only treat two roots per contour
for k = 1:2
  roots(k) = newton(@cyldispepinewt, roots(k), trustr, cyl);
end

% return two roots
roots = roots(1:2);

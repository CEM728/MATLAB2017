function [roots, betas, cen, rad] = disprootsbeta(cyl, nroots)

% Finds roots of the cylinder dispersion relation using argument principle method 
%
% Usage
% roots = disprootsbeta(cyl, nroots)
% [roots, sings, cen, rad] = disprootsbeta(cyl, nroots)
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
ncntrs = ceil((nroots+1)/2);

% singularity at zero missing for m=0
zsing = cyl.orders ~=0;

% get zero of J, corresponds to kperpi*a
jzero = besselzero(cyl.orders, max(4, ncntrs+1));
% include zero at zero argument
if zsing; jzero = [0 jzero]; end;
% calculate corresponding betas
betas = cyl.k.^2.*cyl.epi.*cyl.mu-(jzero./cyl.a).^2;

% integration contours
branchpt = cyl.k.^2.*cyl.ep.*cyl.mu;
[rad, cen] = contours(betas, branchpt);

% all poles are double except at zero
betas = [betas betas(1+zsing:end)];

% root search by contours centered on singularities
count = 0;
roots = zeros([1 2*ncntrs]);
for k = 1:ncntrs
  newroots = search(cyl, betas, roots(1:count), cen(k), rad(k));
  roots(count+1:count+length(newroots)) = newroots;
  count = count + length(newroots);
end

% check for reality
if all(abs(imag(roots)./real(roots)) < eps)
  roots = real(roots);
end

% sort roots
[~, order] = sort(real(roots(1:count)),2,'descend');
roots = roots(order);

% exclude unwanted roots
% not all roots may have been found
nroots = min([nroots count]);
roots = roots(1:nroots);

function [rad, cen, flag] = contours(sings, branchpt)
% branch cut buffer
buff = 1e-7;

% radii and centers are 3/4 distance to adjacent singularities
cen = 0.25.*sings(2:end-1) + 0.375.*sings(1:end-2) + 0.375.*sings(3:end);
rad = abs(0.375.*sings(3:end) - 0.375.*sings(1:end-2));

% first radius is 3/4 distance to second singularity
cen = [sings(1) cen];
rad = [abs(0.75*(sings(2)-sings(1))) rad];

% branch point collision detection
flag = false;
for s = 1:length(rad)
if abs(branchpt-cen(s)) < rad(s)
  % place new center halfway between branchpt and old boundary
  % maintain imaginary offset
  if real(sings(s)) > branchpt
    cen(s) = real(branchpt+cen(s)+rad(s))/2 + 1i*imag(cen(s));
  else
    cen(s) = real(branchpt+cen(s)-rad(s))/2 + 1i*imag(cen(s));
  end
  % radius avoids branch point
  rad(s) = abs(branchpt-cen(s)) - buff;
end
end

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

% detect singularities in contour
incs = abs(cen - sing) < rad;
sings = sing(incs);

% detect previous roots in large contour
incp = abs(cen - proot) < rad;
proots = proot(incp);

% large coarse contour to include at least two roots
phi = 0; obl = 1;

roots = apm(@cyldispbetanewt, cen, rad, obl, phi, sings, proots, crstol, cyl);
roots = sort(roots);

% small fine contour not implemented because roots usually not close to singularity

% polish roots using Newton-Raphson
% only treat at most two largest roots per contour
inds = max([1 length(roots)-1]):length(roots);
for k = inds
  roots(k) = newton(@cyldispbetanewt, roots(k), trustr, cyl);
end

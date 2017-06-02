function sroots = contour(fn, c, r, o, phi, zsings, zroots, tol, varargin)

% Implementation of argument principle method for root search with deflation of singularities and known roots
% Quadrature: equidistant trapezoidal rule (number of evaluations points increased iteratively, switching to adaptive Simpson's rule
% Polishing: Newton's method
%
% Input parameters
% fn - function handle to transcendental equation (derivative/function value)
% c - center of contour
% r - radius of contour
% o - oblateness of contour
% phi - rotation of contour
% zsings - singularities within contour
% zroots - known roots within contour
% tol - tolerance threshold (determins root accuracy)
% varargin - input arguments for function handle

% count number of singularities and known roots in contour
nsings = length(zsings); kroots = length(zroots);

% initial and max number of points for trapezoidal
ipts = 32; maxpts = 512;

% set up evaluation points
t(1:maxpts+1) = linspace(0,1,maxpts+1);
t(maxpts+1) = [];

% create oblate contour
expt = exp(1i*phi).*(cos(2.*pi.*t) + 1i.*o.*sin(2.*pi.*t));
dexpt = exp(1i.*phi).*(o.*cos(2.*pi.*t) + 1i.*sin(2.*pi.*t));
z = c + r.*expt;

% define additional and cumulative points
tot = 1:maxpts/ipts:maxpts;
add = tot;

% pre-allocate variables
fnval = zeros([1 maxpts]);
flag = false;

% trapezoidal rule until count of roots reaches tolerance
for k = 1:log2(maxpts/ipts)+1
  % evaluate ratio of function and derivative
  fnval(add) = dexpt(add)./fn(z(add), varargin{:});
  
  % evaluate error based on count of roots
  cnt = r*mean(fnval(tot));
  err = abs(cnt - round(cnt));

  % end loop conditions
  if err < tol
    flag = true;
    break
  else % intercalate points
    tot = 1:maxpts/ipts/2^k:maxpts;
    add = 1+maxpts/ipts/2^k:maxpts/ipts/2^(k-1):maxpts;
  end
end

% function handle for Simpson's rule
function grandval = intgrand(t)
  expt = exp(1i*phi).*(cos(2.*pi.*t) + 1i.*o.*sin(2.*pi.*t));
  dexpt = exp(1i.*phi).*(o.*cos(2.*pi.*t) + 1i.*sin(2.*pi.*t));
  z = c + r.*expt;
  
  % evaluate ratio of function and derivative
  grandval = dexpt./fn(z, varargin{:});
end

% use adaptive simpson if maxpts reached
if ~flag
  stol = max([1e-2*tol 1e-12]); % different stopping criterion for Simpson
  [int, fnvals, segs] = adaptsimp(@intgrand, 0, 1, [fnval fnval(1)], stol);
  cnt = r*int;
end

% print number of function evaluations
%if flag
%  npts = length(tot)
%else
%  npts = length([fnvals{:}])
%end

% roots to evaluate
nroots = round(cnt) + nsings - kroots;

% evaluate the sum of powers
s = 1:nroots;
if flag
  % sp = r/mean(expt(tot).^s.*fnval(tot))
  sp = r.*mean(bsxfun(@times, bsxfun(@power, expt(tot).', s), fnval(tot).'));
else
  sp = r.*adaptfft(nroots, fnvals, segs);
end

% sum of powers of singularities
if nsings > 0
  % normalized coordinates
  us = (zsings-c)./r;
  % sps = 2*sum(us.^s);
  sps = sum(bsxfun(@power, us.', s), 1);
  
  % deflate singularities
  sp = sp + sps;
end

% sum of powers of known roots
if kroots > 0
  % normalized coordinates
  ur = (zroots-c)./r;
  % spr = 2*sum(us.^s);
  spr = sum(bsxfun(@power, ur.', s), 1);
  
  % deflate roots
  sp = sp - spr;
end

% recursive evaluation symmetric polynomials
% index begins with zeroth order
ep = zeros([1 nroots]);
ep(1) = 1;
for k = 1:nroots
  inds = 1:k;
  ep(k+1) = 1/k * sum((-1).^(inds-1).*ep(k-inds+1).*sp(inds));
end

% define polynomial and find roots
poly = (-1).^(0:nroots) .* ep;
uroots = roots(poly);

% shift roots from unit circle
sroots = sort(c + r.*uroots);
end

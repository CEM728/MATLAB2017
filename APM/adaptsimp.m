function [Q, y, idxQ] = adaptsimp(fn, a, b, iny, tol, varargin)

% default values
dtol = 1e-8; initparts = 8; maxloops = 50;

% test input args
if nargin < 5 || isempty(tol)
  tol = dtol;
end

% initial values, termination constant (incl. quadr factor 15)
Q = 0; 
tol = 15 * tol; 
poleWarning = 0;

% initial evaluation points or use supplied values 
if nargin >= 4 && ~isempty(iny)
  t = linspace(a, b, length(iny));
  y{1} = iny;
else
  t = linspace(a, b, 4*initparts+1);
  y{1} = fn(t, varargin{:}); 
end

% initial interval
dx = t(5) - t(1);

% poles at initial points
if ~isempty(find(~isfinite(max(abs(y{1}))))), poleWarning = 1; end

% initialize interval boundaries and values
A = t(1:4:end-1); yA = y{1}(1:4:end-1);
B = t(5:4:end);   yB = y{1}(5:4:end);
M = t(3:4:end-1); yM = y{1}(3:4:end-1);
L = t(2:4:end-1); yL = y{1}(2:4:end-1);
R = t(4:4:end-1); yR = y{1}(4:4:end-1);

for k = 1:maxloops
  % Simpson formulas (on course and fine grid)
  Q1 = dx/6 .* (yA + 4.*yM + yB);
  Q2 = dx/12 .* (yA + 4.*yL + 2.*yM + 4.*yR + yB);

  % difference of Simpson formulas
  diffQ = Q2 - Q1; diffQ(find(isnan(diffQ))) = 0;

  % intervals which do not fulfill termination criterion
  idx = find(abs(diffQ) > tol);

  % intervals fulfill termination criterion
  idxQ{k} = setdiff(1:length(A), idx);

  % termination
  if isempty(idx)
    Q = Q + sum(Q2+diffQ/15);
    break
  elseif k == maxloops
    Q = Q + sum(Q2+diffQ/15);
    warning('Matlab:adaptsimp:MaxResolution',...
    ['Maximum number of iterations, ' num2str(maxloops) ', reached.'])
    break
  end

  % update quadrature value
  Q = Q + sum(Q2(idxQ{k}) + diffQ(idxQ{k})/15);

  % update A, B, and M
  tk = zeros(1, 2*length(idx));
  tk(1:2:end) = A(idx); tk(2:2:end) = M(idx); A = tk;
  tk(1:2:end) = M(idx); tk(2:2:end) = B(idx); B = tk;
  tk(1:2:end) = L(idx); tk(2:2:end) = R(idx); M = tk;

  % update yA, yB, and yM
  yk = zeros(1, 2*length(idx));
  yk(1:2:end) = yA(idx); yk(2:2:end) = yM(idx); yA = yk;
  yk(1:2:end) = yM(idx); yk(2:2:end) = yB(idx); yB = yk;
  yk(1:2:end) = yL(idx); yk(2:2:end) = yR(idx); yM = yk;

  % update interval length
  dx = dx/2;

  % update L and R by interval bisection
  L = (A + dx/4); R = (A + 3/4*dx);

  % calculate yD and yE
  t = [L R];
  y{k+1} = fn(t); 

  % poles at new points
  if ~isempty(find(~isfinite(max(abs(y{k+1}))))), poleWarning = 1; end

  % assign new values to yD and yE
  yL = y{k+1}(1:end/2); yR = y{k+1}(end/2+1: end);
end

%evals = length([y{:}])

if poleWarning
  warning('Matlab:adaptsimp:PoleDetection',...
    'A detection of a pole; singularity likely.')
end

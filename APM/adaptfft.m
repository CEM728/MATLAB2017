function Q = adaptfft(m, y, idxQ)

% evaluates moments using preexisting function evaluations from adaptsimp

% initial values
Q = 0; p = (1:m).';

% integration interval must be [0 1]: extract spacing
t = linspace(0, 1, length(y{1}));
dx = t(5) - t(1);

% initialize interval boundaries and values
A = t(1:4:end-1); yA = y{1}(1:4:end-1);
B = t(5:4:end);   yB = y{1}(5:4:end);
M = t(3:4:end-1); yM = y{1}(3:4:end-1);
L = t(2:4:end-1); yL = y{1}(2:4:end-1);
R = t(4:4:end-1); yR = y{1}(4:4:end-1);

for k = 1:length(y)
  if ~isempty(idxQ{k})
    % evaluate exponents
    exptA = bsxfun(@power, exp(2i.*pi.*A(idxQ{k})), p);
    exptB = bsxfun(@power, exp(2i.*pi.*B(idxQ{k})), p);
    exptM = bsxfun(@power, exp(2i.*pi.*M(idxQ{k})), p);
    exptL = bsxfun(@power, exp(2i.*pi.*L(idxQ{k})), p);
    exptR = bsxfun(@power, exp(2i.*pi.*R(idxQ{k})), p);

    % Simpson formulas including exponents (on course and fine grid)
    Q1 = dx/6 .* (...
    bsxfun(@times, exptA, yA(idxQ{k}))... 
    + 4.*bsxfun(@times, exptM, yM(idxQ{k}))...
    + bsxfun(@times, exptB, yB(idxQ{k})));

    Q2 = dx/12 .* (...
    bsxfun(@times, exptA, yA(idxQ{k}))...
    + bsxfun(@times, exptL, 4.*yL(idxQ{k}))...
    + bsxfun(@times, exptM, 2.*yM(idxQ{k}))... 
    + bsxfun(@times, exptR, 4.*yR(idxQ{k}))... 
    + bsxfun(@times, exptB, yB(idxQ{k})));
  
    % difference of Simpson formulas
    diffQ = Q2 - Q1;
  
    % update quadrature value
    Q = Q + sum(Q2 + diffQ/15, 2);
  end
 
  if k == length(y), break, end

  % extract intervals that did not fulfill termination criterion
  idx = setdiff(1:length(A), idxQ{k});

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

  % assign new values to yD and yE
  yL = y{k+1}(1:end/2); yR = y{k+1}(end/2+1: end);
end

% row vector output
Q = Q.';

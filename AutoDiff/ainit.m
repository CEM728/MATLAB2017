function varargout = ainit(varargin)

%AINIT Initialize audi variable.
%   x = AINIT(v) is an audi variable with value v and derivative 1 
%   to be used as argument of univariate functions, as in y = f(x). 
%   v is a double array of arbitrary dimension.
%
%   x = AINIT(v,k) is an audi variable with derivatives up to order k, 
%   where all derivatives of order 2 or higher are set to 0. 
%
%   [x1,...,xn] = AINIT(v1,...,vn,[k=1]) is a set of audi variables to
%   be used as arguments of functions depending on n variables, as in 
%   y = f(x1,...,xn).
%
%   Example: <a href="matlab: ahelp(1)">Functions, derivatives, and Newton's methods</a>
%   See also: aeval, adim, aord, asize

n = nargout;
if nargin == n
  k = 1;
else
  k = varargin{end};
end
varargout = cell(1,n);
for i = 1:n
  if ~isequal(size(varargin{1}),size(varargin{i}))
    error('All variables must have equal size.')
  end
  varargout{i} = audi(varargin{i},i,n,k);
end
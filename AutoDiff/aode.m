function varargout = aode(f,x0,y0,k,varargin)

%AODE Taylor polynomial of solution of ODE.
%   p = AODE(f,x0,y0,k) is the degree k Taylor polynomial of the solution
%   of the scalar differential equation y'(x) = f(x,y) with y(x0) = y0.
%   f = @(x,y)... is an anonymous function, and x0,y0 are double scalars.
%   
%   [p1,...,pm] = AODE(f,x0,y0,k) are the degree k Taylor poynomials of the  
%   solution y = [y1;...;ym] of the system y'(x) = f(x,y) in R^m with y(x0) = y0.
%   f = @(x,y1,...,ym) = [...;...;...] is an anonymous function returning vectors,
%   x0 is a double scalar, and y0 is a double vector.
%
%   [p1,...,pm] = AODE(f,x0,y0,k,c1,c2,...) passes additional parameters to
%   the function f = @(x,y1,...,ym,c1,c2,...)...
%
%   Caution: AODE is manipulating the function string of f, thereby substituting 
%   names of arguments by the internal standard names x,y1,...,yn. Occasionally, 
%   this causes errors if a name of an argument also appears in the function string 
%   as part of other expressions. For example, p = aode(@(t,y)atan(t*y),0,1,3) 
%   fails  because "atan(t,y)" is replaced by "axan(x,y1)". It is always save to use 
%   the standard names x,y1,...,yn. Unfortunately, there are no save names for the 
%   parameters, as they are internally replaced by the given values. Just take
%   care to choose names which do not appear elsewhere in the function string.
% 
%   Example: <a href="matlab: ahelp(3)">The mathematical pendulum</a>
%   See also: ataylor, polyval

f = func2str(f);                            % function string to be parsed
f = strrep(f,'y','Y');
n = numel(y0);
c = find(f==',',1);
x = f(3:c-1);
f = f(c+1:end);
for i = 1:n
  c = find((f==',') | (f==')'),1);
  yi= f(1:c-1);
  f = f(c+1:end);
  f = strrep(f,yi,['y' num2str(i)]);        % rename arguments as y1,...,ym
end
for i = 1:numel(varargin)
  c = find((f==',') | (f==')'),1);
  z = f(1:c-1);
  f = f(c+1:end);
  f = strrep(f,z,num2str(varargin{i},20));  % substitute values of parameters
end
if ~isempty(f) && f(1)=='[' && f(end)==']'
  f = [f(2:end-1) ';'];
end
if any(f=='[')
  error('Arrays not allowed in definition of coordinate functions.')
end
f = strrep(f,x,'x');                        % rename first argument as x
F = cell(1,n);
for i = 1:n
  c = find(f==';',1);
  if isempty(c)
    F{i} = f;
  else
    F{i} = f(1:c-1);                        % string of ith coordinate function
    f = f(c+1:end);
  end
end

% call internal taylor
varargout = cell(1,max(1,nargout));
[varargout{:}] = intern_taylor(F,ainit(x0,k),y0);
% Draw a color-coded line by using surface with no facecolor
%
% SYNTAX
% ======
% h = cline(x, y, cdata [, linewidth])
%
% INPUT
% =====
% x                     vector with x-values
% y                     vector with y-values
% cdata                 vector with color-data
% linewidth             optional linewidth
%
% OUPUT
% =====
% h                 Handle to line
%
% Modified by PYC to use surface instead

function h = cline(x, y, c, lw)

z = zeros(size(x));

if nargin == 3
  lw = 1;
end

p = surface([x;x], [y;y], [z;z], [c;c], 'facecol', 'no', 'edgecol', 'interp', 'linew', lw);

if nargout == 1
  h = p;
end

function b = aget(a,varargin)

%AGET Select parts of audi data.
%   b = AGET(a,i) or b = AGET(a,i,j) read the function values and
%   derivatives with indices i or i,j and store it in the audi variable b.
%   It is size(b) = size(a), while asize(b)<=asize(a).
%
%   See also: aset

b = a;
if numel(a) > 1
  for i = 1:numel(a)
    b(i) = aget(a(i),varargin{:});
  end
else
  for i = 1:numel(a.c)
    b.c{i} = a.c{i}(varargin{:});
  end
end
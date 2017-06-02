function a = aset(a,b,varargin)

%ASET Overwrite data.
%   a = ASET(a,b,i) or a = ASET(a,b,i,j) replace the function values and
%   derivatives with indices i or i,j by those in the audi variable b.
%
%   See also: aget

for i = 1:numel(a.c)
  a.c{i}(varargin{:}) = b.c{i};
end
function y = adval(df,varargin)

%ADVAL Internal processing of derivative data. 

if numel(varargin)==1
  varargin{2} = 1;
end
ivar = varargin{end};
x = varargin{ivar};
arg = varargin(1:end-1);
s = dbstack;
[~,f] = s.name;
l = find(f=='.',1,'last');
if ~isempty(l)
  f = f(l+1:end);
end
if numel(x)>1
  y = x;
  for i = 1:numel(x)
    arg{ivar} = x(i);
    try
      y(i) = feval(f,arg{:});
    catch
      y(i) = feval(['audi.' f],arg{:});
    end
  end
  return
end
u = x.c{1};
if x.k>0
  E = eye(x.n);
  [~,J] = audi.mat(x.k-1,x.n);
  v = audi([],0,x.n,x.k-1);
  for m = 1:size(J,1)
    v.c{idx(v,J(m,:))} = x.c{idx(x,J(m,:))};
  end
  w = cell(1,x.n);
  for d = 1:x.n
    w{d} = audi([],0,x.n,x.k-1);
    for m = 1:size(J,1)
      w{d}.c{idx(w{d},J(m,:))} = x.c{idx(x,J(m,:)+E(d,:))};
    end
  end
else
  v = audi([],0,x.n,0);
  w = v;
end
if v.k==0
  v = v.c{1};
end
arg{ivar} = u;
try
  fval = feval(f,arg{:});
catch
  fval = feval(['audi.' f],arg{:});
end
arg{ivar} = v;
dfval = df(arg{:});
E = eye(w{1}.n);
y = audi([],0,w{1}.n,w{1}.k+1);
y.c{1} = fval; 
[~,J] = audi.mat(y.k,y.n);
for d = 1:length(w)
  w{d} = dfval.*w{d};
end
for m = 2:size(J,1)
  j = J(m,:);
  [~,d] = max(j);
  y.c{idx(y,j)} = w{d}.c{idx(w{d},j-E(d,:))};
end

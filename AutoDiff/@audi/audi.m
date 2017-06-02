classdef audi
  %AUDI MATLAB class for automatic differentiation.
  %   The AUDI class establishes automatic differentiation by means of
  %   operator overloading. It provides
  %
  %   * Ordinary and partial derivatives of arbitrary order
  %   * Differential operators, like Laplacian and curl
  %   * Taylor expansions of explicit expressions and solutions of ODEs
  %   * Curvature computation for curves and surfaces
  %   * Resolution of singularities of type "0/0" by l'Hospital's rule
  %
  %   Detailed information can be found in the <a href="matlab: ahelp(0)">Documentation</a> 
  %   and the <a href="matlab: ahelp(-1)">Examples</a>.
  %
  %   Version: 2.0
  %   Date:    May 09, 2016
  %   Author:  Ulrich Reif
  
  properties
    n       % number of variables
    k       % maxial order of differentiation
    c       % cell array of derivative data
  end
  
  methods
    % CONSTRUCTOR ---------------------------------------------------
    function r = audi(a,d,n,k)
      if nargin==0                             % empty audi
        r.c = [];
        return
      end
      if k<0
        error('Order must not be negative.')
      end
      z = zeros(1,n);                          % index vector of zeros
      r.n = n;
      r.k = k;
      r.c = cell(1,nchoosek(n+k,k));
      [r.c{:}] = deal(zeros(size(a)));
      r.c{1} = a;
      if d>0 && k>0
        z(d) = 1;
        r.c{idx(r,z)}(:) = 1;
      end
    end
    % ---------------------------------------------------------------
    
    % USER-DEFINED OVERLOADS ----------------------------------------
    % Functions you want to overload can be entered here.
    %
    % function y = fct(x)
    %   df = @(x) ...
    %   y = adval(df,x);
    % end
    % ---------------------------------------------------------------
    
    % FUNCTIONS -----------------------------------------------------
    function y = sin(x)
      df = @(x) cos(x);
      y = adval(df,x);
    end
    function y = asin(x)
      df = @(x) (1 - x.^2).^(-1/2);
      y = adval(df,x);
    end
    function y = cos(x)
      df = @(x) -sin(x);
      y  = adval(df,x);
    end
    function y = acos(x)
      y = pi/2 - asin(x);
    end
    function y = exp(x)
      df = @(x) exp(x);
      y  = adval(df,x);
    end
    function y = log(x)
      df = @(x) 1./x;
      y  = adval(df,x);
    end
    function y = tan(x)
      df = @(x) 1 + tan(x).^2;
      y  = adval(df,x);
    end
    function y = atan(x)
      df = @(x) 1./(1 + x.^2);
      y  = adval(df,x);
    end
    function y = atan2(x2,x1)
      y = angle(x1 + 1i*x2);
    end
    function y = sinh(x)
      df = @(x) cosh(x);
      y  = adval(df,x);
    end
    function y = asinh(x)
      df = @(x) (x.^2 + 1).^(-1/2);
      y  = adval(df,x);
    end
    function y = cosh(x)
      df = @(x) sinh(x);
      y  = adval(df,x);
    end
    function y = acosh(x)
      df = @(x) (x.^2 - 1).^(-1/2);
      y  = adval(df,x);
    end
    function y = tanh(x)
      df = @(x) 1 - tanh(x).^2;
      y  = adval(df,x);
    end
    function y = atanh(x)
      df = @(x) 1./(1 - x.^2);
      y  = adval(df,x);
    end
    function y = sqrt(x)
      y = x.^(1/2);
    end
    function y = abs(x)
      if isreal(x)
        y = -min(x,-x);
      else
        y = sqrt(x.*conj(x));
      end
    end
    function y = sign(x)
      df = @(x) 0*x;
      y  = adval(df,x);
    end
    function y = polyval(a,x)
      df = @(a,x) polyval(polyder(a),x);
      y = adval(df,a,x,2);
    end
    function y = cross(a,b)
      y = a;
      y(1) = a(2).*b(3) - a(3).*b(2);
      y(2) = a(3).*b(1) - a(1).*b(3);
      y(3) = a(1).*b(2) - a(2).*b(1);
    end
    function y = dot(a,b)
      y = a(:)'*b(:);
    end
    function y = sum(a,dim)
      if length(size(a))>2
        error('sum works for 2d arrays only.')
      end
      if nargin==1
        if size(a,1) == 1;
          dim = 2;
        else
          dim = 1;
        end
      end
      if dim==1
        y = ones(1,size(a,1))*a;
      else
        y = a*ones(size(a,2),1);
      end
    end
    function y = diff(a,n,dim)
      if length(size(a))>2
        error('diff works for 2d arrays only.')
      end
      if nargin<2
        n = 1;
      end
      if nargin<3
        if size(a,1) == 1;
          dim = 2;
        else
          dim = 1;
        end
      end
      if dim == 1
        y = diff(eye(size(a,1)),n,dim)*a;
      else
        y = a*diff(eye(size(a,2)),n,dim);
      end
    end
    function y = det(a,w)
      if nargin==1
        w = 1;
      end
      m = size(a,1);
      if size(a,2)~=m
        error('Matrix must be square.')
      end
      if m == 1
        y = a;
      elseif m==2;
        y = a(1,1)*a(2,2) - a(1,2)*a(2,1);
      else
        if m>4 && w
          warning(['det is likely to be slow for audi matrix of size ' num2str(m) '.'])
        end
        y = 0*a(1);
        for i = 1:size(a,1)
          I = 1:m;
          I(i) = [];
          y = y - (-1)^i*a(i,1)*det(a(I,2:end),0);
        end
      end
    end
    function y = inv(a)
      m = size(a,1);
      if size(a,2)~=m
        error('Matrix must be square.')
      end
      if m == 1
        y = 1./a;
      elseif m==2;
        D = det(a);
          y = [a(2,2) -a(1,2);-a(2,1) a(1,1)]./D;
      else
        if m>3
          warning(['inv is likely to be slow for audi matrix of size ' num2str(m) '.'])
        end
        y = a;
        D = det(a,0);
        for i = 1:m
          I = 1:m;
          I(i) = [];
          for j = 1:m
            J = 1:m;
            J(j) = [];
            y(j,i) = (-1)^(i+j)*det(a(I,J))./D;
          end
        end
      end
    end
    function y = trace(a)
      y = 0*a(1);
      for i = 1:size(a,1)
        y = y + a(i,i);
      end
    end
    function y = norm(a,p)
      if ~isvector(a)
        error('norm is not supported for audi matrices.')
      end
      if numel(a)==1
        y = abs(a);
        return
      end
      if nargin==1
        p = 2;
      end
      if isinf(p)
        if p>0
          y = max(abs(a(1)),norm(a(2:end),p));
        else
          y = min(abs(a(1)),norm(a(2:end),p));
        end
      else
        y = sum(a.^p).^(1/p);
      end
    end
    % ---------------------------------------------------------------
    
    % BINARY OPERATOS -----------------------------------------------
    function y = plus(a,b)
      if ~isa(a,'audi')
        y = plus(b,a);
        return
      end
      if ~isa(b,'audi')
        y = a;
        if numel(a)==1 && numel(b)>1
          for i = 1:numel(y.c)
            y.c{i} = a.c{i} + (i==1)*b;
          end
        elseif isequal(size(a),size(b))
          for i = 1:numel(a)
            yi = y(i);
            yi.c{1} = yi.c{1} + b(i);
            y(i) = yi;
          end
        elseif numel(b)==1
          y = a;
          for i = 1:numel(a)
            yi = y(i);
            yi.c{1} = yi.c{1} + b;
            y(i) = yi;
          end
        else
          error('Cannot add audi and double arrays of different size.')
        end
        return
      end
      if numel(a) == 1
        a = repmat(a,size(b));
      end
      if numel(b) == 1
        b = repmat(b,size(a));
      end
      if numel(a) > 1
        y = a;
        for i = 1:numel(a)
          y(i) = a(i) + b(i);
        end
        return
      end
      [a,b] = adapt(a,b);
      y = a;
      for i = 1:numel(a.c)
        y.c{i} = a.c{i} + b.c{i};
      end
    end
    function y = minus(a,b)
      y = a + (-1.*b);
    end
    function y = times(a,b)
      if ~isa(a,'audi') || (isa(b,'audi') && numel(b)>numel(a))
        y = times(b,a);
        return
      end
      if numel(a) > 1
        if numel(b) == 1
          b = repmat(b,size(a));
        end
        y = a;
        for i = 1:numel(a)
          y(i) = a(i).*b(i);
        end
        return
      end
      if ~isa(b,'audi')
        y = a;
        for i = 1:numel(a.c)
          y.c{i} = a.c{i}.*b;
        end
      else
        [a,b] = adapt(a,b);
        y = a;
        if a.k==0
          y.c{1} = (a.c{1}).*(b.c{1});
          return
        end
        [~,J] = audi.mat(a.k,a.n);
        for m = 1:size(J,1)
          j = J(m,:);
          s = 0*a.c{1};
          for mm = 1:size(J,1)
            jj = J(mm,:);
            if all(jj<=j)
              p = 1;
              for d = 1:a.n
                p = p*nchoosek(j(d),jj(d));
              end
              s = s + p*a.c{idx(a,j-jj)}.*b.c{idx(b,jj)};
            end
          end
          y.c{idx(y,j)} = s;
        end
      end
    end
    function y = rdivide(a,b)
      y = a.*(b.^(-1));
      if numel(a)*numel(b)==1 && isa(a,'audi') && isa(b,'audi') && adim(a)*adim(b)==1
        i = a==0 & b==0;
        if any(i)
          aa = aget(a,i);
          bb = aget(b,i);
          kk  = aord(a);
          for j = 1:kk
            aa.c{j} = aa.c{j+1}/max(1,j);
            bb.c{j} = bb.c{j+1}/max(1,j);
          end
          aa.c{kk+1} = nan(asize(aa));
          bb.c{kk+1} = nan(asize(bb));
          y = aset(y,aa./bb,i);
        end
      end
    end
    function y = ldivide(a,b)
      y = rdivide(b,a);
    end
    function y = mtimes(a,b)
      if (isnumeric(a)&&numel(a)==1) || (isnumeric(b)&&numel(b)==1)
        y = a.*b;
      elseif isnumeric(b)
        y = (mtimes(b',a'))';
      elseif size(a,2)==size(b,1)
        y(size(a,1),size(b,2)) = b(1);
        for i = 1:size(a,1)
          for j = 1:size(b,2)
            y(i,j) = 0*(a(1) + b(1));
            for m = 1:size(a,2)
              y(i,j) = y(i,j) + a(i,m).*b(m,j);
            end
          end
        end
      elseif isnumeric(a) && size(a,2)==size(b.c{1},1)
        y = b;
        for i = 1:numel(b.c)
          y.c{i} = a*b.c{i};
        end
      else
        error('Inner matrix dimensions must agree.')
      end
    end
    function y = mldivide(a,b)
      if size(a,1)==size(a,2)
        y = mtimes(inv(a),b);
      elseif size(a,1)>size(a,2)
        y = inv(a'*a)*(a'*b);
      else
        error('Solution of underdetermined systems is not supported.')
      end
    end
    function y = mrdivide(a,b)
      y = mldivide(b',a')';
    end
    function y = power(a,b)
      if isnumeric(b)
        if b==0
          y = 0*a + 1;
        elseif b==1
          y = a;
        else
          df = @(x,m) m*power(x,m-1);
          y = adval(df,a,b,1);
        end
      else
        y = exp(b.*log(a));
      end
    end
    function y = mpower(a,b)
      if isnumeric(b) && b>=0 && b==round(b)
        if b==0
          y = 0*a + eye(size(a));
        elseif b==1
          y = a;
        elseif b==2
          y = a*a;
        elseif b<8
          y = a*mpower(a,b-1);
        else
          bb = floor(b/2);
          y = mpower(a,b-2*bb)*mpower(mpower(a,bb),2);
        end
      else
        error('Exponent must be positive integer.')
      end
    end
    function y = min(a,b)
      % MIN Minimum of audi scalars.
      %    y = MIN(a,b) returns audi scalar with minimum values.
      %    If a and b are univariate, derivatives are taken from the right.
      %  
      %    Example: If x = ainit(0,2) then
      %                aeval(min(0,x),1) = 0
      %                aeval(min(0,-x^2),2) = -2
      %    See also: audi.max
      [a,b] = adapt(a,b);
      a0 = aeval(a,0);
      b0 = aeval(b,0);
      y = a;
      if isreal(a) && isreal(b)
        if adim(a)==1
          i = a0>b0;
          e = a0==b0;
          for d = 1:aord(a)
            i = i | (e & (aeval(a,d)>aeval(b,d)));
            e = e & (aeval(a,d)==aeval(b,d));
          end
        else
          i = a0>b0;
        end
      else
        i = abs(a0)>abs(b0)
      end
      y = aset(y,aget(b,i),i);
    end
    function y = max(a,b)
      % MAX Maximum of audi scalars.
      %    y = MAX(a,b) returns audi scalar with maximal values.
      %    If a and b are univariate, derivatives are taken from the right.
      %  
      %    Example: If x = ainit(0,2) then
      %                aeval(max(0,x^2),2) = 2
      %                aeval(max(0,-x),1) = 0
      %    See also: audi.min
      y = -min(-a,-b);
    end
    function y = real(a)
      y = a;
      for i = 1:numel(a)
        b = a(i);
        for j = 1:numel(b.c)
          b.c{j} = real(b.c{j});
        end
        y(i) = b;
      end
    end
    function y = imag(a)
      y = a;
      for i = 1:numel(a)
        b = a(i);
        for j = 1:numel(b.c)
          b.c{j} = imag(b.c{j});
        end
        y(i) = b;
      end
    end
    function y = conj(a)
      y = a;
      for i = 1:numel(a)
        b = a(i);
        for j = 1:numel(b.c)
          b.c{j} = conj(b.c{j});
        end
        y(i) = b;
      end
    end
    function y = angle(a)
      y = imag(log(a));
    end
    function y = isreal(a)
      y = true;
      for i = 1:numel(a)
        b = a(i);
        for j = 1:numel(b.c)
          y = y && isreal(b.c{j});
        end
      end
    end
    function y = isfinite(a)
      y = true;
      for i = 1:numel(a)
        b = a(i);
        for j = 1:numel(b.c)
          y = y && isfinite(b.c{j});
        end
      end
    end
    % ---------------------------------------------------------------
    
    % UNARY OPERATORS -----------------------------------------------
    function y = uminus(a)
      y = (-1).*a;
    end
    function y = uplus(a)
      y = a;
    end
    % ---------------------------------------------------------------
    
    % RELATIONAL OPERATORS ------------------------------------------
    function y = le(a,b)
      y = le(aeval(a,0),aeval(b,0));
    end
    function y = ge(a,b)
      y = ge(aeval(a,0),aeval(b,0));
    end
    function y = lt(a,b)
      y = lt(aeval(a,0),aeval(b,0));
    end
    function y = gt(a,b)
      y = gt(aeval(a,0),aeval(b,0));
    end
    function y = eq(a,b)
      y = eq(aeval(a,0),aeval(b,0));
    end
    function y = ne(a,b)
      y = ne(aeval(a,0),aeval(b,0));
    end
    function y = and(a,b)
      y = and(aeval(a,0),aeval(b,0));
    end
    function y = or(a,b)
      y = or(aeval(a,0),aeval(b,0));
    end
    function y = xor(a,b)
      y = xor(aeval(a,0),aeval(b,0));
    end
    function y = not(a)
      y = not(aeval(a,0));
    end
    % ---------------------------------------------------------------
    
    % AUXILIARY ROUTINES --------------------------------------------
    function y = subsref(a,s)
      if ~strcmp(s(1).type,'{}')
        y = builtin('subsref',a,s);
      else
        y = aeval(a,[s(1).subs{:}]);
        if numel(s)>1
          y = builtin('subsref',y,s(2:end));
        end
      end
    end
    function varargout = intern_taylor(D,x0,y0)
      if nargin==1
        a = D;
        if prod(asize(a)) > 1
          error('Evaluation at single points only.')
        elseif size(a,2) > 1
          error('Evaluation for column vectors only.')
        elseif adim(a) > 1
          error('Evaluation for univariate functions only.')
        else
          y = zeros(size(a,1),aord(a)+1);
          for i = 1:size(a,1)
            b = a(i);
            for j = 0:b.k
              y(i,end-j) = b.c{j+1}/factorial(j);
            end
          end
          [i,j] = find(isnan(y));
          if ~isempty(i)
            y = y(:,max(j)+1:end);
          end
        end
        varargout{1} = y;
      else
        if ~iscell(D)
          D = {D};
        end
        N = numel(D);
        for i = 1:N
          s = D{i};
          for j = 1:N
            yold = ['y' num2str(j)];
            ynew = ['audi.odeaux(x,' num2str(j) ',D,y0)'];
            s = strrep(s,yold,ynew);
          end
%           if N==1
%             s = strrep(s,'y','audi.odeaux(x,1,D,y0)');
%           end
          D{i} = str2func(['@(x,b,D,y0) ' s]);
        end
        varargout = cell(1,min(N,nargout));
        for i = 1:min(N,nargout)
          varargout{i} = intern_taylor(audi.odeaux(x0,i,D,y0));
        end
      end
    end
    function i = idx(a,s)
      [C,~] = audi.mat(a.k,a.n);
      s = num2cell(s+1);
      i = C(s{:});
    end
    function s = sub(a,i)
      [~,J] = audi.mat(a.k,a.n);
      s = J(i,:);
    end
    function [a,b] = adapt(a,b)
      %ADAPT Generate compatible audi scalars.
      %   [aa,bb] = ADAPT(a,b) returns two compatible audis.
      %   Input arguments of class double become audi constants.
      %   If a and b are audi scalars, aa and bb have equal order,
      %   aord(aa) = aord(bb) = min(aord(a),aord(b)).
      if isnumeric(a)
        a = 0*b + a;
      elseif isnumeric(b)
        b = 0*a + b;
      elseif a.n ~= b.n
        error('Audi variables must have equal dimension.')
      elseif a.k < b.k
        aa = a;
        for i = 1:numel(a.c)
          aa.c{i} = b.c{idx(b,sub(aa,i))};
        end
        b = aa;
      elseif a.k > b.k
        [b,a] = adapt(b,a);
      end
    end
    function [a1,a2,p] = prepfct(~,d,k)
      if nargin==2
        k = sum(d);
      else
        d = double(d);
        d(d==false) = nan;
      end
      i = isfinite(d);
      j = isnan(d);
      n = numel(d);
      p = cell(1,n);
      [p{i}] = deal(0);
      [p{i}] = ainit(p{i},k);
      [p{j}] = deal(0);
      d = d(i);
      a1 = '';
      a2 = '';
      for i = 1:n
        a1 = [a1 ',x' num2str(i)];
        a2 = [a2 ',x' num2str(i) '+p{' num2str(i) '}'];
      end
      a1 = a1(2:end);
      a2 = a2(2:end);
    end
  end
  % --------------------------------------------------------------------------------------
  methods(Static)
    function [C,J] = mat(k,n)
      persistent CC JJ
      if k==0
        C = 1;
        J = zeros(1,n);
        return
      end
      if all(size(CC)>=[k,n]) && ~isempty(CC{k,n})
        C = CC{k,n};
        J = JJ{k,n};
      else
        if n==1
          C = zeros(1,k+1);
        else
          C = zeros((k+1)+zeros(1,n));
        end
        K = C;
        K(:) = rem(0:numel(C)-1,k+1);
        J = zeros(numel(C),n);
        J(:,1) = K(:);
        for m=2:n
          K = permute(K,[2:n 1]);
          J(:,m) = K(:);
        end
        f = find(sum(J,2)<=k);
        J = J(f,:);
        for i = 1:length(f)
          c = num2cell(J(i,:) + 1);
          C(c{:}) = i;
        end
        CC{k,n} = C;
        JJ{k,n} = J;
      end
    end 
    function f = odeaux(t,b,D,y0)
      if isnumeric(t)
        f = y0(b);
      else
        f  = adval(D{b},t,b,D,y0,1);
      end
    end
  end
end


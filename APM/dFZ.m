function df = dFZ(z)
  % ***********************************************************************
  % 
  %      Computes the derivative of the function "f"
  %      at the point "z".
  % 
  % ***********************************************************************


sqrteps = (eps)^(1/3); % Used to make h.  
h = sqrteps*z; % From numerical recipes, make h = h(xr)
df = (FZ(z+h)-FZ(z-h))./(2*h); % First Derivative
% d2f = (FZ(z+h)+FZ(z-h) - 2*FZ(z))./(h.^2); % Second Derivative
end

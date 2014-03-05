% Filter coefficients in the spectral domain of a Gaussian filter.
%
% Input:
% radius: half-with radius parameter in km.
% maxDegree: maximum degree to compute.
%
% Output:
% wn: (n + 1) x 1 vector with n = maxDegree.
% The vector element wn(n + 1) contains the filter coefficient of degree n.

function wn = filterCoefficientsGaussian(radius, maxDegree)

R = 6378.137; % Earth radius [km]
b = log(2.)/(1-cos(radius/R));

wn = zeros(maxDegree+1,1);
wn(1) = 1;
wn(2) = (1+exp(-2*b))/(1-exp(-2*b)) - 1/b;
for n=2:maxDegree
  wn(n+1) = -(2*n-1)/b * wn(n) + wn(n-1);
  % Stop in the case where the recursion gets numerically unstable
  if wn(n+1)>wn(n) | wn(n+1)<0  
      wn(n+1) = 0;
      break
  end
  
end

% Spherical harmonic synthesis 
% according to 
%
% f(lambda,theta) = ...
%     sum_{n=0...maxDegree} sum_{m=0...n} ...
%    P_nm(cos(theta))*( C_nm*cos(m*lambda) + S_nm*sin(m*lambda) ),
%
% where lambda and theta are spherical longitude and colatitude in radians.
% No additional factors are applied to the input spherical
% harmonic coefficients. To calculate geodetic magnitudes like geoid etc.,
% those factors have to be applied in advance.
%
% Synthax:
%         Grid = SH_Synthesis(lambda, theta, Coeff)
% Input:
% lambda: p x 1 vector containing the longitudes (in radian).
% theta: q x 1 vector containing the co-latitudes (in radian).
% Coeff: matrix containing the spherical harmonic coefficients
% Output:
% Grid: Grid of synthetized values


function [Grid, Grid_Err] = SH_Synthesis_ERR(lambda, theta, Coeff, Coeffer)
maxDegree = size(Coeff,1)-1;
Cnm = Coeff(:,maxDegree+1:end);
Snm = fliplr(Coeff(:,1:maxDegree+1));
Cnmerr = Coeffer(:,maxDegree+1:end);
Snmerr = fliplr(Coeffer(:,1:maxDegree+1));

% cos,sin matrix for all lambdas and all orders
cosm = cos([0:1:maxDegree]'*lambda');
sinm = sin([0:1:maxDegree]'*lambda');

Grid = zeros(length(theta), length(lambda));
Grid_Err = zeros(length(theta), length(lambda));
for k=1:length(theta) % loop over all thetas
  Pnm = legendreFunctions(theta(k), maxDegree);  % all Legendre Functions for one theta
  % compute result for all lambdas (one row) in one step (as matrix multiplications)
   Grid(k,:) =  sum( (Cnm.*Pnm) * cosm + (Snm.*Pnm) * sinm );
   Grid_Err(k,:) = sum(((Cnmerr.*Pnm).*(Cnmerr.*Pnm))* (cosm .* cosm) + ...
       (Snmerr.*Pnm).*(Snmerr.*Pnm) * (sinm.*sinm));
end


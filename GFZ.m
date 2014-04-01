clear;
Rearth = 6378137;
max_degree = 10;
Cut_degree = 10;

type_gfz = 'gfz_05';
pathGFZ='/home/sghelichkhani/Workplace/DATA/GFZ';
%pathGFZ = '/home/sghelichkhani/Workplace/GRACE_DATA/GFZ';

%% Reading Coefficients in
fprintf(['Readin Data ', type_gfz, '\n'])
[CoefGFZ, mjdmidGFZ, yearmidGFZ, regulflagGFZ] = ...
        read_grace_solution_series_20140128(type_gfz, ...
        pathGFZ ,max_degree);
%% Filter with the cut off degree specified at top
gauss_radius = 20000/Cut_degree;
GAUSS_FILTER = ...
filterCoefficientsGaussian(gauss_radius,max_degree);
%% Calculation of Lambda Theta
lambda = pi/180 * (0:1:360)';    % Longitude
theta  = pi/180 * (0:1:180)';       % Co-latitude
% phi2 = pi/2 - theta;                 % Latitude

% The coastlines are loaded and corrected for our geometry
load coast3;
for i=1:length(phi)
    phi(i) = -phi(i)+90;
    if(lam(i) < 0)
         lam(i)=lam(i)+360;
    end
end

%% Conversion to the field on the Sphere
%   GFZ Conversion
for i=1:length(mjdmidGFZ)
    fprintf(['synthesis  ',num2str(i),' of ',num2str(length(mjdmidGFZ)),' File:', type_gfz,'\n'])
    geoGFZ(:,:,i) = ...
    SH_Synthesis(lambda, theta, diag(GAUSS_FILTER)*CoefGFZ(:,:,i));
end
geoGFZ = Rearth*geoGFZ*1e6;
GFZmean10 = mean(geoGFZ,3);


%% PLOTTING
 for i=1:length(geoGFZ(1,1,:))
   figure(1);
   imagesc(geoGFZ(:,:,i)-GFZmean10)
   hold on;
   plot(lam, phi,'k.', 'MarkerSize', 2)
   hold off;
   colorbar;
   str2 = sprintf('TIME =  %f', yearmidGFZ(i));
   title(str2)
   pause;
   clf;
 end


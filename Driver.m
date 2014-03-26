clear;
Rearth = 6378137;
max_degree = 60;
Cut_degree = 15;

type_gfz = 'gfz_05';
%type_gfz = 'jpl_05';
%type_jpl = 'jpl_05';

%pathGFZ = '/home/benutzer/Workspace/DATA/JPL';
%/home/sghelichkhani/Workplace/DATA/GFZ
pathGFZ = '/home/sghelichkhani/Workspace/DATA/TEST'
%pathGFZ = '/home/sghelichkhani/Workspace/DATA/JPL';
%pathCSR = '/home/benutzer/Workplace/DATA/CSR';
%pathJPL = '/home/benutzer/Workplace/DATA/JPL';

%% Reading Coefficients in
fprintf(['Readin Data ', type_gfz, '\n'])
[CoefGFZ, ErrCoef, mjdmidGFZ, yearmidGFZ, regulflagGFZ] = ...
        read_grace_solution_series_20140128_ERR(type_gfz, ...
        pathGFZ ,max_degree);
fprintf(['Readin Data ', type_gfz, '\n'])
[CoefGFZ, mjdmidGFZ, yearmidGFZ, regulflagGFZ] = ...
        read_grace_solution_series_20140128(type_gfz, ...
        pathGFZ ,max_degree);
CoefGFZ(3,61,:)=0;
ErrCoef(3,61,:)=0;
%% Filter with the cut off degree specified at top
gauss_radius = 20000/Cut_degree;
GAUSS_FILTER = ...
filterCoefficientsGaussian(gauss_radius,max_degree);
%% Calculation of Lambda Theta
lambda = pi/180 * [0:1:360]';    % Longitude
theta  = pi/180 * [0:1:180]';       % Co-latitude
phi = pi/2 - theta;                 % Latitude

load Coastlines;
% coast_phi = -coast_phi+90;
% coast_lam = coast_lam+180;
for i=1:length(coast_phi)
    if(coast_lam(i) <= 0)
         coast_lam(i)=coast_lam(i)+360;
    end
end


%% Conversion to the field on the Sphere
%   GFZ Conversion
for i=1:length(mjdmidGFZ)
    fprintf(['synthesis  ',num2str(i),' of ',num2str(length(mjdmidGFZ)),' File:', type_gfz,'\n'])
    [geoGFZ(:,:,i), Error_GFZ(:,:,i)] = ...
    SH_Synthesis_ERR(lambda, theta, diag(GAUSS_FILTER)*CoefGFZ(:,:,i),...
    diag(GAUSS_FILTER)*ErrCoef(:,:,i));
end

geoGFZ = Rearth*geoGFZ*1e6;
Error_GFZ=Rearth*Error_GFZ*1e6;
GFZmean10 = mean(geoGFZ,3);


%% Mean changes
lati = 1;
latf = 181;
loni = 1;
lonf = 361;

trend = zeros(latf-lati+1,lonf-loni+1);

for lat0=lati:latf
    fprintf(['LAT: ',num2str(lat0),'\n'])
    for lon0=loni:lonf
        for i=1:length(geoGFZ(1,1,:))
            gfz(i) = (geoGFZ(lat0,lon0,i)-GFZmean10(lat0,lon0));
        end
        fitgfz = polyfit(yearmidGFZ,gfz',1);       
        trend(lat0,lon0)= fitgfz(1);
    end
end

% The coastlines are loaded and corrected for our geometry
figure;imagesc(lambda*180/pi,phi*180/pi,trend);
axis xy;
hold on;
plot(coast_lam,coast_phi,'k.', 'MarkerSize', 2);
colorbar;
caxis([-3e3 +3e+3])
title(type_gfz)
hold off;

%% PLOT MAP
for i=1:length(yearmidGFZ)
    figure(1);
    imagesc(lambda*180/pi,phi*180/pi, ...
         geoGFZ(:,:,i));
%         Error_GFZ(:,:,i));
    axis xy;
    hold on;
    plot(coast_lam, coast_phi, 'k.','MarkerSize', 2);
    hold off;
    colorbar;
    %caxis([-1e-3, +1e+3])
    title(num2str(yearmidGFZ(i)));
    pause;
    clf;
end






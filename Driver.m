clear;
Rearth = 6378137;
max_degree = 60;
Cut_degree = 15;
type = 'gfz_05';
path = '/home/siavash/Workplace/DATA/GRACE/GFZ';
%% Reading Coefficients in
fprintf(['Readin Data ', type, '\n'])
[sphcoef, errcoef, mjdmid, yearmid, outflags] = ...
        read_grace_solution_series_20140128_ERR(type, ...
        path ,max_degree);
sphcoef(3,61,:)=0;
errcoef(3,61,:)=0;

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

gridgeo = zeros(length(theta), length(lambda),length(mjdmid));
griderr = zeros(length(theta), length(lambda),length(mjdmid));
mulflag=1;
for i=1:length(mjdmid)
    fprintf(['synthesis  ',num2str(i),' of ',num2str(length(mjdmid)),' File:', type,'\n'])
    [gridgeo(:,:,i), griderr(:,:,i)] = ...
    SH_Synthesis_ERR(lambda, theta, diag(GAUSS_FILTER)*sphcoef(:,:,i),...
    (diag(GAUSS_FILTER).*diag(GAUSS_FILTER))*errcoef(:,:,i));
end
if mulflag==1
    gridgeo = Rearth*gridgeo*1e6;
    griderr = Rearth*griderr*1e6;
%    meangeo = mean(gridgeo,3);
%     for i=1:length(mjdmid)
%         gridgeo(:,:,i) = gridgeo(:,:,i) -meangeo;
%     end
    mulflag=0;
end
%% Mean changes
lati = 1;
latf = 181;
loni = 1;
lonf = 361;

trend = zeros(latf,lonf:lonf,2);
trenerr = zeros(latf,lonf:lonf,2);
for lat0=lati:latf
    fprintf(['LAT: ',num2str(lat0),'\n'])
    for lon0=loni:lonf
        nodedat = zeros(length(mjdmid),1);
        nodeerr = zeros(length(mjdmid),length(mjdmid));
        Amatrix = ones(length(mjdmid),2);
        for i=1:length(yearmid)
            nodedat(i,1) = gridgeo(lat0,lon0,i);
            nodeerr(i,i) = 1/griderr(lat0,lon0,i);
            Amatrix(i,2) = yearmid(i);
        end
        varian = inv(Amatrix'*nodeerr*Amatrix);
        trenerr(lat0, lon0, 1) = varian(1,1);
        trenerr(lat0, lon0, 2) = varian(2,2);
        trend(lat0,lon0,:)=varian*Amatrix'*nodeerr*nodedat;
%        fitgfz = polyfit(yearmid,nodedat,1);
    end
end
%% Plotting the Trend
% The coastlines are loaded and corrected for our geometry
figure;imagesc(lambda*180/pi,phi*180/pi,trend(:,:,2));
axis xy;
hold on;
plot(coast_lam,coast_phi,'k.', 'MarkerSize', 2);
colorbar;
%caxis([-1e3, +1e3])
title(type)
hold off;

%% PLOT: Snapshot of the difference
for i=1:length(yearmid)
    figure(2);
    imagesc(lambda*180/pi,phi*180/pi,griderr(:,:,i));
    axis xy;
    hold on;
    plot(coast_lam, coast_phi, 'k.','MarkerSize', 2);
    
    colorbar;
%    caxis([-1e+1, +1e+1])
    title(num2str(yearmid(i)));
    pause;
    clf;
end
%%






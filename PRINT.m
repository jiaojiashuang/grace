clear;
Rearth = 6378137;
max_degree = 100;
Cut_degree = 5;

type_gfz = 'gfz_05';
type_csr = 'csr_05';
type_jpl = 'jpl_05';
pathGFZ = '/home/benutzer/Workplace/DATA/GFZ';
pathCSR = '/home/sghelichkhani/Workplace/DATA/CSR';
pathJPL = '/home/sghelichkhani/Workplace/DATA/JPL';

%% Reading Coefficients in
fprintf(['Readin Data ', type_gfz, '\n'])
[CoefGFZ, mjdmidGFZ, yearmidGFZ, regulflagGFZ] = ...
        read_grace_solution_series_20140128(type_gfz, ...
        pathGFZ ,max_degree);
fprintf(['Readin Data ', type_csr, '\n'])
[CoefCSR, mjdmidCSR, yearmidCSR, regulflagCSR] = ...
        read_grace_solution_series_20140128(type_csr, ...
        pathCSR ,max_degree);
fprintf(['Readin Data ', type_jpl, '\n'])
[CoefJPL, mjdmidJPL, yearmidJPL, regulflagJPL] = ...
        read_grace_solution_series_20140128(type_jpl, ...
        pathJPL ,max_degree);
%% Filter with the cut off degree specified at top
gauss_radius = 20000/Cut_degree;
GAUSS_FILTER = ...
filterCoefficientsGaussian(gauss_radius,max_degree);
%% Calculation of Lambda Theta
lambda = pi/180 * [0:1:360]';    % Longitude
theta  = pi/180 * [0:1:180]';       % Co-latitude
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
%   CSR Conversion
for i=1:length(mjdmidCSR)
    fprintf(['synthesis  ',num2str(i),' of ',num2str(length(mjdmidCSR)),' File:', type_csr,'\n'])
    geoCSR(:,:,i) = ...
    SH_Synthesis(lambda, theta, diag(GAUSS_FILTER)*CoefCSR(:,:,i));
end
geoCSR = Rearth*geoCSR*1e6;
CSRmean10 = mean(geoCSR,3);
%   JPL Conversion
for i=1:length(mjdmidJPL)
    fprintf(['synthesis  ',num2str(i),' of ',num2str(length(mjdmidJPL)),' File:', type_jpl,'\n'])
    geoJPL(:,:,i) = ...
    SH_Synthesis(lambda, theta, diag(GAUSS_FILTER)*CoefJPL(:,:,i));
end
geoJPL = Rearth*geoJPL*1e6;
JPLmean10 = mean(geoJPL,3);

%% Mean changes
lati = 1;
latf = 181;
loni = 1;
lonf = 361;
for lat0=lati:10:latf
    for lon0=loni:10:lonf
        fprintf(['LAT: ',num2str(lat0),' LON: ',num2str(lon0),'\n'])
        for i=1:length(geoGFZ(1,1,:))
            gfz(i) = geoGFZ(lat0,lon0,i)-GFZmean10(lat0,lon0);
        end
        for i=1:length(geoCSR(1,1,:))
            csr(i) = geoCSR(lat0,lon0,i)-CSRmean10(lat0,lon0);
        end
        for i=1:length(geoJPL(1,1,:))
            jpl(i) = geoJPL(lat0,lon0,i)-JPLmean10(lat0,lon0);
        end
        fitgfz = polyfit(yearmidGFZ,gfz',1);
        fitcsr = polyfit(yearmidCSR',csr,1);
        fitjpl = polyfit(yearmidJPL',jpl,1);
        
        hFig=figure(1);
        clf;
        plot(yearmidGFZ', gfz, 'Color', 'b', 'LineWidth',1);
        set(hFig, 'Position', [0 0 2500 1000], 'Visible', 'off')
        hold on;
        plot(yearmidCSR', csr, 'Color', 'r', 'LineWidth',1);
        hold on;
        plot(yearmidJPL', jpl, 'Color', 'k', 'LineWidth',1);
        hold on;
        plot(yearmidGFZ', fitgfz(1)*yearmidGFZ +fitgfz(2), ...
        'b', 'LineStyle', '-.', 'LineWidth', 2.0);
        hold on;
        plot(yearmidCSR', fitcsr(1)*yearmidCSR +fitcsr(2), ...
        'r', 'LineStyle', '-.', 'LineWidth', 2.0 );
        hold on;
        plot(yearmidJPL', fitjpl(1)*yearmidJPL +fitjpl(2), ...
        'k', 'LineStyle', '-.', 'LineWidth', 2.0 );
        hold on;
        Name=sprintf('TRENDS/Lat%iLon%i.png', lat0, lon0);
        title_fig=sprintf('Lat%iLon%i', lat0, lon0);
        title(title_fig);
        
        legend('GFZ','CSR','JPL', 'Location', 'NorthWest')
        hold off;
        saveas(hFig,Name,'png')
        
    end
end


%% PLOTTING
 for i=1:length(geoGFZ(1,1,:))
   figure(i,'Visible','on');
   imagesc(geoGFZ(:,:,i)-GFZmean10)
   hold on;
   plot(lam, phi,'k.', 'MarkerSize', 2)
   colorbar;
   str2 = sprintf('TIME =  %f', yearmidGFZ(i));
   title(str2)
   pause;
   clf;
 end
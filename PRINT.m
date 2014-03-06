clear;
type_data = 'csr_05';
path = '/home/sghelichkhani/Workplace/GRACE_DATA/CSR';
max_degree = 10;
gauss_radius = 1;

[GRACE_COEF, mjdmid, yearmid, regulflag] = ...
read_grace_solution_series_20140128(type_data, ...
path ,max_degree);

%%
GAUSS_FILTER = ...
filterCoefficientsGaussian(gauss_radius,max_degree);
%%
lambda = pi/180 * [0:1:360]';    % Longitude
theta  = pi/180 * [0:1:180]';       % Co-latitude
phi2 = pi/2 - theta;                 % Latitude
%%
for i=1:length(mjdmid)
    fprintf(['synthesis  ',num2str(i),' of ',num2str(length(mjdmid)),'\n'])
    GRACEgeo(:,:,i) = ...
        SH_Synthesis(lambda, theta, diag(GAUSS_FILTER)*GRACE_COEF(:,:,i));
end
%% Converting the parameter to meters
Rearth = 6378137;
GRACEgeo = Rearth*GRACEgeo*1e6;
GraceMean = mean(GRACEgeo,3);
for i=1:length(mjdmid)
    fprintf(['Deviation ',num2str(i),' of ',num2str(length(mjdmid)),'\n'])
    GRACEgeo(:,:,i) = GRACEgeo(:,:,i) - GraceMean;
end
%% PLOTING SECTION

load coast3;
% Here comes lam phi;
for i=1:length(phi)
    phi(i) = -phi(i)+90;
    if(lam(i) < 0)
         lam(i)=lam(i)+360;
    end
end
%%
figure(2);
for i=1:length(GRACEgeo(1,1,:))
   imagesc(GRACEgeo(:,:,i))
   hold on;
   plot(lam, phi,'k.', 'MarkerSize', 2)
   colorbar;
   str2 = sprintf('TIME =  %f', yearmid(i));
   title(str2)
   pause;
   clf;
end

%%
lat0 = 140;
lon0 = 350;

for i=1:length(GRACEgeo(1,1,:))
    yy(i) = GRACEgeo(lat0,lon0,i);
end
% 
% xx = linspace(1,length(yy), length(yy));
% p = polyfit(xx,yy,1);
% figure;
 plot(yearmid, yy, 'Color', 'g');
 %xlim([min(yearmid), max(yearmid)])
% title(type_data)
% hold on


% plt1 = plot(xx, p(1)*yy+p(2),'Color','r');
% plt1 = title('LAT, LONG =  %i , %i \n %f y + x %f', ...
%    lat0, lon0, p(1), p(2));
% plt1 = title(num2str(p(1), num2str(p(1))));
% plt1 = legend('GEOID VARIATIONS', 'POLYNOMIAL FIT');












pcolor(log10(abs(xx)))

%% Load the Coastlines

load coast3
for i=1:length(phi)
    if(lam(i) < 0)
        lam(i)=lam(i)+360;
    end
end

hold on
plot(lam,90-phi,'k');

%% 
 
 imagesc(GRACEgeo_mean2);
 

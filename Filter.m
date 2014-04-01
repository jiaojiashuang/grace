Cut_degree=5;
max_degree=100;
gauss_radius = 2000/Cut_degree;
GAUSS_FILTER = ...
filterCoefficientsGaussian(gauss_radius,max_degree);

figure;plot([1:length(GAUSS_FILTER)], GAUSS_FILTER);
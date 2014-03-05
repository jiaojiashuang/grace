function year = mjd2year(mjd)

% Calculation of decimal year AD from Modified Julian Day
% Syntax: year = mjd2year(mjd)
% mjd may be a scalar, vector, or array
% year will be ot the same dimension.

year = 2000*ones(size(mjd)) + ( mjd - 51544*ones(size(mjd)) ) / 365.25;

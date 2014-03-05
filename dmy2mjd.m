function mjd = dmy2mjd( dat )

% Conversion from (Day, Month, Year) format to Modified Julian Day.
% Syntax: 
%		mjd = dmy2mjd( dmy_matrix )
%
% dmy_matrix is a 3-column-matrix.
% where the columns are day, month, and year.
% mjd is a colulmn vector

mjd = zeros(size(dat,1),1);

for k = 1:size(dat,1)
  Y = dat(k,3); M = dat(k,2); D = dat(k,1);
  A=367*Y;
  B = fix(7*(Y+fix((M+9)/12))/4);
  C = fix( 3*(fix((Y+(M-9)/7)/100) +1 ) /4 );
  E = fix(275*M/9) + D + 1721028.5;
  JD = A-B-C+E;
  mjd(k) = JD - 2400000.5;
end


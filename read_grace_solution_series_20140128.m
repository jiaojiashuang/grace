function [MON, mjdmid, yearmid, regulflag] = ...
    read_grace_solution_series_20140128(releaseflag, releasepath, maxdeg)

% Read series of GRACE solutions (or AOD background fields)
% and write them to triangular (SC) format, with the 3rd dimension being
% the time.
% Syntax: [MON, mjdmid, yearmid, regulflag] = ...
%    read_grace_solution_series_20130729(releaseflag, releasepath, maxdeg)
% 
% releaseflag is a switch between releases
% 'gfz_05'
% 'csr_05'
% 'jpl_05'
% 'csr_05_gac'
% 'csr_05_gad'
% 'itg2010'
%
% releasepath is the path where the data can be found
% maxdeg is the maximum degree up to which the coefficients shall be
% stored.
%
% MON is the 3D matrix of the coefficients
% mjdmid are the Modified Julian Days of the midpoints of the solution
% intervals
% yearmid ist the same in decimal years
% regulaflag is a vector of logicals indicating whether the solutions are
% regularized or not. This only applies to the GFZ RL05 solutions.


%% Create list of available solutions

switch releaseflag
    case 'gfz_05'
        list=dir([releasepath,filesep,'GSM-2_*_EIGEN_*-_005a']);
    case 'gfz_05_60_alternativ'
        list=dir([releasepath,filesep,'GSM-2_*_EIGEN_*-60_5001']);
    case 'gfz_05_90_alternativ'
        list=dir([releasepath,filesep,'GSM-2_*_EIGEN_*-_5001']);
    case 'csr_05'
        list=dir([releasepath,filesep,'GSM-2_*_UTCSR_0060_0005']);
    case 'csr_05_96'
        list=dir([releasepath,filesep,'GSM-2_*_UTCSR_0096_0005']);
    case 'jpl_05'
        list=dir([releasepath,filesep,'GSM-2_*_JPLEM_0000_0005']);
    case 'itg2010'
        list=dir([releasepath,filesep,'ITG-Grace*.gfc']);
    case 'csr_05_gac'
        list=dir([releasepath,filesep,'GAC-2_*_UTCSR_0000_0005']);
    case 'csr_05_gad'
        list=dir([releasepath,filesep,'GAD-2_*_UTCSR_0000_0005']);
    otherwise
        error('Release flag unknown')
end
if size(list,1) == 0
    error('No file found')
end

%% Initialize output variables
perzahl = size(list,1);
mjdmid = zeros(perzahl, 1);
%regulflag = zeros(perzahl, 1);
MON = zeros( maxdeg+1, 2*maxdeg+1, perzahl);

%% Determine whether GFZ solutions are regularized and write this to
% regulflag

regulflag = false(perzahl,1);
if strcmp(releaseflag, 'gfz_05')
    for per = 1:perzahl
        regulflag(per) = strcmp( list(per).name(35:36), 'K2' );
    end
end

%% Read solutions

if ~strcmp(releaseflag,'itg2010') % non-ITG solutions, in GRACE SDS format
    for per = 1:perzahl
        monthID = list(per).name(7:26);
        mjdmid(per) = monthID2mjdmid_20111222(monthID);
        filename = fullfile(releasepath,list(per).name);
        fid = fopen(filename);
        s=fgets(fid);
        while s>0
            if strcmp(s(1:6),'GRCOF2')
                zwi = str2num(s(8:54));
                if zwi(1)<=maxdeg
                    MON(zwi(1)+1, maxdeg+1-zwi(2), per) = zwi(4);
                    MON(zwi(1)+1, maxdeg+1+zwi(2), per) = zwi(3);
                end
            end
            s=fgets(fid);
        end
        fclose(fid);
    end
else    % ITG solutions, in GFC format
    for per = 1:perzahl
        per
        mjdmid(per) = yearmon2mjdmid(str2num(list(per).name(15:18)), str2num(list(per).name(20:21)));
        filename = fullfile(releasepath,list(per).name);
        fid = fopen(filename);
        s=fgets(fid);
        while(strncmp(s, 'end_of_head', 11) == 0),
            s=fgets(fid);
        end
        while s>0
            if strcmp(s(1:3),'gfc')
                zwi = str2num(s(5:52));
                if zwi(1)<=maxdeg
                    MON(zwi(1)+1, maxdeg+1-zwi(2), per) = zwi(4);
                    MON(zwi(1)+1, maxdeg+1+zwi(2), per) = zwi(3);
                end
            end
            s=fgets(fid);
        end
        fclose(fid);
    end
end

yearmid = mjd2year(mjdmid);


    
    

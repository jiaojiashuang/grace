function [mjdmid,mjd_anf,mjd_end] = monthID2mjdmid_20111222(monthID)

    jahr_anf = monthID(1:4);
    tag_anf = monthID(5:7);
    jahr_end = monthID(9:12);
    tag_end = monthID(13:15);
    
    mjd_anf = dmy2mjd([1,1, str2num(jahr_anf)]) + str2num(tag_anf)-1;
    mjd_end = dmy2mjd([1,1, str2num(jahr_end)]) + str2num(tag_end)-1;
    mjdmid = mean([mjd_anf,mjd_end])+0.5; 

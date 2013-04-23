function optipotential_save( optipot )
%OPTIPOTENTIAL_SAVE Summary of this function goes here
%   Detailed explanation goes here

global optdatadir;

save([optdatadir,filesep,'optipot_',optipot.problemname,'.mat'],'optipot');


end


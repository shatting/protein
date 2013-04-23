function [] = runlga( infile, options )
%RUNLGA Summary of this function goes here
%   Detailed explanation goes here

global proteinroot;

replaceinfile([proteinroot,filesep,'scripts',filesep,'runlga.txt.tpl'],{'#{infile}#'},{infile});
replaceinfile([proteinroot,filesep,'stephan',filesep,'lga',filesep,'lga_bin',filesep,'runlga.sh.tpl'],{'#{options}#'},{options});

old = pwd;
cd(proteinroot);

!runlga

cd(old);

end

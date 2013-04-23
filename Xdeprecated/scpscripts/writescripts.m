function writescripts
%generates script files, depends on set_globals

persistent ampl_set;

if (ampl_set)    
    disp('ampl files already set.');
    if (~ask('set again?'))
        return
    end
end

% get globals previously set
global proteinroot;
global ampl_username;
global tpldir;
global keyfile;
global local_ampldir;
global local_amploutdir;
global scriptdir;
global scriptdir_relative;

global remote_ampldir;
global remote_amploutdir;
global remote_proteinroot;


dprintf('global "ampl_username" is set to "%s".',ampl_username);

% .bats
replaceinfile([tpldir,filesep,'optimize.bat.tpl'],{'#{keyfile}#','#{scriptdir_relative}#'},{keyfile,scriptdir_relative},proteinroot);
replaceinfile([tpldir,filesep,'sync.bat.tpl'],{'#{keyfile}#','#{lscriptdir}#'},{keyfile,scriptdir_relative},proteinroot);

% .scripts
replaceinfile([tpldir,filesep,'data_up.txt.tpl'],{'#{username}#','#{rampldir}#'},{ampl_username,remote_ampldir},scriptdir);
replaceinfile([tpldir,filesep,'data_down.txt.tpl'],{'#{username}#','#{outdir}#'},{ampl_username,remote_amploutdir},scriptdir);
replaceinfile([tpldir,filesep,'execute.txt.tpl'],{'#{username}#','#{rampldir}#'},{ampl_username,remote_ampldir},scriptdir);
replaceinfile([tpldir,filesep,'sync.txt.tpl'],{'#{username}#','#{rprotdir}#'},{ampl_username,remote_proteinroot},scriptdir);

% .knitro script
replaceinfile([tpldir,filesep,'testknitro.tpl'],{'#{username}#','#{rampldir}#'},{ampl_username,remote_ampldir},local_ampldir);
%replaceinfile([ampldir,'test1.run.tpl'],{'#{username}#','#{rampldir}#','#{outdir}#'},{username,rampldir,ramploutoutdir});
    
ampl_set = 1;

end
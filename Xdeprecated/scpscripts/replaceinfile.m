function outfile = replaceinfile(tplfile, tokens, replaces, outputdir)
%outfile = replaceinfile(tplfile, tokens, replaces)
%replaces tokens{i} by replaces{i} in file tplfile and writes to a file
%named tplfile w/o the .tpl extension in (optional) outputdir
   
    [pathstr, name, ext, versn] = fileparts(tplfile);
    
    if (nargin >= 4)
        outfile = [outputdir,filesep,name];
    else
        outfile = [pathstr,filesep,name]; % same dir
    end
    
    if (~strcmp(ext,'.tpl'))
        dprintf('Error: file "%s" has extension "%s", not the required ".tpl"',tplfile,ext);
    else
        dprintf('Writing file "%s" to "%s"',tplfile,outfile);
    end

    fid=fopen(tplfile,'r');
    s = fread(fid, '*char')';
    for i=1:length(tokens),
        s = strrep(s,tokens{i},replaces{i});
    end
    fclose(fid);
       
    fid=fopen(outfile,'w'); 
    fprintf(fid,escapeString(s));
    fclose(fid);
    
end

function str = escapeString(str)
    str = strrep(str,'%','%%');
    str = strrep(str,'\','\\');
    str = strrep(str,'''','''''');    
end
function [ mat ] = pdb2mat( pdbfile, onlycalphas )
%PDB2MAT extract atom information from pdbfile 
%   [ mat ] = pdb2mat( pdbfile )
% mat.atom(1)
%      name: {'N'}
%     amino: {'MET'}
%    aminon: 1
%     coord: [44.0050 -3.1250 9.0380]
%       vol: 1
%     therm: 25.6400
    
if nargin<2, onlycalphas = 0; end
    
fid = fopen(pdbfile,'r');

mat = struct;
mat.amino = struct;
mat.amino.atom = struct;
i = 1;
p = -1;
while 1
    tline = fgetl(fid);
    if ~ischar(tline) || strcmp(tline,''),   break,   end
    if (strcmp(tline(1:6),'HEADER')),
        mat.header = tline;
    end
    
    if (strcmp(tline(1:4),'ATOM')),
        % atom record
        %C = textscan(tline,'ATOM %d8 %s %s %d8 %n %n %n %n %n %*s %*n');        
        atomr = pdbrecord(tline);       

        if (strcmp(atomr.name,'CA') || ~onlycalphas)
            
            if onlycalphas,
                atomr.i = atomr.aanum;
            end
            % save atom record
            mat.atomr(i) = atomr;
            i=i+1;
        end
    end
    %disp(tline)
end
fclose(fid);


end

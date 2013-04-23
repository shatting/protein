function [ pdbr ] = pdbrecord( line )
%PDBRECORD Summary of this function goes here
%   Detailed explanation goes here
% Field |    Column    | FORTRAN |                                         
%   No. |     range    | format  | Description                                   
% ---------------------------------------------------------------------------
%    1. |    1 -  6    |   A6    | Record ID (eg ATOM, HETATM)       
%    2. |    7 - 11    |   I5    | Atom serial number                            
%    -  |   12 - 12    |   1X    | Blank                                         
%    3. |   13 - 16    |   A4    | Atom name (eg " CA " , " ND1")   
%    4. |   17 - 17    |   A1    | Alternative location code (if any)            
%    5. |   18 - 20    |   A3    | Standard 3-letter amino acid code for residue 
%    -  |   21 - 21    |   1X    | Blank                                         
%    6. |   22 - 22    |   A1    | Chain identifier code                         
%    7. |   23 - 26    |   I4    | Residue sequence number                       
%    8. |   27 - 27    |   A1    | Insertion code (if any)                       
%    -  |   28 - 30    |   3X    | Blank                                         
%    9. |   31 - 38    |  F8.3   | Atom's x-coordinate                         
%   10. |   39 - 46    |  F8.3   | Atom's y-coordinate                         
%   11. |   47 - 54    |  F8.3   | Atom's z-coordinate                         
%   12. |   55 - 60    |  F6.2   | Occupancy value for atom                      
%   13. |   61 - 66    |  F6.2   | B-value (thermal factor)                    
%    -  |   67 - 67    |   1X    | Blank                                         
%   14. |   68 - 70    |   I3    | Footnote number                               
% ---------------------------------------------------------------------------
%ATOM      9  N   ASN     2      41.607  -0.820  10.254  1.00 17.71      110L 128
%ATOM      1   CA MET     1       0.000   0.000   0.000  1.00


pdbr.id = strtrim(line(1:6));
pdbr.i = str2double(strtrim(line(7:11)));               %    2. |    7 - 11    |   I5    | Atom serial number                            
pdbr.name = strtrim(line(13:16));                       %    3. |   13 - 16    |   A4    | Atom name (eg " CA " , " ND1")   
pdbr.altloc = strtrim(line(17));                        %    4. |   17 - 17    |   A1    | Alternative location code (if any)            
pdbr.aaname = strtrim(line(18:20));                     %    5. |   18 - 20    |   A3    | Standard 3-letter amino acid code for residue 
pdbr.ch = strtrim(line(22));                            %    6. |   22 - 22    |   A1    | Chain identifier code                         
pdbr.aanum = str2double(strtrim(line(23:26)));          %    7. |   23 - 26    |   I4    | Residue sequence number                       
pdbr.icode = strtrim(line(27));                         %    8. |   27 - 27    |   A1    | Insertion code (if any)                       
pdbr.x = str2double(strtrim(line(31:38)));              %    9. |   31 - 38    |  F8.3   | Atom's x-coordinate                         
pdbr.y = str2double(strtrim(line(39:46)));              %   10. |   39 - 46    |  F8.3   | Atom's y-coordinate                         
pdbr.z = str2double(strtrim(line(47:54)));              %   11. |   47 - 54    |  F8.3   | Atom's z-coordinate                         
pdbr.occ = str2double(strtrim(line(55:60)));            %   12. |   55 - 60    |  F6.2   | Occupancy value for atom                      
pdbr.bval = str2double(strtrim(line(61:66)));           %   13. |   61 - 66    |  F6.2   | B-value (thermal factor)                    
pdbr.footn= strtrim(line(68:70));                       %   14. |   68 - 70    |   I3    | Footnote number 

end

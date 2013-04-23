function [modfile] = writelgainput( mat1, name1 ,mat2, name2, modfile )
%WRITELGAINPUT Summary of this function goes here
%   Detailed explanation goes here

if nargin<5, modfile = [name1, '.', name2]; end

str1 = mat2mol2(mat1,name1);
str2 = mat2mol2(mat2,name2);


fid = fopen(modfile,'w');
fprintf(fid,str1);
fprintf(fid,'\n');
fprintf(fid,str2);
fclose(fid);

global stephanroot;
movefile(modfile,[stephanroot,filesep,'lga',filesep,'lga_bin',filesep,'MOL2']);


end

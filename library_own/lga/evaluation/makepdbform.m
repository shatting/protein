function makepdbform(makestuffup,x,y,seq1,seq2)
% converts the cartesian coordinate chain x to PDB form
% chain can be just Calphas, or N,Ca,C,O,N,Ca,C,O,...
% seq is the atom sequence (for now amino acid sequence, bc
% only calphas are being used)
% x is the model, y is the template (usually, but they can be completely
% different sequences)
% x,y have to be original x/10, since our data bank is 10*calphas
%have to be incredibly careful with spacing-- pdb require exactly 80
%columns

sizex = size(x,1); %number of rows of x -- number of amino acids
sizey = size(y,1); %usually sizex = sizey, since x is the template for y


if nargin<5
    seq2 = seq1;
    if sizex ~= sizey
        error('model and template amino acid chains have different lengths');
    end
end


fid=fopen('pdbform.pdb','w');

makepdbstyle(x,seq1,fid,'model',makestuffup);
makepdbstyle(y,seq2,fid,'target',makestuffup);

fclose(fid);


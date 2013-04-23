%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% aaname3.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [name,letter,polarity,all]=aaname(aa)
% return name and letter of amino acid with label aa
% and polarity according to my classification by 3D analysis
% all contains label, name and letter in the form 'aa=NAM=L:HP'
%
% for polarity scale see hydro.scale
% HH hydrophobic
% H0 hydrophobic but exceptional
% P0,P1 neutral or polar but exceptional
% PP neutral or polar
%       ALA     Alanine          1      A       PP      tiny
%       ARG     Arginine         2      R       P0      long,thin
%       ASN     Asparagine       3      N       PP      short
%       ASP     Aspartic acid    4      D       PP      short
%       CYS     Cysteine         5      C       HH      SH-group,short
%       GLN     Glutamine        6      Q       PP      long,thin
%       GLU     Glutamic acid    7      E       PP      long,thin
%       GLY     Glycine          8      G       PP      no Cb atom
%       HIS     Histidine        9      H       PP      aromatic
%       ILE     Isoleucine      10      I       HH      short
%       LEU     Leucine         11      L       HH      short
%       LYS     Lysine          12      K       PP      long,thin
%       MET     Methionine      13      M       HH      long,thin
%       PHE     Phenylalanine   14      F       HH      aromatic
%       PRO     Proline         15      P       P1      CB-C' bond
%       SER     Serine          16      S       PP      short
%       THR     Threonine       17      T       PP      (short)
%       TRP     Tryptophan      18      W       H0      aromatic
%       TYR     Tyrosine        19      Y       H0      aromatic
%       VAL     Valine          20      V       HH      short
%
function [name,letter,polarity,all]=aaname1(aa)

    if aa== 1, name='Alanine PP';letter='A';polarity='PP';all=' 1=ALA=A:PP';
elseif aa== 2, name='Arginine PO';letter='R';polarity='P0';all=' 2=ARG=R:P0';
elseif aa== 3, name='Asparagine PP';letter='N';polarity='PP';all=' 3=ASN=N:PP';
elseif aa== 4, name='Aspartic acid PP';letter='D';polarity='PP';all=' 4=ASP=D:PP';
elseif aa== 5, name='Cysteine HH';letter='C';polarity='HH';all=' 5=CYS=C:HH';
elseif aa== 6, name='Glutamine PP';letter='Q';polarity='PP';all=' 6=GLN=Q:PP';
elseif aa== 7, name='Glutamic acid PP';letter='E';polarity='PP';all=' 7=GLU=E:PP';
elseif aa== 8, name='Glycine PP';letter='G';polarity='PP';all=' 8=GLY=G:PP';
elseif aa== 9, name='Histidine PP';letter='H';polarity='PP';all=' 9=HIS=H:PP';
elseif aa==10, name='Isoleucine HH';letter='I';polarity='HH';all='10=ILE=I:HH';
elseif aa==11, name='Leucine HH';letter='L';polarity='HH';all='11=LEU=L:HH';
elseif aa==12, name='Lysine PP';letter='K';polarity='PP';all='12=LYS=K:PP';
elseif aa==13, name='Methionine HH';letter='M';polarity='HH';all='13=MET=M:HH';
elseif aa==14, name='Phenylalanine HH';letter='F';polarity='HH';all='14=PHE=F:HH';
elseif aa==15, name='Proline P1';letter='P';polarity='P1';all='15=PRO=P:P1';
elseif aa==16, name='Serine PP';letter='S';polarity='PP';all='16=SER=S:PP';
elseif aa==17, name='Threonine PP';letter='T';polarity='PP';all='17=THR=T:PP';
elseif aa==18, name='Tryptophan HO';letter='W';polarity='H0';all='18=TRP=W:H0';
elseif aa==19, name='Tyrosine HO';letter='Y';polarity='H0';all='19=TYR=Y:H0';
elseif aa==20, name='Valine HH';letter='V';polarity='HH';all='20=VAL=V:HH';
elseif aa==21, name='end';letter='-';polarity='--';all='21=end=-:--';
elseif aa==22, name='AS?';letter='B';polarity='PP';all='22=AS?=B:PP';
elseif aa==23, name='GL?';letter='Z';polarity='PP';all='23=GL?=Z:PP';
                    % GLN or GLU, not GLY
elseif aa==24, name='any';letter='X';polarity='--';all='24=any=X:--';
                    % undetermined
elseif aa==28, name='///';letter='!';polarity='--';all='28=any=!:--';
                    % chain break
else           aa,error('bad amino acid');
end;

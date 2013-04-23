function [Vhydro,fs,Vs] = hydropotential(data,p,pot,wrongseq,tryoutseqs,bond)
% checks to see what the hydropot does, and if it is at all accurate
% so, Vhydro should be lower, if the hydrophobics are on the inside
% and the hydrophilics are on the outside, I think
% really this is just for me to see how it should be typed up before 
% I put it into ampl
% p is the protein, pot is the structure with potential parameters
% wrongseq: calculate Vhydro using different differen aa sequence. bonds stay
%              the same, so Vhydrop should be higher
% be careful that the protein chosen to provide the wrong sequence is
% longer than the actual protein!
% Vhydro is only needed output. fs, Vs give different values depending on
% options
% using tryoutseqs shows that Vhydro is usually smallest for the correct
% sequence, but the difference is very slight. how do I make it bigger?
% bond: find potential from a specific bond matrix

if nargin < 4
    wrongseq = 0;
end
if nargin < 5
    tryoutseqs = 0;
end
if nargin < 3
    load hydropots;
end
if nargin < 6
    bond = 0;
end

Vhydro = 0;

% a function needed to calculate the feature value, made in 
% step1 and describes the average distance of amino acids from
% the center of their protein, based on size of protein
% independent variable is sizex-1
pol=[1.1417e+003,-2.8793e+003, 1.1208e+004];
fit= @(x) pol(1).*x.^(2/3)+pol(2).*x.^(1/3)+pol(3); 
if bond == 0
   bond=double(data{p}.bond);
end

sizebond = size(bond,1);

%try with the wrong protein!! Vhydro should be bigger!
%if wrongstruct > 0
%  bond = double(data{wrongstruct}.bond);
%  bond = bond(1:sizebond,:);
%end

x=calphas(bond)/10;
sizex=size(x,1);
xmean=[0,0,0];
for i=1:sizex
    xmean=xmean+x(i,:);
end
xmean=xmean/sizex;
fitp = fit(sizex - 1);

sequ = data{p}.seq;
if wrongseq > 1
    sequ = data{wrongseq}.seq;
    sequ = sequ(1:sizesex);
end

fs = zeros(sizex,1);
Vs = zeros(sizex,1);



for j = 1:sizex
    aa = sequ(j);
    if aa >= 21
        continue;
    end

    dist = (norm(x(j,:)-xmean))^2; 
    f = dist/fitp; % x is actually the value of feat14(this occurrence of aa in the entire data set)
    fs(j) = f;
    %pot.L = R^(-T)
    Vhydro = Vhydro + (pot(aa).L*(f - pot(aa).mu))'*(pot(aa).L*(f - pot(aa).mu)) - pot(aa).ent/10;    
    Vs(j) = Vhydro;
end

Vhydro = Vhydro * 100000;

if tryoutseqs > 0 %try out some sequences
    Vs = zeros(1,tryoutseqs+1);
    Vs(1) = Vhydro;
    numseq = 0;
    
    while numseq < tryoutseqs
        sequ = data{p}.seq;
        if size(sequ,1) < sizex
            continue;
        end
        numseq = numseq + 1;
        
        for j = 1:sizex
        aa = sequ(j);
          if aa >= 21
            continue;
          end

        dist = (norm(x(j,:)-xmean))^2; 
        f = dist/fitp; % x is actually the value of feat14(this occurrence of aa in the entire data set)
        fs(numseq) = f;
        %pot.L = R^(-T)
        Vhydro = Vhydro + (pot(aa).L*(f - pot(aa).mu))'*(pot(aa).L*(f - pot(aa).mu)) - pot(aa).ent/10;    
        
        Vs(numseq) = Vhydro; %Vs is different now
        
        end %end for j = 1:sized 
    end  % end for p 
end
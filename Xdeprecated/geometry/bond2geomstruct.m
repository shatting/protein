function [ geomstruct ] = bond2geomstruct( bondorchainstruct )
% was bond2geomstruct.m
% [ geomstruct ] = bond2geomstruct( bondorchainstruct )
%   returns a geom struct with geometry info (c,c+,t,t+,..)
%   either augments the chain struct or generates a new one from bond
% [ geomversion ] = bond2geomstruct()
%   returns current version of geometry calculations

if (nargin == 0)
    geomstruct = 1;
    return;
end

if (isfield(bondorchainstruct,'bond'))
    geomstruct = bondorchainstruct;
    bond = geomstruct.bond;
else
    geomstruct = struct;    
    bond = bondorchainstruct;
    geomstruct.bond = bond;
end

naa = size(bond,1)+1;
% global rotational invariants
bb=double(bond);% bond vector sequence
svdbb=svd(bb); % bond matrix singular values

% length and normalized bonds r
len=sqrt(bb(:,1).^2+bb(:,2).^2+bb(:,3).^2); % |b_i|
rr=bb./[len len len];  % r_i

% cos(alpha_i) = r_i^T * r_i+1
ind1=1:naa-2;
ind2=2:naa-1;
cc=rr(ind1,1).*rr(ind2,1)+rr(ind1,2).*rr(ind2,2)...
        +rr(ind1,3).*rr(ind2,3); 

% sin(alpha_i)
ss=sqrt(1-cc.^2); 
sss=ss(1:end-1).*ss(2:end)+realmin; %s_i*s_i+1

% r_i^T * r_i+2
ind3=1:naa-3;
ind4=3:naa-1;
ccc=rr(ind3,1).*rr(ind4,1)+rr(ind3,2).*rr(ind4,2)...
        +rr(ind3,3).*rr(ind4,3); 

% torsions
tcos = cc(ind3).*cc(ind3+1) - ccc; % c_i * c_i+1 - r_i^T * r_i+2^T
tsin = zeros(naa-3,1);
for t=1:naa-3,
    tsin(t)=det(rr(t:t+2,:));
end;
tcos = tcos./sss;
tsin = tsin./sss;
beta = atan2(tsin,tcos);
        
geomstruct.r = rr;
geomstruct.c = cc;
geomstruct.s = ss;
geomstruct.t = tcos;
geomstruct.tp = tsin;
geomstruct.beta = beta;
geomstruct.svdbb = svdbb;

end

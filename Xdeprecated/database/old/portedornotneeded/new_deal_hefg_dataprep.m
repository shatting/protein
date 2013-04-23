% this is now chainstruct_addgeom
function [ dat, naa, res, svdbb, afreq, res24 ] = new_deal_hefg_dataprep( dat )
%HEFG_GEOM2CLASS Summary of this function goes here
%   Detailed explanation goes here

if (nargin<=2),
    mklen=1; mkcos=1; mktors=1; mktorsold=1;
end

res24 = inf;
afreq=zeros(24,1);
% amino acid frequencies
aa=dat.seq;% amino acid sequence
naa=length(aa);
res=dat.res;
% if res(ch)==0, res(ch)=inf;data{ch}.res=inf; end;
for t=1:naa,
    afreq(aa(t))=afreq(aa(t))+1;
    if t>naa-4, break; end;
end;
if max(aa)>23, res24=min(res24,res); end;

% global rotational invariants
bb=double(dat.bond);% bond vector sequence
svdbb=svd(bb); % bond matrix singular values

% create additional information
if mklen | mkcos | mktors,
    % make length data:           len(t)=norm(bond(t,:),2);
    len=sqrt(bb(:,1).^2+bb(:,2).^2+bb(:,3).^2);
    if mklen, dat.len=uint8(round(5*len)); end;
end;
if mkcos | mktors,
    ind1=1:naa-2;
    ind2=2:naa-1;
    rr=bb./[len len len];  % normalized bond vectors
    cc=rr(ind1,1).*rr(ind2,1)+rr(ind1,2).*rr(ind2,2)...
        +rr(ind1,3).*rr(ind2,3);
    if mkcos, dat.cos=int8(round(100*cc)); end;
end;
if mktors || mktorsold,
    ind3=1:naa-3;
    ind4=3:naa-1;
    ccc=rr(ind3,1).*rr(ind4,1)+rr(ind3,2).*rr(ind4,2)...
        +rr(ind3,3).*rr(ind4,3);

    % implicit rotation here! -cos->sin, sin->cos means pi/2, ie.
    % old_angle=new_angle-pi/2
    tsin=ccc-cc(ind3).*cc(ind3+1);
    tcos=zeros(naa-3,1);
    for t=1:naa-3,
        tcos(t)=det(rr(t:t+2,:));
    end;
    
    % get torsion classes
    atsin=abs(tsin);
    atcos=abs(tcos);
    tcl=zeros(naa-3,1);
    tqu=zeros(naa-3,1);
    
    % helix class H=1
    if diagbds, 
        ind=(tcos>=atsin); % [-pi/4,pi/4]
    else
        ind=(tcos>=0 & tsin<=0); % [0,pi/2]
    end
    tcl(ind)=1;
    tqu(ind)=tsin(ind)./atcos(ind); % [-1,1]
    
    % extended class E=2
    if diagbds, 
        ind=(tsin>=atcos); % [pi/4,3 pi/4] 
    else
        ind=(tcos<=0 & tsin>=0); % [pi/2,pi]
    end
    tcl(ind)=2;
    tqu(ind)=tcos(ind)./atsin(ind); % [1,-1]
    
    % irregular class F=3
    if diagbds, 
        ind=(-tcos>=atsin); % [3 pi/4, pi] u [-pi, -3 pi/4]
    else
        ind=(tcos>=0 & tsin>=0); %[-pi,-pi/2]
    end
    tcl(ind)=3;
    tqu(ind)=tsin(ind)./atcos(ind); % [-1,1]
    
    % irregular class G=4
    if diagbds, 
        ind=(-tsin>=atcos); % [-3 pi/4,-pi/4]
    else
        ind=(tcos<=0 & tsin<=0); %[-pi/2,0]
    end
    tcl(ind)=4;
    tqu(ind)=tcos(ind)./atsin(ind); % [1,-1]

    ss = sqrt(1-cc.^2); % 1..naa-2
    ss1 = ss(1:naa-3);
    ss2 = ss(2:naa-2);
    sss = ss1.*ss2;
    % round quotients
    tqu=round(100*tqu);
    dat.tors=int8([tqu,tcl]);
    dat.torsold = [-tsin tcos]; % these are exactly [t,t'] as in dataanalysis pjs
    dat.cosold = cc; % these are exactly the cos(alpha_i)
end

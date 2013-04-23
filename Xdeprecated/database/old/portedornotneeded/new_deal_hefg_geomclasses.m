%dont need that, functionality went into HEFG_classifier
function [ s ] = new_deal_hefg_geomclasses( bond )
%NEW_DEAL_HEFG_GEOMCLASSES Summary of this function goes here
%   Detailed explanation goes here

bond=double(bond);

naa = size(bond,1) + 1;

len=sqrt(bond(:,1).^2+bond(:,2).^2+bond(:,3).^2);
rr=bond./[len len len];  % normalized bond vectors

ind=1:naa-2;
cc=rr(ind,1).*rr(ind+1,1)+rr(ind,2).*rr(ind+1,2)...
        +rr(ind,3).*rr(ind+1,3);

ind=1:naa-3;
ccc=rr(ind,1).*rr(ind + 2,1)+rr(ind,2).*rr(ind+2,2)...
    +rr(ind,3).*rr(ind+2,3);

tsin=ccc-cc(ind).*cc(ind+1);
tcos=zeros(naa-3,1);
for t=1:naa-3,
    tcos(t)=det(rr(t:t+2,:));
end;

% get torsion classes
atsin=abs(tsin);
atcos=abs(tcos);
tcl=zeros(naa-3,1);
%tqu=zeros(naa-3,1);
% helix class H=1
ind=(tcos>=atsin); % [-pi/4,pi/4]
tcl(ind)=1;
%tqu(ind)=tsin(ind)./atcos(ind); % [-1,1]
% extended class E=2
ind=(tsin>=atcos); % [pi/4,3 pi/4]
tcl(ind)=2;
%tqu(ind)=tcos(ind)./atsin(ind); % [1,-1]
% irregular class F=3
ind=(-tcos>=atsin); % [3 pi/4, pi] u [-pi, -3 pi/4]
tcl(ind)=3;
%tqu(ind)=tsin(ind)./atcos(ind); % [-1,1]
% irregular class G=4
ind=(-tsin>=atcos); % [-3 pi/4,-pi/4]
tcl(ind)=4;
%tqu(ind)=tcos(ind)./atsin(ind); % [1,-1]

s = tcl;

end
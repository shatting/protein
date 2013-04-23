function geom = bond2geo(bond)
% special for new HEFG classes

naa = size(bond,1) + 1;
% global rotational invariants
bb=double(bond); % bond vector sequence

% create additional information
len=sqrt(bb(:,1).^2+bb(:,2).^2+bb(:,3).^2);
ind1=1:naa-2;
ind2=2:naa-1;
rr=bb./[len len len];  % normalized bond vectors
cc=rr(ind1,1).*rr(ind2,1)+rr(ind1,2).*rr(ind2,2)... % rr_i'*rr_i+1
  +rr(ind1,3).*rr(ind2,3);

%cos=int8(round(100*cc));

tors=zeros(naa-3,2);
ind1=1:naa-3;
ind2=3:naa-1;
ccc=rr(ind1,1).*rr(ind2,1)+rr(ind1,2).*rr(ind2,2)...
   +rr(ind1,3).*rr(ind2,3);
tsin=ccc-cc(ind1).*cc(ind1+1);  %tsin = t1
tcos=zeros(naa-3,1);   %tcos = t2
for t=1:naa-3,
  tcos(t)=det(rr(t:t+2,:));
end;

atsin=abs(tsin);  %tsin = t1
atcos=abs(tcos);  %tcos = t2
tcl=zeros(naa-3,1);
tqu=zeros(naa-3,1);
% helix class H=1
    ind=(tcos>=atsin); % [-pi/4,pi/4]
    tcl(ind)=1;
    tqu(ind)=tsin(ind);%./atcos(ind); % [-1,1]
% extended class E=2
    ind=(tsin>=atcos); % [pi/4,3 pi/4]
    tcl(ind)=2;
    tqu(ind)=tcos(ind);%./atsin(ind); % [1,-1]
% irregular class F=3
    ind=(-tcos>=atsin); % [3 pi/4, pi] u [-pi, -3 pi/4]
    tcl(ind)=3;
    tqu(ind)=tsin(ind);%./atcos(ind); % [-1,1]
% irregular class G=4
    ind=(-tsin>=atcos); % [-3 pi/4,-pi/4]
    tcl(ind)=4;
    tqu(ind)=tcos(ind);%./atsin(ind); % [1,-1]

ss=sqrt(1-cc.^2);
t0=ss(1:end-1).*ss(2:end)+realmin;
tors(:,1)=tsin./t0;  % don't need this for ampl!!! q is already the sin
tors(:,2)=tcos./t0;

geom.co =  cc(1:naa-3);
geom.cpo = cc(2:naa-2);
geom.to = tors(:,1);
geom.tpo = tors(:,2);
geom.tcl = tcl;
geom.tqu = tqu; %I should be able to find to, tpo from tqu! (using tcl)
%geom.t0 = t0;

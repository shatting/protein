function [ geom ] = bond2geometry( bond )
%BOND2GEOMETRY converts bond sequence to geometry struct
% [ geom ] = bond2geometry( bond )
% input: bond(l,:) = l-th bond vector
% output:        
%       [ notation: b... normalized bond vectors ]
%       [           s... sqrt(1-c.^2)            ]
% geom with
%       c: [362x1 int8]     b_i'*b_i+1 [-100,100], i=1:nfrag = cos(alpha_i)
%      co: [362x1 double]     - ,, -   [-1,1]
%      cp: [362x1 int8]     b_i+1'*b_i+2 [-100,100], i=1:nfrag
%     cpo: [362x1 double]     - ,, -   [-1,1]
%       t: [362x1 int8]     (b_i'*b_i+2 - c_i*c_i+1)/(s_i*s_i+1) = cos(beta_i)
%      to: [362x1 double]     - ,, -
%      tp: [362x1 int8]     det(b_i,b_i+1,b_i+2)/(s_i*s_i+1) = sin(beta_i)
%     tpo: [362x1 double]   
%     phi: [362x1 int8]     = beta_i/pi=atan2(tp,t)/pi   [-1,1]
%    phio: [362x1 double]
%
% note: cp actually means c plus 1, and indeed c(2:end) == cp(1:end-1)
%       whereas tp means tprime = t' and no such thing like for c holds
% note2: *o means *"original", those without o are the integerized normalized
%       numbers
% note3: t is actually -t from the formulas given in the project seminar

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

cos=int8(round(100*cc));

tors=zeros(naa-3,2);
ind1=1:naa-3;
ind2=3:naa-1;
ccc=rr(ind1,1).*rr(ind2,1)+rr(ind1,2).*rr(ind2,2)...
   +rr(ind1,3).*rr(ind2,3);
t1=ccc-cc(ind1).*cc(ind1+1);
t2=zeros(naa-3,1);
for t=1:naa-3,
  t2(t)=det(rr(t:t+2,:));
end;
ss=sqrt(1-cc.^2);
t0=ss(1:end-1).*ss(2:end)+realmin;
tors(:,1)=t1./t0;
tors(:,2)=t2./t0;
tors(:,3)=atan2(tors(:,1),tors(:,2))/pi;

% dont divide by pi since :,1 and :,2 are sinuses not angles
rtors=int8(round(100*tors));

geom.c = cos(1:naa-3); % integer c \in [-100,100]
geom.co = cc(1:naa-3); % original c \in [-1,1]

geom.cp = cos(2:naa-2); 
geom.cpo = cc(2:naa-2);

geom.t = rtors(:,1);
geom.to = tors(:,1);

geom.tp = rtors(:,2);
geom.tpo = tors(:,2);

geom.phi = rtors(:,3); % integer phi \in [-100,100]\pi
geom.phio = tors(:,3)*pi; % original phi in radians \in [-pi,pi]
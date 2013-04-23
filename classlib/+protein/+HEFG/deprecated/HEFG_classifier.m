function [ hefg, tqu ] = HEFG_classifier( geomstructORttpORfragdata, diagbds )
%[  hefg, tqu  ] = HEFG_classifier( geomstructORttp, diagbds ) 
% get HEFG classes.
% INPUT:   
%           either      a geomstruct with at least .t and .tp
%           or          an (nfrag x 2)-matrix containing (sin(beta),cos(beta))
%           (optional)  diagbds
% OUTPUT:  
%           hefg    an (nfrag x 1)-matrix containing (HEFG)
%           tqu     sin/cos quotients
    
if (nargin<2)
    diagbds = 1;
end

if (isstruct(geomstructORttpORfragdata))
    if (isfield(geomstructORttpORfragdata,'featurenames')) % fragdata
       feats = fragdata_getfeaturesbyname(geomstructORttpORfragdata,{'t','tp'});
       tsin = feats(:,2);
       tcos = feats(:,1);
    else
       tsin = geomstructORttpORfragdata.tp; % geomstruct
       tcos = geomstructORttpORfragdata.t;
    end           
else % features
    tsin = geomstructORttpORfragdata(:,2);
    tcos = geomstructORttpORfragdata(:,1);
end

atsin=abs(tsin);
atcos=abs(tcos);

tcl=zeros(size(tsin));
tqu=tcl;

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

% round quotients
tqu=round(100*tqu);

hefg=tcl;

end


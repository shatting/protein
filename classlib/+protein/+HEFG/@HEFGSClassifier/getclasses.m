function [ hefg, tqu ] = getclasses(obj, dataset )
%[  hefg, tqu  ] =  classify(obj, dataset )
% Simple bound classifier for secondary structure HEFG classes.
% INPUT:   
%           dataset with at least {'t','tp'} or {'beta'}
% OUTPUT:  
%           hefg    an (nfrag x 1)-matrix containing (HEFG)
%           tqu     sin/cos quotients
    
diagbds = obj.diagbds;

if (dataset.haveallfeatures({'beta'}))
    tcossin = dataset.getdata({'beta'});
    tcos = cos(tcossin);
    tsin = sin(tcossin);    
else
    tcossin = dataset.getdata({'t','tp'});
    tcos = tcossin(:,1);
    tsin = tcossin(:,2);      
end

atsin=abs(tsin);
atcos=abs(tcos);

hefg=zeros(size(tsin));
%tqu=tcl;

% helix class H=1
if diagbds, 
    ind=(tsin>=atcos); % [pi/4,3 pi/4] 
else
    ind=(tcos>=0 & tsin>=0); % [0,pi/2)
end
hefg(ind)=1;
%tqu(ind)=tsin(ind)./atcos(ind); % [-1,1]

% extended class E=2
if diagbds, 
    ind=(-tcos>=atsin); % [3 pi/4, pi] u [-pi, -3 pi/4]
else
    ind=(tcos<=0 & tsin<=0); %[-pi,-pi/2)
end
hefg(ind)=2;
%tqu(ind)=tcos(ind)./atsin(ind); % [1,-1]

% irregular class F=3
if diagbds, 
    ind=(-tsin>=atcos); % [-3 pi/4,-pi/4]
else
    ind=(tcos<=0 & tsin>=0); % [pi/2,pi)
end
hefg(ind)=3;
%tqu(ind)=tsin(ind)./atcos(ind); % [-1,1]

% irregular class G=4
if diagbds, 
    ind=(tcos>=atsin); % [-pi/4,pi/4]
else
    ind=(tcos>=0 & tsin<=0); %[-pi/2,0]
end
hefg(ind)=4;
%tqu(ind)=tcos(ind)./atsin(ind); % [1,-1]

% round quotients
%tqu=round(100*tqu);

end



function [ obj ] = fromdatabond(obj, bond)
% was bond2geomstruct.m
    % [ geomstruct ] = bond2geomstruct( bondorchainstruct )
    %   returns a geom struct with geometry info (c,c+,t,t+,..)
    %   either augments the chain struct or generates a new one from bond
    % [ geomversion ] = bond2geomstruct()
    %   returns current version of geometry calculations                        

    obj.naa = size(bond,1)+1;
    obj.nfrags = size(bond,1)-2;

    % global rotational invariants
    bb = double(bond);
    obj.svdbb = svd(bb); % bond matrix singular values

    % length and normalized bonds r
    len=sqrt(bb(:,1).^2+bb(:,2).^2+bb(:,3).^2); % |b_i|
    rr=bb./[len len len];  % r_i

    % cos(alpha_i) = r_i^T * r_i+1
    ind1=1:obj.naa-2;
    ind2=2:obj.naa-1;
    cc=rr(ind1,1).*rr(ind2,1)+rr(ind1,2).*rr(ind2,2)...
            +rr(ind1,3).*rr(ind2,3); 

    % sin(alpha_i)
    ss=sqrt(1-cc.^2); 
    sss=ss(1:end-1).*ss(2:end)+realmin; %s_i*s_i+1

    % r_i^T * r_i+2
    ind3=1:obj.naa-3;
    ind4=3:obj.naa-1;
    ccc=rr(ind3,1).*rr(ind4,1)+rr(ind3,2).*rr(ind4,2)...
            +rr(ind3,3).*rr(ind4,3); 

    % torsions
    tcos = cc(ind3).*cc(ind3+1) - ccc; % c_i * c_i+1 - r_i^T * r_i+2^T
    tsin = zeros(obj.nfrags,1);
    for f=1:obj.nfrags,
        tsin(f)=det(rr(f:f+2,:));
    end;
    obj.t = tcos./sss;
    obj.tp = tsin./sss;
    obj.beta = atan2(obj.tp,obj.t);
    
    beta = atan2(tsin,tcos); %% next 3 lines to make HEFG classification consistent between t/tp and beta
    obj.t = cos(beta);
    obj.tp = sin(beta);
    %ATAN2(Y,X) is the four quadrant arctangent of the real parts of the
    %elements of X and Y.

    obj.c = cc;
    obj.s = ss;    
    obj.len = len;

end
function [ sseq, gseq ] = new_deal_hefg_interface( ndstruct, bond, seq )
% [ sbond, gseq ] = new_deal_hefg_interface( ndstruct, bond, seq )
%           get protein classes
% [ conf ] = new_deal_hefg_interface( ndstruct )
%           get only confusion matrix
% [ conf, cov ] = new_deal_hefg_interface( ndstruct, geomdb )
%           get confusion matrix and covariance model (slow)
% IN:
%   ndstruct ... new_deal_data struct, eg. from "load new_deal_ginny"
%   bond     ... bond sequence, eg. from data{1}.bond
%   seq      ... aa sequence, eg. from data{1}.seq
%   data     ... data from RT*.mat or DSbond.mat
% OUT: (see interface.pdf)
%   sseq    ... s sequence
%   gseq     ... \gamma sequence
%   conf(i,j)... C_{S\Gamma}(i,j)
%   cov(k)   ... covariance model for (c,c+,tsincos) and data with
%                sclass=i, gammaclass=j and sys10toX(k,[ns ng])==[i j]
%                tsincos = sin(beta) iff i = 1,3
%                tsincos = cos(beta) iff i = 2,4
%                struct with fields: 
%                [sys10toX(k,[ns ng])==[i j], sysXto10([i j],[ns ng])==k]
%       .freq ... number of elements (=conf(i,j))
%       .mean ... <(c,c+,tsincos)>_D = \mu_{ij}
%       .L    ... L_{ij} = inv(chol(cov)')
%       .ent  ... e_{ij} = log(det(L)) = log(sqrt(det(C^-1)))
%       .cov  ... <(c,c+,tsincos)*(c,c+,tsincos)'>_D = L^-1*L^-T
%  

cli = ndstruct.hefg_info;
% get gamma classes clustering info
if (isfield(cli,'cluster_result')),
    mg = cli.cluster_info.merge(:,1:cli.cluster_result.cut);    
else
    mg = [];
end

if nargin==3, % class sequence for single chain
    % s classes of chain
    bondgeom = bond2geomstruct(bond);
    sseq = HEFG_classifier(bondgeom,true);
    
    % gamma classes of chain
    pot = cli.pot;
    protfeats= cli.featfun(seq2fragseq(seq));

    if (cli.options.autropy),
      pr = cli.ent;
    else      
      gdatam = cli.cl(:,end-1); % use priors of correct iteration
      pr = getfreq(gdatam)/length(gdatam);
    end
    
    %[ gseq] = suff_potvals( pot, protfeats, pr, cli.options ); 
    [ gseq ] = suff_data2v( protfeats, pot, pr, cli.options );
    
    % apply gamma clustering
    gseq = applymerge(gseq,mg);    
else     
    sdata = cli.cl(:,1);
    gdata = cli.cl(:,end);
    % apply gamma clustering
    gdata = applymerge(gdata,mg); 
    
    % conf
    sseq = getfreq([sdata gdata]);
    if nargin == 1 || nargout == 1, return; end % we wanted only confusion matrix
    
    % rename input
    data = bond;
    
    % cov
    ns = max(sdata);
    ng = max(gdata);
    
    %TODO: port from here on
    [frags, tors, cos] = getfrags_hefg(data); % get features
    %frags = [scaled_cos, scaled_cos+, scaled sin/cos or cos/sin, tclass,
    %seq(1:4)], tors = [t=cos(beta), t'=sin(beta)], cos = [c, c+]
    
    % cov with true values..
    %torsfeat = tors;
    
    % cov with class dependent values
    i13 = frags(:,4) == 1 | frags(:,4) == 3;
    i24 = frags(:,4) == 2 | frags(:,4) == 4;
    torsfeat(i13,1) = tors(i13,2); % sin in classes 1,3
    torsfeat(i24,1) = tors(i24,1); % cos in classes 2,4
    
    n = size(frags,1);   
    feats = [cos torsfeat];
    suff = suffstat(ns*ng, zeros(size(feats,2),1));
    suff = suffstat(suff, feats, n, sysXto10([sdata gdata],[ns ng]));

    cov = suff_num2pot(suff,0);

    % rename output
    gseq = cov;

    % here is how its used: +-cos/sin etc..
    if 1,
        i1=frags(:,4) == 1;
        i2=frags(:,4) == 2;
        i3=frags(:,4) == 3;
        i4=frags(:,4) == 4;
                    
        %1
        sind = tors(i1,2);

        angl = asin(sind);  
        cosd = sqrt(1-sind.^2);
        plotsincos(sind,cosd,angl,1);

        %2
        cosd = tors(i2,1);

        angl = acos(cosd);  
        sind = sqrt(1-cosd.^2);
        plotsincos(sind,cosd,angl,2);

        %3
        sind = tors(i3,2);

        angl = asin(sind);  
        cosd = -sqrt(1-sind.^2);
        plotsincos(sind,cosd,angl,3);
        
        %4
        cosd = tors(i4,1);

        angl = acos(cosd);  
        sind = -sqrt(1-cosd.^2);
        plotsincos(sind,cosd,angl,4);
    end
end

    function plotsincos(sin,cos,angl,sub)
        %sin red        
        subplot(1,4,sub);        
        plot(angl,sin,'.r',angl,cos,'.');
    end

end

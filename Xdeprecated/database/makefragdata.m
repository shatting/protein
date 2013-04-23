function [fragdata] = makefragdata(featurenames)
% [fragdata] = makefragdata(featurenames)
% loads geom db and generates a fragdata structure
% INPUT:
%       featurenames (optional) note that the current HEFG and HEL27 classifiers need at least
%                             {'c','cp','t','tp','beta'}. if omitted those
%                             default field names will be used. if =[], all
%                             reasonable fields will be used (defined in
%                             geomstruct_extractfields.m)
% OUTPUT:    
%       fragdata with fields
%           .data   (nfrags x length(featurenames)) features of all fragments from geomdb with aa <= 20
%           .hefg   
%           .hel27
%           .seq
%           .featurenames

loadgeomdb;

if (nargin == 0)
    featurenames = {'c','cp','t','tp','beta'};
end

fragdata = geomdb2fragdata(geomdb,featurenames);
fragdata.hefg = HEFG_classifier(fragdata,true);
fragdata.hel27 = HEL27_classifier(fragdata);
%fragdata.seq = uint8(fragdata.seq);
%fragdata.readme = 'attention! some fields are uint8. convert before use if needed!';

s = whos('bytes','fragdata');
s = s.bytes;
dprintf('\nfragdata is %.1f mb in size, that is %.0f bytes per fragment',s/1024/1024,s/length(fragdata.hefg));



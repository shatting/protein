function [ fielddata, fields ] = geomstruct_getfeaturesbyname( geomstruct, fields, idx )
% [ features, fields ] = geom2features( geomstruct, fields, idx )
% extract fields from geomstruct and returns a matrix containing (all/idx) fragments' 
% features 
% INPUT:    geomstruct   
%           fields       cell array of fieldnames, all fields if omitted or
%                        empty. special treatment of cp/sp: if c is in
%                        fieldnames, cp will be added although it is no
%                        fieldname. same applies to s. this is due to the
%                        fact that c and s dont directly correspond to
%                        fragments (cp and sp need to be shifted by 1)
%           idx          fragment indices
% OUTPUT:
%           fielddata    (nfrags or |idx| x nfields) matrix
%           fields       fields actually used. this is useful if no fields
%                        were given
% uses: (none)


if (nargin < 2 || isempty(fields)), % all fields except name, bond, svdbb
    fields = setdiff(fieldnames(geomstruct),{'name','res','r','bond','svdbb','seq'});
    fields = [fields,'cp','sp'];
    fields = sort(fields);    
end
if (nargin < 3), % all fragments
    % either beta, t or z need to be present    
    if (isfield(geomstruct,'beta'))
        idx = 1:length(geomstruct.beta);
    elseif (isfield(geomstruct,'t'))
        idx = 1:length(geomstruct.t);
    elseif (isfield(geomstruct,'z'))
        idx = 1:length(geomstruct.t);
    else
        error('length could not be determined');
    end
end

nfeats = length(fields);
fielddata = zeros(length(idx),nfeats);
    
for i=1:length(fields),
    if strcmp(fields{i},'cp'),
        temp = geomstruct.c(2:end);
    elseif strcmp(fields{i},'sp'),
        temp = geomstruct.s(2:end);
    else
        temp = geomstruct.(fields{i});
    end
    fielddata(:,i) = temp(idx);
end

end


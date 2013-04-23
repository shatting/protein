function [ output ] = factory( input, x0 )
%FACTORY put in anything, get out the most.
% input:
%   - Dataset
%   - Coords (Chain)
%   - Geometry
%   - bond
%   - coords
%   - chainstruct
% output:
%   - Coords

if (nargin < 2)
    x0 = [0 0 0];
end

output = [];

if (isobject(input))
    cn = class(input);
    if strcmp(cn,'geom.Chain')
        output = geom.Chain(input);
    elseif strcmp(cn,'geom.Coords')
        output = geom.Coords(input.bond,input.x0);
    elseif strcmp(cn,'geom.Geometry')
        output = geom.Coords(geom.coords2bond(input.getcoords),x0);
    elseif strcmp(cn,'geom.Dataset')
        % TODO: implement  
        error('nyi');
        output = geom.Coords(input.getdata({'c','cp','t','tp'}),input.x0);
    else
        
    end
elseif (isstruct(input) && isfield(input,'bond')) %chainstruct
	output = geom.Coords(input.bond,x0);
elseif (size(input,2)==3) % bonds/coords
    if nargin>1 % bond
        output = geom.Coords(input,x0);
    else % coords
        output = geom.Coords(geom.coords2bond(input),input(1,:));
    end
end

if (isempty(output))
	error('geom.factory: type of data could not be determined');
end

end


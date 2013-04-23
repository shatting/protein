function [ rmsds ] = rmsd_progressive( coordsx,coordsy )
%RMSD_PROGRESSIVE Summary of this function goes here
%   Detailed explanation goes here

if(strcmp(class(coordsx),'geom.Coords') || strcmp(class(coordsx),'geom.Chain')),
    coordsx = coordsx.coords;
end
if(strcmp(class(coordsy),'geom.Coords') || strcmp(class(coordsy),'geom.Chain')),
    coordsy = coordsy.coords;
end

for i=1:size(coordsx,1),
    rmsds(i) = sqrt(1/i*sum(sum((coordsx(1:i,:)-coordsy(1:i,:)).^2)));
end

end


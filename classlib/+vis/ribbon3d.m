function [] = ribbon3d( coords, npoints, cl )
%RIBBON3D ribbon3d( coords, npoints )
%   coords(i,:)     ith point
%   npoints         points between data points
%   cl (optional)   fragment classes
n = size(coords,1);

pp = spline(1:n, coords');
interpp = ppval(pp,1:1/npoints:n)';

%figure('Color',[0 0 0]);
plot3(interpp(:,1),interpp(:,2), interpp(:,3),'LineWidth',5);
hold on;
plot3(coords(:,1),coords(:,2), coords(:,3),'or');
axis off;

if (nargin >= 3)
    map = class_colormap(max(cl));
    for i=1:n-3,
       p = (coords(i+1,:)  + coords(i+2,:))/2;       
       text(p(1), p(2), p(3),num2str(cl(i)),'Color',map(cl(i),:));
    end
else
    for i=1:n,
       p = coords(i,:);       
       text(p(1), p(2), p(3),num2str(i));
    end
end

hold off;

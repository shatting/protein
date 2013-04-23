radii = rand(size(torsold,1),1);
radii = sqrt(repmat(radii,1,2));
torsoldr = torsold.*radii;
markersize = 0.05;
colmap = colormap(class_colormap(4));
subx = 2;
suby = 4;

subplot(suby,subx,[1 2 3 4]);
plot(torsoldr(hefg==1,1),torsoldr(hefg==1,2),'.','MarkerSize',markersize,'Color',colmap(1,:))
hold on;
plot(torsoldr(hefg==2,1),torsoldr(hefg==2,2),'.','MarkerSize',markersize,'Color',colmap(2,:))
plot(torsoldr(hefg==3,1),torsoldr(hefg==3,2),'.','MarkerSize',markersize,'Color',colmap(3,:))
plot(torsoldr(hefg==4,1),torsoldr(hefg==4,2),'.','MarkerSize',markersize,'Color',colmap(4,:))
title('cos(\beta)/sin(\beta), radius radomized');
hold off;

subplot(suby,subx,[5]);
class_barhist( atan2(torsold(:,2),torsold(:,1)), hefg, 40, 0, 4, 0, 1, {'H','E','F','G'}, 0);
title('atan2(sin(\beta),cos(\beta))');

subplot(suby,subx,[6]);
%class_barhist( hefg, hefg, 4, 0, 4, 0, 1, {'H','E','F','G'}, 0);
class_barhist(hefg,hefg,1:4,0,4,0,1,[], 0);
title('class counts');

subplot(suby,subx,7);
class_barhist( cosold(:,1), hefg, 20, 0, 4, 0, 1, [], 0);
title('cos(\alpha), absolute');

subplot(suby,subx,8);
class_barhist( cosold(:,1), hefg, 20, 0, 4, 1, 1, [], 0);
title('cos(\alpha), relative');

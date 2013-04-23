cl = hefg_class_info.cl(:,end);
ncl = max(cl);

radii = rand(size(torsold,1),1);
radii = sqrt(repmat(radii,1,2));
torsoldr = torsold.*radii;
markersize = 0.05;
colmap = colormap(class_colormap(ncl));
subx = 2;
suby = 4;


conf = getfreq([hefg cl]);
clab = new_deal_hefg_prednames(conf);

subplot(suby,subx,[1 2 3 4]);
for c = 1:ncl,
    plot(torsoldr(cl==c,1),torsoldr(cl==c,2),'.','MarkerSize',markersize,'Color',colmap(c,:))
    hold on;
end
title('cos(\beta)/sin(\beta), radius radomized');
hold off;

subplot(suby,subx,[5]);
class_barhist( atan2(torsold(:,2),torsold(:,1)), cl, 40, 0, ncl, 0, 1, clab, 0);
title('atan2(sin(\beta),cos(\beta))');

subplot(suby,subx,[6]);
%class_barhist( hefg, hefg, 4, 0, 4, 0, 1, {'H','E','F','G'}, 0);
class_barhist(cl,cl,20,0,ncl,0,1,[], 0);
title('class counts');

subplot(suby,subx,7);
class_barhist( cosold(:,1), cl, 20, 0, ncl, 0, 1, [], 0);
title('cos(\alpha), absolute');

subplot(suby,subx,8);
class_barhist( cosold(:,1), cl, 20, 0, ncl, 1, 1, [], 0);
title('cos(\alpha), relative');

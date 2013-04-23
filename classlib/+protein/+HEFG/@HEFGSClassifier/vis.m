function vis(diagbds)

ttp = data.FeatureDB.db.getdata({'t','tp'});
c = data.FeatureDB.db.getdata({'c'});
cp = data.FeatureDB.db.getdata({'cp'});

sclassifier = class.HEFGSClassifier(diagbds);
hefg = sclassifier.classify(data.FeatureDB.db);

radii = rand(size(ttp,1),1);
radii = sqrt(repmat(radii,1,2));

ttp_scaled = ttp.*radii;

markersize = 0.01;
colmap = colormap(class_colormap(4));
subx = 2;
suby = 2;

figure;
%subplot(suby,subx,[1 2 3 4]);
plot(ttp_scaled(hefg==1,1),ttp_scaled(hefg==1,2),'.','MarkerSize',markersize,'Color',colmap(1,:))
hold on;
plot(ttp_scaled(hefg==2,1),ttp_scaled(hefg==2,2),'.','MarkerSize',markersize,'Color',colmap(2,:))
plot(ttp_scaled(hefg==3,1),ttp_scaled(hefg==3,2),'.','MarkerSize',markersize,'Color',colmap(3,:))
plot(ttp_scaled(hefg==4,1),ttp_scaled(hefg==4,2),'.','MarkerSize',markersize,'Color',colmap(4,:))
%title('\beta, radius radomized');
axis off;
axis equal;
hold off;

figure
subplot(suby,subx,[1]);
vis.class_barhist( data.FeatureDB.db.getdata({'beta'}), hefg, 50, 0, 4, 0, 1, {'H','E','F','G'}, 0);
title('Histogram of \beta');

subplot(suby,subx,[2]);
%class_barhist( hefg, hefg, 4, 0, 4, 0, 1, {'H','E','F','G'}, 0);
vis.class_barhist(hefg,hefg,1:4,0,4,0,1,[], 0);
title('s-class Frequencies');

subplot(suby,subx,3);
vis.class_barhist( c, hefg, 50, 0, 4, 0, 1, [], 0);
title('Histogram of cos(\alpha)');

subplot(suby,subx,4);
vis.class_barhist( c, hefg, 50, 0, 4, 1, 1, [], 0);
title('Scaled Histogram of cos(\alpha)');
axis equal

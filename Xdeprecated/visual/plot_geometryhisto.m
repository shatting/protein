function plot_geometryhisto(geometry);
figure;

subplot(4,2,1);
hist(geometry(:,1),100);
title 'c_i';

subplot(4,2,2);
hist(geometry(:,2),100);
title 'c_{i+1}';

subplot(4,2,3);
hist(geometry(:,3),100);
title 't';

subplot(4,2,4);
hist(geometry(:,4),100);
title 't''';

subplot(4,2,[5 6 7 8]);
scatterplot(geometry(:,[3,4]),1,'','',5000, 0.1);
title 't versus t''';
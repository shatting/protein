function classification_assessment(cl, cs, ts, seq, nclasses, minnperclass)
% classification_assessment(cl, cs, ts, seq, nclasses, minnperclass)
% needed: cl..classes or [target pred]
%         cs = [c cp], ts=[t t'], seq

% if (size(ts,2)==2 && nargin < 5),
%     error('please use classification_assessment(cl, cs, ts, seq) with size(ts,2) == 1 (new t)');
% else
newt = geom.oldt2newt(ts);
    
if (size(cl,2)) == 2,
    target = cl(:,1);
    cl = cl(:,2);
end

n = length(cl);
ncl = max(cl);
% end
for i=1:20,
    xlab(i,:) = protein.utils.aaname(i);
end

if isempty(target)
    for i=1:ncl,
        clab{i} = int2str(i);
    end
else
   conf = suffclass.utils.getfreq([target cl]);
   clab = protein.utils.prednames(conf);
end
    
freq = -sort(-freqs(cl));

if (nargin <5)
    nclasses = ncl;
elseif nargin >= 6,
    nclasses = max(find(freq>minnperclass));
end

geomfeats = [cs ts];

% class hist --------------------------
figure;
vis.class_barhist(cl,cl,1:ncl,0,ncl,0,1,clab);
title('class sizes and colours');


frags = cell(20,20,20,20);
centers = zeros(ncl,2);

% torsion class centers
colmap = colormap(class_colormap(ncl));
figure;
hold on;
for c = 1:ncl,
    cotors = ts(cl==c,:);
    centers = mean(cotors,1);
    plot(centers(1),centers(2),'*','Color',colmap(c,:));
end
title('class centers');
legend(clab);
xlim([-1 1]);
ylim([-1 1]);
hold off;


%torsion frag centers
if 0,
    for i=1:n,
        if mod(i,10000) == 0, disp(i); end;
        frags{seq(i,1),seq(i,2),seq(i,3),seq(i,4)} = [frags{seq(i,1),seq(i,2),seq(i,3),seq(i,4)};[i,geomfeats(i,:),cl(i)]];
    end
    fragtype = sysXto10(seq,[20,20,20,20])';
    fragtypeuq =unique(fragtype); 
    sequq = sys10toX(fragtypeuq,[20 20 20 20]);
    uq = length(fragtypeuq);
    fragcenters = zeros(uq,6); % c, cp, ot, ot', cl, n
    j = 1;

    for i=1:uq,
        if ~mod(i,10000), disp(i); end;
        seqi = sequq(i,:);
        fragcell = frags{seqi(1), seqi(2),seqi(3),seqi(4)};
        
        fragcenters(j,1:4) = mean(fragcell(:,2:5),1);
        fragcenters(j,5) = fragcell(1,6);
        fragcenters(j,6) = size(fragcell,1);
        j=j+1;

    end
    figure;
    scatterplot(fragcenters(:,[3 4]),fragcenters(:,5),'','',100000,0.1)
    title('fragment centers, class color coded');
end

% %scatterplot(data,class,classname,compname,ndraw, markersize, legend);
%r = rand(n,1);
%r = [r r];
%scatterplot((ts+2*(rand(n,1)-0.5) ).*r,cl,'','',1000); % rand works only
%easy with new t

%class_barhist( data, cl, nbins, integerdata, nbiggest, relative )
%class_hist( data, cl, suby, subx, subidx, nbins, xlim, ylim, integerdata,nbiggest )

% figure;
% subplot(3,2,[1 2]);
% class_hist(cs(:,1), cl, 3, 2, [1 2], 200,[],[],0,nclasses);
% title('c histogram, absolute');
% subplot(3,2,[3 4]);
% class_hist(cs(:,2), cl, 3, 2, [3 4], 200,[],[],0,nclasses);
% title('c+ histogram, absolute');
% subplot(3,2,[5 6]);
% class_hist(newt, cl, 3, 2, [5 6], 200,[],[],0,nclasses);
% title('t histogram, absolute');

figure;
subplot(3,2,[1 2]);
vis.class_barhist(cs(:,1), cl, 200,0,nclasses,0,1,clab)
title('c histogram, absolute stacked');
subplot(3,2,[3 4]);
vis.class_barhist(cs(:,2), cl, 200,0,nclasses,0,1,[])
title('c+ histogram, absolute stacked');
subplot(3,2,[5 6]);
vis.class_barhist(newt, cl, 200,0,nclasses,0,1,[])
title('t histogram, absolute stacked');

figure;
subplot(3,2,[1 2]);
vis.class_barhist(cs(:,1), cl, 200,0,nclasses,1,1,clab)
title('c histogram, relative stacked');
subplot(3,2,[3 4]);
vis.class_barhist(cs(:,2), cl, 200,0,nclasses,1,1,[])
title('c+ histogram, relative stacked');
subplot(3,2,[5 6]);
vis.class_barhist(newt, cl, 200,0,nclasses,1,1,[])
title('t histogram, relative stacked');

figure;
subplot(3,2,[1 2]);
vis.class_barhist(cs(:,1), cl, 200,0,nclasses,1,1,clab,1)
title('c histogram, density stretched');
subplot(3,2,[3 4]);
vis.class_barhist(cs(:,2), cl, 200,0,nclasses,1,1,[],1)
title('c+ histogram, density stretched');
subplot(3,2,[5 6]);
vis.class_barhist(newt, cl, 200,0,nclasses,1,1,[],1)
title('t histogram, density stretched');


% plot absolute aa frequency of 10 top classes
figure;
subplot(4,1,1);
vis.class_barhist(seq(:,1), cl, 1:20,0,nclasses,0,0.8,clab)
title('aa_1 class distribution, absolute');
set(gca,'XTick',[1:20]);
set(gca,'XTickLabel',xlab);
set(gca,'xminortick','off');


subplot(4,1,2);
vis.class_barhist(seq(:,2), cl, 1:20,0,nclasses,0,0.8,[])
title('aa_2 class distribution, absolute');
set(gca,'XTick',[1:20]);
set(gca,'XTickLabel',xlab);
set(gca,'xminortick','off');

subplot(4,1,3);
vis.class_barhist(seq(:,3), cl, 1:20,0,nclasses,0,0.8,[])
title('aa_3 class distribution, absolute');
set(gca,'XTick',[1:20]);
set(gca,'XTickLabel',xlab);
set(gca,'xminortick','off');

subplot(4,1,4);
vis.class_barhist(seq(:,4), cl, 1:20,0,nclasses,0,0.8,[])
title('aa_4 class distribution, absolute');
set(gca,'XTick',[1:20]);
set(gca,'XTickLabel',xlab);
set(gca,'xminortick','off');

% -------------------

figure;
subplot(4,1,1);
vis.class_barhist(seq(:,1), cl, 1:20,0,nclasses,1,0.8,clab)
title('aa_1 class distribution, relative');
set(gca,'XTick',[1:20]);
set(gca,'XTickLabel',xlab);
set(gca,'xminortick','off');

subplot(4,1,2);
vis.class_barhist(seq(:,2), cl, 1:20,0,nclasses,1,0.8,[])
title('aa_2 class distribution, relative');
set(gca,'XTick',[1:20]);
set(gca,'XTickLabel',xlab);
set(gca,'xminortick','off');

subplot(4,1,3);
vis.class_barhist(seq(:,3), cl, 1:20,0,nclasses,1,0.8,[])
title('aa_3 class distribution, relative');
set(gca,'XTick',[1:20]);
set(gca,'XTickLabel',xlab);
set(gca,'xminortick','off');

subplot(4,1,4);
vis.class_barhist(seq(:,4), cl, 1:20,0,nclasses,1,0.8,[])
title('aa_4 class distribution, relative');
set(gca,'XTick',[1:20]);
set(gca,'XTickLabel',xlab);
set(gca,'xminortick','off');

end
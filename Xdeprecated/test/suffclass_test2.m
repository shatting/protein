addpath '../suffclass';

tic
numclasses = 4;
numgroupsperclass = 10;
numpointspergroup = 50;
numinitialclassesperclass = 1;

groupsizefactor = 0.5;

% setup classes
classsize = [2,2,1.5,0.1];
classcenters = [[-1,-1];[2,2];[1,-3];[0,0]];

groupcenters = rand(numgroupsperclass,2,numclasses) - 0.5;

for i=1:numclasses,
   groupcenters(:, :, i) = groupcenters(:,:, i)*classsize(i) + repmat(classcenters(i,:),numgroupsperclass,1);
end

p = zeros(numclasses,numgroupsperclass,numpointspergroup,2);
%p = [];

%generate group points
figure(1);
clf;
hold on;
colmap = colormap(hsv(numclasses));
for i=1:numclasses,
    
    classcol = colmap(i,:);
    
    for j=1:numgroupsperclass,        
        p(i,j,:,:) = randn(numpointspergroup,2)*abs(randn)*groupsizefactor + repmat(groupcenters(j,:,i),numpointspergroup,1);
        
        % plot group points
        plot(reshape(p(i,j,:,1),numpointspergroup,1),reshape(p(i,j,:,2),numpointspergroup,1),'.','Color',classcol);       
    end
    
    % plot group centers for class
    plot(groupcenters(:,1,i),groupcenters(:,2,i),'ok');
end
title 'raw data';
hold off;


% create suff stats
suffs = cell(numclasses*numgroupsperclass,1);

emptystat = suffstat(1,zeros(2,1));

index = 1;
for i=1:numclasses,
    for j=1:numgroupsperclass,
        grouppoints = reshape(p(i,j,:,:),numpointspergroup,2);
%         figure(100 + j);
%         hold on;
        suffs{index} = suffstat(emptystat, grouppoints, numpointspergroup, ones(numpointspergroup,1));
        
%         pot = suffpot(suffs{index});      
%         plot(grouppoints(:,1),grouppoints(:,2),'.');        
%         ellipse(pot.mean,1,pot.cov,'-k');
%         hold off;
        index = index +1;
    end
end


classes = zeros(numclasses*numgroupsperclass,1);
% set some initial classes
for i=1:numclasses,
    add = [1:numinitialclassesperclass];
    classes((i-1)*numgroupsperclass + add) = i;    
end

% plot initial assignments
suffclass_test_plot2(2,suffs,classes,0,1);
title 'seed points + group covariance';

% classify
[pot,cl,confus] = suffcov(suffs,classes);

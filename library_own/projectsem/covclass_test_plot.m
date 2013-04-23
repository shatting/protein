function [ ] = covclass_test_plot(i, p, cl, pot )
%COVCLASS_TEST_PLOT Summary of this function goes here
%   Detailed explanation goes here

numpoints = size(p,1) / 2;

class1 = p(cl==1,:);
class2 = p(cl==2,:);

figure(i);
clf;

hold on;
plot(class1(:,1),class1(:,2),'+r')
plot(class2(:,1),class2(:,2),'*g')

proto1 = p(1,:);
proto2 = p(numpoints*2,:);
plot(proto1(1),proto1(2),'or');
plot(proto2(1),proto2(2),'og');

ellipse(pot.mean(:,1),1,pot.cov{1},'r');
ellipse(pot.mean(:,2),1,pot.cov{2},'g');
hold off;
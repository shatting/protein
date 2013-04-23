% trying to find out why circle mean isnt the same as angle mean

bins = -pi:pi/16:pi;

x=randn(100,1);
x=x/max(abs(x))*pi;

subplot(2,2,1);
mx = mean(x)

ct = histc(x,bins);
bar(bins,ct);
hold on;
y = get(gca,'YLim');
plot([mx mx],y*1.1);
hold off;
xlim([-pi,pi]);

x=x+pi;

x(x>pi) = x(x>pi) - 2*pi;
mxshift = mx + pi;
if mxshift > pi, mxshift=mxshift-2*pi; end

mxshift

s = sin(x);
c = cos(x);
ms = mean(s)
mc = mean(c)

subplot(2,2,3);
plot(x,s,'r*');
hold on;
plot([-pi,pi],[ms ms],'r');
plot([-pi,pi],[mc mc],'g');
plot(x,c,'g+');
hold off;
xlim([-pi,pi]);

m = [ms mc];
mnorm = norm(m)
m = m/mnorm;
w = atan2(m(1),m(2))

err = min(norm(mxshift-w),norm(norm(mxshift-w)-2*pi))

subplot(2,2,2);
ct = histc(x,bins);
bar(bins,ct);
hold on;
y = get(gca,'YLim');
plot([mxshift mxshift],y*1.1,'.-');
plot([w w],y*1.1,'*-r');
hold off;
xlim([-pi,pi]);


subplot(2,2,4);
plot(c,s,'.',m(2),m(1),'*r',mc,ms,'or',0,0,'+r');
xlim([-1,1]);
ylim([-1,1]);

function [bond,t,tp,t0,as,bs] = geo2bond(geo)
% special for new HEFG classes

c = [geo.co;geo.cpo(end)];

s = sqrt(1-c.^2);

%need to find to, tpo using tqu and tcl!
t0=s(1:end-1).*s(2:end)+realmin;
t = zeros(size(t0));
tp = zeros(size(t0));
for i = 1:size(geo.co,1)
    q = geo.tqu(i);
    t0h = t0(i);
    if geo.tcl(i) == 1
        % in this case, tqu = tsin = to
        t(i) = q/t0h;
        tp(i) = sqrt(1 - (q/t0h)^2);
        %[t tp t*geo.t0(i) tp*geo.t0(i) q q*geo.t0(i) geo.to(i) geo.tpo(i) geo.to(i)*geo.t0(i) geo.tpo(i)*geo.t0(i)]
    elseif  geo.tcl(i) == 2
        % in this case, tqu = tcos = tpo
         tp(i) = q/t0h;
         t(i) = sqrt(1 - (q/t0h)^2);
        %[t tp t*geo.t0(i) tp*geo.t0(i) q q*geo.t0(i) geo.to(i) geo.tpo(i) geo.to(i)*geo.t0(i) geo.tpo(i)*geo.t0(i)]

    elseif  geo.tcl(i) == 3
        % in this case, tqu = tsin = to 
         t(i) = q/t0h;
         tp(i) = -sqrt(1 - (q/t0h)^2);

        %[t tp t*geo.t0(i) tp*geo.t0(i) q q*geo.t0(i) geo.to(i) geo.tpo(i) geo.to(i)*geo.t0(i) geo.tpo(i)*geo.t0(i)]
    elseif  geo.tcl(i) == 4
        % in this case, tqu = tcos = tpo
         tp(i) = q/t0h;
         t(i) = -sqrt(1 - (q/t0h)^2);

        %[t tp t*geo.t0(i) tp*geo.t0(i) q q*geo.t0(i) geo.to(i) geo.tpo(i) geo.to(i)*geo.t0(i) geo.tpo(i)*geo.t0(i)]
    end
end

bond = zeros(size(c,1),3);

bond(1,1) = 1;
bond(2,1) = c(1);
bond(2,2) = s(1);
as = zeros(3,3,length(c)-1);
bs = zeros(length(c)-1,3);

for i=1:length(c)-1,
    a(1,:) = bond(i,:)*crossmat(bond(i+1,:)'); % r_i x r_{i+1}
% var a11(i in ILIST):= bond2(i)*bond3(i+1)
% var a12(i in ILIST):=  
    a(2,:) = bond(i,:);
    a(3,:) = bond(i+1,:);
    as(:,:,i) = [a(1,:);a(2,:);a(3,:)];
    b = [tp(i)*s(i)*s(i+1); c(i)*c(i+1)+ t(i)*s(i)*s(i+1); c(i+1)];
    bs(i,:) = b;
    bond(i+2,:) = (a\b)';
    bond(i+2,:) = bond(i+2,:)/norm(bond(i+2,:));
end


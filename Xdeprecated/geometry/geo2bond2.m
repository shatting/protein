function [bond,bond2,t2,t,tp2,tp,t0,bs,b2s] = geo2bond2(geometry)
% special for new HEFG classes

c = [geometry.co;geometry.cpo(end)];

s = sqrt(1-c.^2);

%need to find to, tpo using tqu and tcl!
t0=s(1:end-1).*s(2:end)+realmin;  %don't need this -- I'm
%essentially dividing and then multiplying by this (in definition of b), so I can just skip it 
t2 = zeros(size(c,1)-1,1);
tp2 = zeros(size(c,1)-1,1);
t = zeros(size(t0));
tp = zeros(size(t0));

for i = 1:size(geometry.co,1)
    q = geometry.tqu(i);
    t0h = t0(i);
    if geometry.tcl(i) == 1
        % in this case, tqu = tsin = to
        t(i) = q/t0h;
        tp(i) = sqrt(1 - (q/t0h)^2);
        t2(i) = q;
        tp2(i) = sqrt(1 - q^2);
        %[t tp t*geometry.t0(i) tp*geometry.t0(i) q q*geometry.t0(i) geometry.to(i) geometry.tpo(i) geometry.to(i)*geometry.t0(i) geometry.tpo(i)*geometry.t0(i)]
    elseif  geometry.tcl(i) == 2
        % in this case, tqu = tcos = tpo
         tp(i) = q/t0h;
         t(i) = sqrt(1 - (q/t0h)^2);
         tp2(i) = q;
         t2(i) = sqrt(1 - q^2);
        %[t tp t*geometry.t0(i) tp*geometry.t0(i) q q*geometry.t0(i) geometry.to(i) geometry.tpo(i) geometry.to(i)*geometry.t0(i) geometry.tpo(i)*geometry.t0(i)]
        % need to find something equal to ti*t0 (t0(i) = s(i)*s(i+1) without using t0
         [t(i) t0(i) t(i)*t0(i) t2(i) sqrt(1 - q^2)]
    elseif  geometry.tcl(i) == 3
        % in this case, tqu = tsin = to 
        t(i) = q/t0h;
        tp(i) = -sqrt(1 - (q/t0h)^2);
        t2(i) = q;
        tp2(i) = -sqrt(1 - q^2);
        %[t tp t*geometry.t0(i) tp*geometry.t0(i) q q*geometry.t0(i) geometry.to(i) geometry.tpo(i) geometry.to(i)*geometry.t0(i) geometry.tpo(i)*geometry.t0(i)]
    elseif  geometry.tcl(i) == 4
        % in this case, tqu = tcos = tpo
         tp(i) = q/t0h;
         t(i) = -sqrt(1 - (q/t0h)^2);
         tp2(i) = q;
         t2(i) = -sqrt(1 - q^2);
        %[t tp t*geometry.t0(i) tp*geometry.t0(i) q q*geometry.t0(i) geometry.to(i) geometry.tpo(i) geometry.to(i)*geometry.t0(i) geometry.tpo(i)*geometry.t0(i)]
    end
end

bond = zeros(size(c,1),3);
bond2 = zeros(size(c,1),3);
bs = zeros(length(c)-1,3);
b2s = zeros(length(c)-1,3);
%b3s = zeros(length(c)-1,3);

bond(1,1) = 1;
bond(2,1) = c(1);
bond(2,2) = s(1);
bond2(1,1) = 1;
bond2(2,1) = c(1);
bond2(2,2) = s(1);
for i=1:length(c)-1,
    a(1,:) = bond(i,:)*crossmat(bond(i+1,:)'); % r_i x r_{i+1}
% var a11(i in ILIST):= bond2(i)*bond3(i+1)
% var a12(i in ILIST):=  
    a(2,:) = bond(i,:);
    a(3,:) = bond(i+1,:);
    a2(1,:) = bond2(i,:)*crossmat(bond2(i+1,:)'); % r_i x r_{i+1}
    a2(2,:) = bond2(i,:);
    a2(3,:) = bond2(i+1,:);
    b= [tp(i)*s(i)*s(i+1); c(i)*c(i+1)+ t(i)*s(i)*s(i+1); c(i+1)];
    b2 = [tp2(i); c(i)*c(i+1)+ t2(i); c(i+1)];
   % b3 = [tp2(i); c(i)*c(i+1)- t2(i); c(i+1)];
    bs(i,:) = b;
    b2s(i,:) = b2;
    b3s(i,:) = b3;
    bond(i+2,:) = (a\b)';
    bond(i+2,:) = bond(i+2,:)/norm(bond(i+2,:));
    bond2(i+2,:) = (a\b2)';
    bond2(i+2,:) = bond2(i+2,:)/norm(bond2(i+2,:));
end
function d = percentsequenceidentity(a,b)
% finds percent (actually, percent/100) of sequence elements that are
% identical
% not very advanced: if a is shorter than b (or vice versa) truncates b to
% length of a. but doesn't need to be advanced, bc only meant to compare
% class sequences of the same length

la = length(a);
lb = length(b);

if la < lb
   b = b(1:la);
   lb = la;
end
if lb < la
    a = a(1:lb);
    la = lb;
end

u = find(a == b);
d = length(u)/la;
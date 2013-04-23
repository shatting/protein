function HEL  = classify( obj, dataset, idx )
%TODO: generealize this to boundclassifier

if (~dataset.haveallfeatures(obj.requiredfeaturenames))
    error('not all required features present in dataset');
end
if (nargin < 3 || isempty(idx))
    idx = 1:dataset.getndata;
end

geom = dataset.getdata(obj.requiredfeaturenames,idx);

% convert to match old bounds
geom = geom * 100;
geom(:,3) = geom(:,3)/pi;

maxh = 20; % was 20
maxl = 50; % was 50

% -52 - H - 20 ] L - 50 ] E - 100
Hc=(geom(:,1:2)<=maxh);
Ec=(geom(:,1:2)>maxl);
Lc=(geom(:,1:2)>maxh & geom(:,1:2)<=maxl);
Nc=~(Hc|Ec|Lc);
HELc3x3=Hc+2*Lc+3*Ec+100*Nc; %[1..3]x[1..3]
if 0,
% the boundaries for the torsions are currently fixed 
% but should depend on the cosine class
% -100 - E - -85 ] L [ -35 - H - 0 ] E - 100
    Et=geom(:,3)<=-85 | geom(:,3)>0;
    Ht=geom(:,3)>=-35 & geom(:,3)<=0;
    Lt=~(Ht|Et);
    HELt3=Ht+2*Lt+3*Et; %[1..3]
else
   HELt3 = zeros(size(geom,1),1) + 3; % not explicitly assigned class is 3
   t = geom(:,3);
   c = HELc3x3(:,1);
   cp = HELc3x3(:,2);
   
   % last column==inf = 2 classes, else 3 classes last has to be wraparound
   % classes must not be wider than 100
   bds = obj.bounddef;
   
   % for each combination of c, cp classes, assign t class
   for ci=1:3,
    for cip=1:3,
       ind = c==ci & cp==cip;
       bdsrow = bds(cip + 3*ci - 3,:);
       
       HELt3(ind & (t >= bdsrow(1) & t < bdsrow(2))) = 1;
       if (isinf(bdsrow(3)))
           HELt3(ind & (t >= bdsrow(2) | t < bdsrow(1))) = 2;
       else
           HELt3(ind & (t >= bdsrow(2) & t < bdsrow(3))) = 2;
       end
    end
   end
end

HEL=HELc3x3(:,1)+3*HELc3x3(:,2)+9*HELt3-12;

%HEL27_ccpt = [HELc3x3(:,1) HELc3x3(:,2) HELt3];

if (any(HEL>27 | HEL==0)),
    disp('warning, c-outliers, 28 classes in HEL27! assigning class 28.');
    HEL(HEL>27 | HEL==0)=28;
end
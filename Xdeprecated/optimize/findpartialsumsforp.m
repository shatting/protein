function [Vmin,Vminf,U,W,mins] = findpartialsumsforp(p,confus,class_infos,data,newclasses)
%finds the partial minimums for the sums from the left (W) and the sums
%from the right (U)
%I checked, and abs(Vmin - Vminf) < .00000000001 for all the p's (I didn't
%check for smaller numbers, I'm easily satisfied). may 6, 2008
%the tables T and S should actually be saved, and not recalculated each
%time
% mins = [minR,minL] --> to compare the two versions of the minimum
% sequence
% minsR is the minimum sequence calculated from dynoptR
% minsL is the minimum sequence calculated from dynoptL

if ((p == 98) || (p == 99) || (p == 100)) %these proteins have more than 20% aa's greater than 20
    disp('this is a ridiculous protein, filled with unrealistic amino acids');
    Vmin = 0; Vminf = 0; U = 0; W = 0; mins = 0;
    return;
end

confus = confus(:,:,end);

%tables for dynoptR (T) and dynoptL (S)
if newclasses == 0
    a = [1,2,3];
    T = [a;a;a];
    T = [T,T,T];
    %S = [T,T,[a;a;a]];
    S = T;
    %S = [a',a',a',a',a',a',a',a',a'];
    xs = [1 2 3 10 11 12 19 20 21]; %used to convert from S to si (see makeVforp)
else % newclasses = 1
    T = 1; % [1 2 3 4]'; 
    S = T;
end

[Vr,Vl] = makeVforp(data{p}.seq,confus,class_infos,p,newclasses);
[pr,yr,Vmin,U]=dynoptR(T,Vr);
[pl,yl,Vminf,W]=dynoptL(S,Vl);
if newclasses == 0
    minsR = zeros(size(data{p}.seq,1) - 3,1); %minimum sequence using dynoptR
    minsL = zeros(size(data{p}.seq,1) - 3,1); %minimum sequence using dynoptL
    for i = 1:size(data{p}.seq,1)-3
        minsR(i) = pr(i) + 3*yr(i) - 3;
        minsL(i) = 3*pl(i+1) + xs(yl(i)) - 3;
    end
else %newclasses = 1
    minsR = yr'; 
    minsL = yl';
end
mins = [minsR,minsL];
%disp('is the minimum sequence correct? hopefully!');
%disp('this is what the potential of the min. sequence should be:');
Vmin;
%disp('and this is what it is');
V = 0;
if newclasses == 0
    gamma = seq2gammaclasses(data{p}.seq,class_infos);
else
    [sb,gamma] = new_deal_hefg_interface(class_infos,data{p}.bond,data{p}.seq );
    clear sb;  
end
    
for i = 1:size(data{p}.seq,1)-3
    if gamma(i) > 0
       V = V - log(confus(gamma(i),minsR(i)));
    end
end
V;
if abs(V - Vmin) < .00000001
%    disp('woo-hoo!');
end

% disp('min(Vr end)');
% min(Vr(:,:,end),[],2)
% 
% disp('min(Vr 1)');
% min(Vr(:,:,1),[],2)
% 
% disp('min(Vl end)');
% min(Vl(:,:,end),[],2)
% 
% disp('min(Vl 1)');
% min(Vl(:,:,1),[],2)

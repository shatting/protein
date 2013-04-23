function [ g_classes ] = seq2gammaclasses( seq, new_deal_data, opt )
%SEQ2GAMMACLASSES [HEL27] get predicted (gamma-) class sequence of aa sequence seq
% eg. seq2gammaclasses(data{1}.seq, new_deal_data)
% if opt == 1 or opt not given, then 0s are inserted where there are
% fragments f with max(f) > 20

if nargin < 3, opt = 1; end

c_class_info = new_deal_data.c;
cp_class_info = new_deal_data.cp;
final_class_info = new_deal_data.final;

[n, x] = size(seq);

% function [cpred,prob,datapot]=classpot(pot,data,N,pr,options);
%          classify data(1:N,:) using the potential pot and a prior pr  

fragseq = seq2fragseq(seq);

% remove fragments with aa>20
m = max(fragseq,[],2);
fragseq = fragseq(m <= 20,:);
if isempty(fragseq), % data{98,99,100}
    if opt,
        g_classes = zeros(n,1);
    else
        g_classes = [];
    end
    return 
end

% c
pot = c_class_info.pot;
data = c_class_info.featfun(fragseq);
pr = c_class_info.freqreg;
options = c_class_info.options;
c_classes = suff_data2v(data,pot,pr,options)';

% c+1
pot = cp_class_info.pot;
data = cp_class_info.featfun(fragseq);
pr = cp_class_info.freqreg;
options = cp_class_info.options;
cp_classes = suff_data2v(data,pot,pr,options)';

% final
pot = final_class_info.pot;
data = feat_pairs([fragseq c_classes cp_classes]);
pr = final_class_info.freqreg;
options = final_class_info.options;
g_classes = suff_data2v(data,pot,pr,options)';

% reinsert excluded fragments
if (opt),
    rem = find(m>20);
    for i=1:length(rem),
        g_classes = [g_classes(1:rem(i)-1) 0 g_classes(rem(i):end)];
    end    
end

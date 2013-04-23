function [ g_classes ] = getclasses( obj, dataset )
%SEQ2GAMMACLASSES [HEL27] get predicted (gamma-) class sequence of aa sequence seq
% eg. seq2gammaclasses(data{1}.seq, new_deal_data)
% if opt == 1 or opt not given, then 0s are inserted where there are
% fragments f with max(f) > 20

%TODO: the aa>20 thingy
opt = 1;
new_deal_data = obj.potential;

c_class_info = new_deal_data.c;
cp_class_info = new_deal_data.cp;
final_class_info = new_deal_data.final;

n = dataset.getndata;

% remove fragments with aa>20
% m = max(fragseq,[],2);
% fragseq = fragseq(m <= 20,:);
% if isempty(fragseq), % data{98,99,100}
%     if opt,
%         g_classes = zeros(n,1);
%     else
%         g_classes = [];
%     end
%     return 
% end

fragseq = dataset.getdata(obj.requiredfeaturenames);

% c
pot = c_class_info.pot;
data = c_class_info.featfun(fragseq);
options = c_class_info.options;
c_classes = suff_data2v(data,pot,[],options);

% c+1
pot = cp_class_info.pot;
data = cp_class_info.featfun(fragseq);
options = cp_class_info.options;
cp_classes = suff_data2v(data,pot,[],options);

% final
pot = final_class_info.pot;
data = feat_pairs([fragseq c_classes cp_classes]);
options = final_class_info.options;
g_classes = suff_data2v(data,pot,[],options);

% reinsert excluded fragments
% if (opt),
%     rem = find(m>20);
%     for i=1:length(rem),
%         g_classes = [g_classes(1:rem(i)-1) 0 g_classes(rem(i):end)];
%     end    
% end


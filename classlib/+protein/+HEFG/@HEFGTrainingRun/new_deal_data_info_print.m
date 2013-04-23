function [] = new_deal_data_info_print( new_deal_data )
%NEW_DEAL_INFO_INFO Summary of this function goes here
%   Detailed explanation goes here

dprintf('---------------- Potential information --------------------');
if (isfield(new_deal_data,'hefg_info'))
    opt = new_deal_data.hefg_info.options;
    dprintf('Potential generated %s.',new_deal_data.timedate);
    dprintf('HEFG, diagbounds=%i.',new_deal_data.diagbds);
    dprintf('%i of %i iterations',size(new_deal_data.hefg_info.cl,2),opt.max_iterations);
    
    dprintf('term/new/small: %.2f%%/%.2f%%/%.2f%%, autropy: %i, supervised: %i',opt.term_percentage*100,opt.new_classes*100,opt.remove_small*100,opt.autropy,opt.supervised);
    dprintf('%i gamma classes',max(new_deal_data.hefg_info.cl(:,end)));
    dprintf('s-gamma confusion matrix:');
    tightmat(getfreq([new_deal_data.hefg_info.cl(:,1) new_deal_data.hefg_info.cl(:,end)]));
else
    dprintf('HEL27');
end
dprintf('-----------------------------------------------------------');

end


function [] = class_info( c )
%CLASS_INFO Summary of this function goes here
%   Detailed explanation goes here
    
    type = c.type
    %featfun = c.featfun
    ndata = size(c.cl,1)
    ncl_start = max(c.cl(:,1))
    ncl_end = max(c.cl(:,end))
    unlab_1 = c.cl(:,1) == 0;
    nlab_start = sum(c.cl(:,1) > 0)
    nunlab_start = sum(unlab_1)
    dprintf('confus (excl initially unlabelled)')
    tightmat(c.confus)
    if sum(unlab_1) > 0,
      dprintf('initially unlabeled');
      tightmat(getfreq([c.cl(unlab_1,1)+1 c.cl(unlab_1,end)]));
    end
    dprintf('freq')
    tightmat(freqs(c.cl(:,end))')
    niters = size(c.cl,2)    
    time = c.time
    

    options = c.options
end

function [confus, diagp] = showconfus(cl1, cl2, ncl_i, dprintstr)
  confus = suffclass.utils.getfreq([cl1 cl2]);  
  diagp = sum(diag(confus))/sum(confus(:));

  dprintf(dprintstr);
  dprintf('diag=%.2f%%',diagp*100);

  if size(confus,2) < ncl_i,
        confus(1,ncl_i) = 0;
  end      

  suffclass.display.tightmat([[sum(confus(:)) , -1 ,sum(confus,1)] ; zeros(1,size(confus,2)+2)-1 ; [sum(confus,2), zeros(size(confus,1),1)-1,confus]]);

end

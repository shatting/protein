cl_hefg = double(new_deal_data.hefg_info.cl(:,1));
cl_nocl = double(new_deal_data.hefg_info.cl(:,end));
cl_cl = double(new_deal_data.hefg_info.cluster_result.cl);

print = 0;
disp('--------------------------------------------------------');

disp('compare unclustered vs clustered');
q = kldivqualitymeasure(cl_hefg, cl_nocl, print);
dprintf('unclustered: %f',q);

q = kldivqualitymeasure(cl_hefg, cl_cl, print);
dprintf('clustered: %f',q);

disp('--------------------------------------------------------');

disp('compare sorted vs unsorted (unclustered)');
q = kldivqualitymeasure(cl_hefg, cl_nocl, print);
dprintf('both unsorted: %f',q);

q = kldivqualitymeasure(sortcl(cl_hefg), cl_nocl, print);
dprintf('target sorted: %f',q);

q = kldivqualitymeasure(cl_hefg, sortcl(cl_nocl), print);
dprintf('pred sorted: %f',q);

q = kldivqualitymeasure(sortcl(cl_hefg), sortcl(cl_nocl), print);
dprintf('both sorted: %f',q);

disp('--------------------------------------------------------');

disp('compare sorted vs unsorted (clustered)');
q = kldivqualitymeasure(cl_hefg, cl_cl, print);
dprintf('both unsorted: %f',q);

q = kldivqualitymeasure(sortcl(cl_hefg), cl_cl, print);
dprintf('target sorted: %f',q);

q = kldivqualitymeasure(cl_hefg, sortcl(cl_cl), print);
dprintf('pred sorted: %f',q);

q = kldivqualitymeasure(sortcl(cl_hefg), sortcl(cl_cl), print);
dprintf('both sorted: %f',q);

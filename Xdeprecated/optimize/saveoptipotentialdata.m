function saveoptipotentialdata(usehefgclasses,probname)
% [] = saveoptidata(usehefgclasses,problemname)

global classdatadir;

global geomdb;
if (isempty(geomdb))
    loadgeomdb; % geomdb
end

if ~usehefgclasses
  load([classdatadir,filesep,'HEFG_gamma.mat']); % new_deal_data  
else    
  load([classdatadir,filesep,'HEL_gamma.mat']); % new_deal_data
end

% load maxhelix; % avgnumhelices, avglengthhelices are in there
% these depend only on HEL27_anal08.m, so i decided to not put
% them in new_deal
% [avgnumhelices,avglengthhelices] = numheliceslengthn(data);   

if ~usehefgclasses  
  %geom = fragdb2linear(fragdb); % {'cos(i)'  'cos(i+1)'  'tors(i,3)'  'len(i)'  'len(i+1)'  'res'}
  %all_s = HEL27_anal08(geom); % TODO: we get empty classes. obvious, since some ccp combinations only have 2 t classes
  
  HEL27_getfeats; % HEL27feats, hel27
  
  all_s = hel27;
  all_gamma = double(new_deal_data.final.cl(:,end));
  
  z = HEL27_transform2z(HEL27feats,all_s);
  optifeats = [HEL27feats(:,1:2), z];
  potential = opticov( optifeats, all_s, all_gamma );
    
  confus = getfreq([all_gamma,all_s]); % we want gamma in the rows  
else
  %cd ../stephan/allesneu/neum;
  %load RTdata.mat;
  %cd ..
  %feats = getfrags_hefg(data);
  %feats=[cos, cosp, tors, hefgclass, seq(i:i+3)]
  %geom = feats(1:2);
  %z = feats(3);
  %all_s = feats(4);
  %all_gamma = feats(5);
  [confus,potential] = new_deal_hefg_interface(new_deal_data, data );
  confus = confus'; %gamma should be in the rows, HEFG in the columns    
end

save([probname,'_opt_potential.mat'],potential,confus);

%[avgnumhelices,avglengthhelices] = numheliceslengthn(data); 
%class_info.maxhelices = avgnumhelices;

% prediction of Vmin based on protein length, for use in sequencesforp
% made with potdifferences
%class_info.Vbarprediction = @(x) -7.33*x + 58.7 - 35;
%class_info.Vbarprediction = Vbarprediction;

nums = size(confus,2); %max(all_s);
numgamma = size(confus,1); %max(all_gamma);

% gets actual bounds 
%[clower,cupper,cplower,cpupper,zboth,zlower,zupper] = findsbounds(data,usehefgclasses); %removes values of infinity bc they're annoying
%save opt_sbounds clower cupper cplower cpupper zboth zlower zupper;%
% ---> why not use..
% sbounds = struct;
% for s = 1:nums,
%     ind = all_s==s;
%     if ~sum(ind) % empty class. seems useless to set 
%         sbounds.clower(s) = -100;
%         sbounds.cupper(s) = 100;
%         sbounds.cplower(s) = -100;
%         sbounds.cpupper(s) = 100;
%         sbounds.zlower(s) = -inf;
%         sbounds.zupper(s) = inf;
%         sbounds.zboth(s) = inf;
%     else    
%         sbounds.clower(s) = min(geom(ind,1));
%         sbounds.cupper(s) = max(geom(ind,1));
%         sbounds.cplower(s) = min(geom(ind,2));
%         sbounds.cpupper(s) = max(geom(ind,2));
%         sbounds.zlower(s) = min(z(ind));
%         sbounds.zupper(s) = max(z(ind));
%         sbounds.zboth(s) = max(abs(sbounds.zlower(s)),abs(sbounds.zupper(s)));   
%     end
%     % TODO: maybe that c/cp compatibility thing
% end
% save opt_sbounds sbounds;

all_sgamma = sysXto10([all_s,all_gamma],[nums,numgamma]);
sgbounds = struct;
for s = 1:nums*numgamma,
    ind = all_sgamma==s;
    if ~sum(ind) % empty class. seems useless to set 
        sgbounds.clower(s) = -100;
        sgbounds.cupper(s) = 100;
        sgbounds.cplower(s) = -100;
        sgbounds.cpupper(s) = 100;
        sgbounds.zlower(s) = -inf;
        sgbounds.zupper(s) = inf;
        sgbounds.zboth(s) = inf;
    else    
        sgbounds.clower(s) = min(optifeats(ind,1));
        sgbounds.cupper(s) = max(optifeats(ind,1));
        sgbounds.cplower(s) = min(optifeats(ind,2));
        sgbounds.cpupper(s) = max(optifeats(ind,2));
        sgbounds.zlower(s) = min(optifeats(ind,3));
        sgbounds.zupper(s) = max(optifeats(ind,3));
        sgbounds.zboth(s) = max(abs(sgbounds.zlower(s)),abs(sgbounds.zupper(s)));   
    end
end
save([probname,'_opt_sgbounds.mat'],sgbounds);

% find class avgs for initial values
% sclassavgs = classaverages(usehefgclasses,nums,data);
% save opt_sclassavgs sclassavgs;
% ----> why not use the mu's from potential?

%disp('creating potparams.dat --> the potential file for current R, mean, and ent');
%makepotparamsdotdat([local_ampldir,filesep,'_potparams','.dat'], potential,nums,numgammas);
% <--- saved in opt_potential.mat

%individualopt = 0;
%makeboundparamsdotdat([local_ampldir,filesep,'_boundparams','.dat'],nums,individualopt,clower,cupper,cplower,cpupper,zlower,zupper);
% <--- saved in opt_sbounds.mat

% 3-frags!!
%tf = threefrags_interface(0); 
% make potparams for 3frags - called ft2.dat. already made now with only
% 10 classes, bc 400 classes caused the error 'no license'. ugh.
%newpair2cl= condensetf(tf.pair2cl);
%tf.pair2cl = newpair2cl;
%save thfragsstuff tf;
% make3fragparams('_ft2.dat',tf);
%load thfragsstuff;

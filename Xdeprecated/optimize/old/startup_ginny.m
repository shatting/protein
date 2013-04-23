disp('---------- ginny --------------');
disp('startup: adding ginny directories.');
addpath(pwd);

disp('startup: adding globals.');

global ginnyroot;
ginnyroot = pwd;

if ~ask('do ginny stuff?'), return, end;

%global usehefgclasses;
global local_ampldir;

%using new or old classes
usehefgclasses = false;

if (~exist('loadedginnysmats','var'))    
    disp('loading .mat files');
    
    if ~usehefgclasses
      load new_deal_ginny_old; % new_deal_data
      load fragdb;
    else    
      load new_deal_ginny; % new_deal_data
    end
    
    load RT071127; % data

    loadedginnysmats = 1;
   % load maxhelix; % avgnumhelices, avglengthhelices are in there
    % these depend only on HEL27_anal08.m, so i decided to not put
    % them in new_deal
    % [avgnumhelices,avglengthhelices] = numheliceslengthn(data);   
end

global alphas;
global potential;
global class_info; %what is needed to find s and gamma class sequences for a protein

if ~usehefgclasses
  [bds, alphas] = HEL27_definition;    
  geom = fragdb2linear(fragdb);
  all_s = HEL27_anal08(geom);
  all_gamma = double(new_deal_data.final.cl(:,end));
  z = makezforopticov(geom,all_s,alphas);
  potential = opticov( [geom(:,1) geom(:,2), z], all_s, all_gamma );
  class_info = new_deal_data;
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
  [confus,pot] = new_deal_hefg_interface(new_deal_data, data );
  confus = confus'; %gamma should be in the rows, HEFG in the columns  
  potential = pot;
  class_info = new_deal_data;  %this is not actually class_info, but it is what is needed to find s and gamma class sequences for a protein
end

save geompotential potential;

%[avgnumhelices,avglengthhelices] = numheliceslengthn(data); 
%class_info.maxhelices = avgnumhelices;

% prediction of Vmin based on protein length, for use in sequencesforp
% made with potdifferences
%class_info.Vbarprediction = @(x) -7.33*x + 58.7 - 35;
%class_info.Vbarprediction = Vbarprediction;

%potential = opticov(geom(:,1:3),all_s,all_gamma);
if ~usehefgclasses
   confus = getfreq([all_s ,all_gamma]);
   confus = confus';  % my confus has gamma in the rows, stef's is the other way around
end

nums = size(confus,2); %max(all_s);
global numgammas;
numgammas = size(confus,1); %max(all_gamma);

% gets actual bounds 
[clower,cupper,cplower,cpupper,zboth,zlower,zupper] = findsbounds(data,usehefgclasses); %removes values of infinity bc they're annoying
save sbounds clower cupper cplower cpupper zboth zlower zupper;
%load sbounds

% find class avgs for initial values
classavgs = classaverages(usehefgclasses,nums,data);
save classavg classavgs;

disp('creating potparams.dat --> the potential file for current R, mean, and ent');
makepotparamsdotdat([local_ampldir,filesep,'_potparams','.dat'], potential,nums,numgammas);

individualopt = 0;
makeboundparamsdotdat([local_ampldir,filesep,'_boundparams','.dat'],nums,individualopt,clower,cupper,cplower,cpupper,zlower,zupper);

global potparamswzupperlower; % 1 if zupper lower used.
potparamswzupperlower = 1;


% 3-frags!!
%tf = threefrags_interface(0); 
% make potparams for 3frags - called ft2.dat. already made now with only
% 10 classes, bc 400 classes caused the error 'no license'. ugh.
%newpair2cl= condensetf(tf.pair2cl);
%tf.pair2cl = newpair2cl;
%save thfragsstuff tf;
% make3fragparams('_ft2.dat',tf);
load thfragsstuff;


%run some tests
p = 323;
disp(['p = ',num2str(p)]);

if 0,
    disp('run optimization once, to see if everything is working!');
    runopt(data,class_info,confus,'tester',p,1,1000,0)
end

if 1,
    [rmses,bests,bestpre] = optprogram(p,confus,usehefgclasses,class_info,data,'startuptest');
end

% trial of ampl program with all short p's:
if 0,
    shorttrial
    disp('sometimes shorttrial has errors- I think that they are due to proteins');
    disp('with aa > 20 or confus = 0 or other such problems');
end
function [ optidata ] = initialize( optidata, dataset, sclassifier, gclassifier, options )

if (nargin < 5)
    options = struct;
end

if (~isfield(options,'Vmaxpercentile')), options.Vmaxpercentile = 0; end;

tic;

dprintf('generating optimization potential information');

dprintf('%s  getting s classes..',showtime);
optidata.sclass = protein.SClassification(dataset,sclassifier);

dprintf('%s  getting g classes..',showtime);
optidata.gclass = class.Classification(dataset,gclassifier);

optidata.sgconfus = suffclass.utils.getfreq([optidata.gclass.cl, optidata.sclass.cl]); % we want gamma in the rows  

dprintf('%s  getting sg classes..',showtime);
sgclassifier = class.CombinationClassifier(sclassifier,gclassifier);
optidata.sgclass = class.Classification(dataset,sgclassifier);

dprintf('%s  transforming z..',showtime);
z = dataset.ztransform(optidata.sclass);

dprintf('%s  generating sg covariance model..',showtime);      
optidata.sgcovdataset = data.SimpleDataset([dataset.getdata({'c','cp'}) z],{'c','cp','z'});
optidata.sgcovmodel = optidata.sgcovdataset.getcovmodel({'c','cp','z'},optidata.sgclass);

dprintf('%s  finding sg bounds..',showtime);
sgbounds = struct;
for s = 1:optidata.sgclass.classifier.getncl,
    
    ind = optidata.sgclass.cl==s;
    
    if ~sum(ind) % empty class. seems useless to set 
        sgbounds.clower(s) = -1;
        sgbounds.cupper(s) = 1;
        sgbounds.cplower(s) = -1;
        sgbounds.cpupper(s) = 1;
        sgbounds.zlower(s) = -inf;
        sgbounds.zupper(s) = inf;
        sgbounds.zboth(s) = inf;
    else    
        sgbounds.clower(s) = min(optidata.sgcovdataset.getdata({'c'},ind));
        sgbounds.cupper(s) = max(optidata.sgcovdataset.getdata({'c'},ind));
        sgbounds.cplower(s) = min(optidata.sgcovdataset.getdata({'cp'},ind));
        sgbounds.cpupper(s) = max(optidata.sgcovdataset.getdata({'cp'},ind));
        sgbounds.zlower(s) = min(optidata.sgcovdataset.getdata({'z'},ind));
        sgbounds.zupper(s) = max(optidata.sgcovdataset.getdata({'z'},ind));
        sgbounds.zboth(s) = max(abs(sgbounds.zlower(s)),abs(sgbounds.zupper(s)));   
    end
    
end
optidata.sgbounds = sgbounds;

optidata.Vmaxpercentile = options.Vmaxpercentile;
if (options.Vmaxpercentile)
    dprintf('%s  finding sg Vmax..',showtime);
    for i=1:length(optidata.sgcovmodel),
        covi = optidata.sgcovmodel(i);
        indsg = optidata.sgclass.cl == i;

        sgfeats = optidata.sgcovdataset.getdata({'c','cp','z'},indsg);

        sortV = sort(opti.pointpotential(covi,sgfeats'));
        maxidx = ceil(options.Vmaxpercentile * length(sortV));
        optidata.sgmaxV(i) = sortV(maxidx) + covi.ent;
    end
end

dprintf('%s  done!',showtime);

% should save this, as it does not depend on classification
dprintf('%s  finding hydro bounds..',showtime);
%[ pot, chi ] = constr_hydro_getpot();

dprintf('%s  done!',showtime);

% load maxhelix; % avgnumhelices, avglengthhelices are in there
% these depend only on HEL27_anal08.m, so i decided to not put
% them in new_deal
% [avgnumhelices,avglengthhelices] = numheliceslengthn(data);   

%[avgnumhelices,avglengthhelices] = numheliceslengthn(data); 
%class_info.maxhelices = avgnumhelices;

% prediction of Vmin based on protein length, for use in sequencesforp
% made with potdifferences
%class_info.Vbarprediction = @(x) -7.33*x + 58.7 - 35;
%class_info.Vbarprediction = Vbarprediction;

% 3-frags!!
%tf = threefrags_interface(0); 
% make potparams for 3frags - called ft2.dat. already made now with only
% 10 classes, bc 400 classes caused the error 'no license'. ugh.
%newpair2cl= condensetf(tf.pair2cl);
%tf.pair2cl = newpair2cl;
%save thfragsstuff tf;
% make3fragparams('_ft2.dat',tf);
%load thfragsstuff;


%disp('creating potparams.dat --> the potential file for current R, mean, and ent');
%makepotparamsdotdat([local_ampldir,filesep,'_potparams','.dat'], potential,nums,numgammas);
% <--- saved in opt_potential.mat

%individualopt = 0;
%makeboundparamsdotdat([local_ampldir,filesep,'_boundparams','.dat'],nums,individualopt,clower,cupper,cplower,cpupper,zlower,zupper);
% <--- saved in opt_sbounds.mat

end


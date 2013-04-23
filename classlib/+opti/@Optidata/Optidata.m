classdef Optidata < handle
    % opti.Optidata
    %   Container for Optimization data    
    
    properties
        name
        
        sclass
        gclass
        
        sgclass
        sgconfus
        sgbounds
        sgcovdataset
        sgcovmodel
        sgmaxV
        Vmaxpercentile
    end
    
    methods
        
        function obj = Optidata( dataset, sclassifier, gclassifier, problemname, predict_options)            
            
            obj.name = problemname;
            
            if (nargin < 5)
                options = struct;
            end
            
            obj.initialize(dataset, sclassifier, gclassifier, predict_options);
            
        end
        
        function [ this ] = initialize( this, dataset, sclassifier, gclassifier, predict_options )

            if (nargin < 5)
                predict_options = struct;
            end

            if (~isfield(predict_options,'Vmaxpercentile')), predict_options.Vmaxpercentile = 0; end;

            tic;

            dprintf('generating optimization potential information');

            dprintf('%s  getting s classes..',showtime);
            this.sclass = protein.SClassification(dataset,sclassifier);

            dprintf('%s  getting g classes..',showtime);
            this.gclass = class.Classification(dataset,gclassifier);

            this.sgconfus = suffclass.utils.getfreq([this.gclass.cl, this.sclass.cl]); % we want gamma in the rows  

            dprintf('%s  getting sg classes..',showtime);
            sgclassifier = class.CombinationClassifier(sclassifier,gclassifier);
            this.sgclass = class.Classification(dataset, sgclassifier);

            dprintf('%s  transforming z..',showtime);
            z = dataset.ztransform(this.sclass);

            dprintf('%s  generating sg covariance model..',showtime);      
            this.sgcovdataset = data.SimpleDataset([dataset.getdata({'c','cp'}) z],{'c','cp','z'});
            this.sgcovmodel = this.sgcovdataset.getcovmodel( [], this.sgclass, predict_options.use_sgmatcov );

            dprintf('%s  finding sg bounds..',showtime);
            sgbounds = struct;
            for s = 1:this.sgclass.classifier.getncl,

                ind = this.sgclass.cl==s;

                if ~sum(ind) % empty class. seems useless to set 
                    sgbounds.clower(s) = -1;
                    sgbounds.cupper(s) = 1;
                    sgbounds.cplower(s) = -1;
                    sgbounds.cpupper(s) = 1;
                    sgbounds.zlower(s) = -inf;
                    sgbounds.zupper(s) = inf;
                    sgbounds.zboth(s) = inf;
                else    
                    sgbounds.clower(s) = min(this.sgcovdataset.getdata({'c'},ind));
                    sgbounds.cupper(s) = max(this.sgcovdataset.getdata({'c'},ind));
                    sgbounds.cplower(s) = min(this.sgcovdataset.getdata({'cp'},ind));
                    sgbounds.cpupper(s) = max(this.sgcovdataset.getdata({'cp'},ind));
                    sgbounds.zlower(s) = min(this.sgcovdataset.getdata({'z'},ind));
                    sgbounds.zupper(s) = max(this.sgcovdataset.getdata({'z'},ind));
                    sgbounds.zboth(s) = max(abs(sgbounds.zlower(s)),abs(sgbounds.zupper(s)));   
                end

            end
            this.sgbounds = sgbounds;

            this.Vmaxpercentile = predict_options.Vmaxpercentile;
            if (predict_options.Vmaxpercentile)
                dprintf('%s  finding sg Vmax..',showtime);
                for i=1:length(this.sgcovmodel),
                    covi = this.sgcovmodel(i);
                    indsg = this.sgclass.cl == i;

                    sgfeats = this.sgcovdataset.getdata({'c','cp','z'},indsg);
                    % TODO: dont use pointpotential here, use
                    % potential_values
                    sortV = sort(opti.pointpotential(covi,sgfeats'));
                    maxidx = ceil(predict_options.Vmaxpercentile * length(sortV));
                    this.sgmaxV(i) = sortV(maxidx) + covi.ent;
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


    end
    
end


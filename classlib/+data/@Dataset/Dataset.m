classdef Dataset < handle
    % data.Dataset < handle
    %   Lazy evaluation data carrier interface. Subclasses provide data
    %   through implementing initialize().
    %
    % PROPERTIES
    %   .featurenames*  [1 x n_features] string cell holding feature ids
    %   .data*          [n x n_features] data matrix
    % METHODS
    %       display()
    %       clone = .getdatasetclone(names,idx)
    %           returns a cloned dataset of type data.SimpleDataset, 
    %           optionally filtered with a string cell 'names' and an data 
    %           entry index 'idx'.
    %       getdata()
    %           initializes the dataset, if necessary
    %       [odata,names] = getdata(names,idx)
    %           returns the properties, optionally filtered    
    %       featurenames = getfeaturenames()
    %       n = getnfeatures()
    %       n = getndata()
    %       b = haveallfeatures(names)
    %           returns 1 if all features listed in names are present        
    %       reorder(neworder)
    %           reorders features according to stringcell 'neworder'
    %       appendfeature(fdata,fname)
    %           adds a new feature
    %       removefeature(fname,dontthrowerror)
    %           removes feature, throws error iff feature not present
    %
    %       covmodel = getcovmodel({featurenames} | [], class.Classification)
    %           builds a covariance model from 'featurenames' or 
    %           all features if featurenames = []
    %           
    % Known subclasses:
    %   data.FeatureDB
    %   data.FeatureDBDataset
    %   data.SimpleDataset
    %   geom.Chain
    %   geom.Geometry
    
    properties (SetAccess=private, GetAccess=private) 
        % we dont want subclasses messing around with these
        data 
        featurenames        
    end
    
    methods
        
        function display(obj)                                    
            mc = metaclass(obj);
            % find out if we are dealing with a cell or an array
            datas = {obj.data};
            fnames = {obj.featurenames};
            for i=1:length(datas),
                if (isempty(datas{i}))
                    dprintf('Empty or uninitialized dataset, instance of %s.',mc.Name);
                else                
                    l = min(size(datas{i},1),5);            
                    dprintf('   %s    ',stringcell2string(fnames{i},'      '));
                    suffclass.display.tightmat(datas{i}(1:l,:),' % .3f  ',' % .f  ');                
                    dprintf('Dataset with %i items and %i features, instance of %s.',size(datas{i},1),length(fnames{i}),mc.Name);                
                    dprintf('{%s}',stringcell2string(fnames{i},', '));
                end
            end
        end
        
        function clone = getdatasetclone(obj,names,idx)
            if (nargin == 2) % obj.getdatasetclone(names)
                [d,n] = obj.getdata(names); % do it this way to get names iff == []
                clone = data.SimpleDataset(d,n);
            elseif (nargin == 3) % obj.getdatasetclone(names,idx)
                [d,n] = obj.getdata(names,idx); % do it this way to get names iff == []
                clone = data.SimpleDataset(d,n);
            else % obj.getdatasetclone()
                clone = data.SimpleDataset(obj.getdata,obj.featurenames);
            end
        end  
        
        function [odata,names] = getdata(dataset,names,idx)
            if (isempty(dataset.data))
                dataset.initialize();
            end
                        
            if (nargout > 0)
                if (nargin == 3)
                    [odata,names] = dataset.getdatabynameandidx(names,idx);
                elseif (nargin == 2)
                    [odata,names] = dataset.getdatabynameandidx(names);
                else
                    odata = dataset.data;
                    names = dataset.featurenames;
                end
            end
        end
                
        function value = getfeaturenames(dataset)
            dataset.getdata();                                               
            value = dataset.featurenames;            
        end
        
        function n = getnfeatures(dataset)
            dataset.getdata();
            n = length(dataset.featurenames);
        end
        
        function n = getndata(dataset)
            dataset.getdata();
            n = size(dataset.data,1);
        end               
        
        function b = haveallfeatures(dataset,names)
            % initialize
            dataset.getdata();
            
            b = true;
            for i=1:length(names),
                b = b && (stringcellindexof(dataset.featurenames,names{i}) > 0);
            end
        end
        
        function reorder(obj,neworder)
            obj.data = obj.getdatabynameandidx(neworder);
            obj.featurenames = neworder;
        end
        
        function appendfeature(obj,fdata,fname)
            if (size(fdata,1) ~= size(obj.data,1))
                error('added feature data length does not match dataset size.');
            end            
            if (~ischar(fname))
                error('feature name supplied is not a char array');
            end
            if (obj.haveallfeatures({fname})),                              
                error('feature named "%s" already present in dataset.',fname);                
            end
            obj.data(:,end+1) = fdata;
            obj.featurenames = [obj.featurenames,fname];
        end
        
        function removefeature(obj,fname,dontthrowerror)
            if (obj.haveallfeatures({fname}))      
                remidx = stringcellindexof(obj.featurenames,fname);
                obj.data(:,remidx) = [];
                obj.featurenames = stringcellremove(obj.featurenames,remidx); 
            else
                if nargin<3 || ~dontthrowerror
                    error('feature named "%s" is not in dataset.',fname);
                end
            end
        end
                        
        function covmodel = getcovmodel(obj, featurenames, classification, debug_usematcov)
            
            if (isempty(featurenames)),
                featurenames = obj.getfeaturenames();
            end
            
            covdata = double(obj.getdata(featurenames));
            nfeats = length(featurenames);
            
            % make numerical suff
            suff = suffclass.suffstat(classification.classifier.getncl, zeros(nfeats,1));
            
            opts = struct;
            if (debug_usematcov)                
                opts.debug_usematcov = 1;
            else                
                opts.debug_usematcov = 0;
                suff.adddata( covdata, classification.cl );
            end
                        
            opts.num_reg_weight = 0;
            
            pot = suffclass.potential(suff, opts, covdata, classification.cl);            

            % remove regularization? does not really make sense here            
            covmodel = pot.numpot;
            
        end
        
        function z = ztransform(obj,sclassification,inplace)
            % z_i = tan((beta_i - alpha_s{i})/2)            
            if (isa(sclassification,'protein.SClassifier'))
                cl = sclassification.get_classes(obj);
                betabar = sclassification.gettorsioncenters;
            elseif (isa(sclassification,'class.Classification'))
                cl = sclassification.cl;
                betabar = sclassification.classifier.gettorsioncenters;
            else
                error('could not determine type of input.');
            end
            
            beta = obj.getdata({'beta'});
             
            if (size(betabar,1) < size(betabar,2))
                betabar = betabar';
            end
            
            betabarcl = betabar(cl);
            
            minus = beta - betabarcl;                        
            % minus(minus<-pi/4) = minus(minus<-pi/4)+2*pi; does not change
            % the zs
%             if (max(abs(minus)) > pi/4)
%                 figure;
%                 hist(minus,100);
%                 idx = find(abs(minus)>pi/4)
%                 beta(idx)
%                 sclassification.cl(idx)
%                 %error('too large angles found.');
%             end
            
            z = tan( minus/2 );   
            
            if (nargin > 2 && inplace)
                obj.appendfeature(z,'z');
                obj.appendfeature(cl,'clsz');
            end
        end
        
        % beta = betatransform(obj,sseq,sclassifier)
        % beta = betatransform(obj,sclassifier,inplace)            
        function beta = betatransform(obj,sseq,sclassifier)            
            % beta_i = 2*atan(z_i) + alpha_s{i};
            inplace = 0;
            if (isa(sseq,'protein.SClassifier'))
                if (nargin > 2)
                    inplace = sclassifier;
                end
                sclassifier = sseq;
                sseq = obj.getdata({'clsz'});
            end
            z = obj.getdata({'z'});
            
            betabar = sclassifier.gettorsioncenters;
            
            if (size(betabar,1) < size(betabar,2))
                betabar = betabar';
            end            
            
            betabarcl = betabar(sseq);
            beta = 2*atan(z) + betabarcl; 
            beta(beta>pi) = beta(beta>pi)-2*pi; % just to be correct
            beta(beta<-pi) = beta(beta<-pi)+2*pi; % just to be correct
            
            if (inplace)
                obj.appendfeature(beta,'betafromz');
            end
        end
        
        % should be faster than betatransform
        function [t,tp] = ttptransform_rational(obj,sseq,sclassifier)
            inplace = 0;
            if (isa(sseq,'protein.SClassifier'))
                if (nargin > 2)
                    inplace = sclassifier;
                end
                sclassifier = sseq;
                sseq = obj.getdata({'clsz'});
            end
            
            betabar = sclassifier.gettorsioncenters;            
            if (size(betabar,1) < size(betabar,2))
                betabar = betabar';
            end                        
            sinbetabar = sin(betabar);
            cosbetabar = cos(betabar);                        
            z = obj.getdata({'z'});      
            
            [t,tp] = geom.ttptransform_rational(z, sseq, cosbetabar, sinbetabar);
            
            if (inplace)
                obj.appendfeature(t,'tfromz');
                obj.appendfeature(tp,'tpfromz');
            end
        end
     
        % [x] = getknitrofeatures(obj) uses default features {'c','cp','z'}
        % [x] = getknitrofeatures(obj,zfeaturetouse)
        %           zfeaturetouse .. feature name of z feature
        function [x] = getknitrofeatures(obj,zfeaturetouse)
            if (nargin < 2), zfeaturetouse = {'z'}; end;
            
            featurestouse = [{'c','cp'},zfeaturetouse];
            if (~obj.haveallfeatures(featurestouse))
                error('not all features required for knitro are present in dataset:%s',stringcell2string(featurestouse));
            end
            f = obj.getdata(featurestouse);                        
            
            x = [f(:,1);f(end,2);f(:,3)];                
        end
        
        % 
        function [havefiltered,filteridx] = filter(obj,filternames,filterfn)
                        
            filterdata = obj.getdata(filternames);
            filteridx = filterfn(filterdata);            
            havefiltered = any(filteridx==0);
            
            % fix a bug introduced through inter-item dependencies on c-cp
            % if we remove one or multiple consecutive items in 
            % data, cp of the previous
            % item is not the same as c in the next, as would be required
            % to make it compatible to getknitrofeatures, so we replace the
            % previous items cp with the followings c (not the other way
            % around, which is equally reasonable).
            % TODO: PROBLEM: right after finishing it, i occured to me that 
            % all this code isnt needed: just rewrite the whole cp from the 
            % filtered c.. but oh well it works so im leaving it here for now
            
            % we dont need at least 2 items left
            remainingndata = size(obj.data,1)-sum(filteridx==0);
            % we also dont have to adjust first or last
            remainingndatatoadjust =  sum(filteridx==0) - (filteridx(1)==0) - (filteridx(end)==0);                        
            
            if havefiltered && remainingndata>=2 && remainingndatatoadjust > 0 ...
                    && ismember('c',obj.featurenames) && ismember('cp',obj.featurenames)
                
                cidx = stringcellindexof(obj.featurenames,'c');
                cpidx = stringcellindexof(obj.featurenames,'cp');
                csourceidx = find(~filteridx)+1; % take c from following
                cptargetidx = find(~filteridx)-1; % put into previous
                
                % no source/ target if target were out of bounds
                if (cptargetidx(1) == 0),
                    csourceidx = csourceidx(2:end);
                    cptargetidx = cptargetidx(2:end);
                end
                % same for source
                if (csourceidx(end) > size(obj.data,1))
                    csourceidx = csourceidx(1:end-1);
                    cptargetidx = cptargetidx(1:end-1);
                end
                
                % now look for consecutive filtered items
                % we only want the highest sources and the lowest targets
                % TODO: do this nifty, without for loop, somehow
                csourceidxtemp = csourceidx;
                cptargetidxtemp = cptargetidx;
                for i=length(csourceidx)-1:-1:1,
                    if (csourceidx(i+1) == csourceidx(i)+1)
                        csourceidxtemp(i) = csourceidxtemp(i+1);
                    end
                    if (cptargetidx(end-i) == cptargetidx(end-i+1)-1)
                        cptargetidxtemp(end-i+1) = cptargetidxtemp(end-i);
                    end
                end
                csourceidx = unique(csourceidxtemp); % gets sorted which is no prob
                cptargetidx = unique(cptargetidxtemp);
                
                % now copy sources to targets
                obj.data(cptargetidx,cpidx) = obj.data(csourceidx,cidx);
                
                % and filter. 
                obj.data = obj.data(filteridx,:);
            else
                obj.data = obj.data(filteridx,:);
            end                        
            
            filteridx = ~filteridx; % seems more natural
        end
        
    end

    methods (Access=protected)
        
        % callback from initialize
        function initialized(dataset,data,featurenames)
            dataset.data = data;
            dataset.featurenames = featurenames;
        end
    end
    
    methods (Access=private)
        % this returns an (|names| x |idx|) matrix containing the features
        % both names and idx may be empty, i which case all possible values
        % will be returned
        function [feats,names] = getdatabynameandidx(dataset,names,idx)
            % initialize
            dataset.getdata();
            
            if (nargin < 2 || isempty(names))
                names = dataset.featurenames;
            end
            
            if (nargin < 3 || isempty(idx))
                idx = 1:size(dataset.data,1);
            end
            
            if (~iscell(names)) % 'featurename' instead of {'featurename'}
                names = {names};
            end
            
            idxf = zeros(1,length(names));
            for i=1:length(names),
                iii = stringcellindexof(dataset.featurenames,names{i});
                if iii<1
                    error('dataset does not have requested feature: %s',names{i});
                end
                idxf(i) = iii;
            end

            feats = dataset.data(idx,idxf);
        end
        
    end
           
    % subclasses must implement this to fill the .data property and call
    % initialized. if featurenames is given, data must be re-initialized
    % with the new features
    % dataset.initialize() fills the dataset with all possible features    
	methods (Abstract,Access=protected)        
        initialize(dataset)        
    end
        
end


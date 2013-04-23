classdef FeatureDB < data.Dataset
    % data.FeatureDB < data.Dataset
    %
    % Singleton class providing all data from data.GeomDB as one Dataset.
    %
    % METHODS
    %   see data.Dataset
    
    %% Singleton behaviour
    properties (GetAccess=private, Constant)
        dbsing = data.FeatureDB();
    end
        
    methods(Static=true)
        function outdb = db()
            outdb = data.FeatureDB.dbsing;
        end
    end
    %% Properties
    properties        
                    
    end
    
    %% Methods
    
    methods (Access=protected)
        % implement abstract dataset method and make it private
        % [ fragdata ] = geomdb2features( geomdb, featfun )
        % split up chains into fragments and store into this data.Dataset
        function [ obj ] = initialize( obj  )        
                                                  
            geomdb = data.GeomDB.db;
            dprintf('collecting fragment features from %i fragments',geomdb.nfrag);
            % TODO: take smallest common subset of features over all chains
            fieldnames = geomdb.chains{1}.getfeaturenames();

            features = zeros(geomdb.nfrag,length(fieldnames)); 

            count = 0;

            for ch=1:geomdb.nch,

                if (mod(ch,100) == 0),
                    disp(sprintf('collecting features from chain %i',ch));
                end

                chain = geomdb.chains{ch};     

                fragseq = chain.getdata({'aa1','aa2','aa3','aa4'});

                maxfragaa = max(fragseq,[],2);
                idx = maxfragaa<=20;
                nidx = find(idx); % 

                if (nidx)
                    nchfrag = double(sum(idx));
                    % transpose into frag index
                    fidx = count+1:count+nchfrag;

                    features(fidx,:) = chain.getdata(fieldnames, nidx );        

                else
                   dprintf('skipped chain %i, no fragments <= 20', ch);
                end
                count = max(fidx);

            end %end for ch = 1:datasize

            % truncate
            features(max(fidx)+1:end,:) = [];

            obj.initialized(features, fieldnames);

            end
    end
    
    methods (Access=private)       
        % make constructor private, too
        function obj = FeatureDB()
            obj = initialize(obj); % fill with all possible features
        end       
    end        
            
end


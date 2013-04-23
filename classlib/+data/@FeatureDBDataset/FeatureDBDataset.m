classdef FeatureDBDataset < data.Dataset
    % data.FeatureDBDataset < data.Dataset 
    %   Provides a view (features,dbidx) of data.FeatureDB.
    %
    % currently unused.
    
    properties
        featuredbnames
        featuredbindex % from where the data comes from in featuredb
    end
    
    methods (Access = protected)
        function initialize(obj)           
            fdata = data.FeatureDB.db().getdata(obj.featuredbnames,obj.featuredbindex);
            
            obj.initialized(fdata,obj.featuredbnames);
        end 
    end
    
    methods

        % construct an instance of FeatureDBDataset, containing the
        % features described by featurenames and indices from dbindex
        function obj = FeatureDBDataset(featurenames,dbidx)
            if (nargin < 1 || isempty(featurenames))                
                obj.featuredbnames = data.FeatureDB.db().getfeaturenames;
            else
                obj.featuredbnames = featurenames;
            end  
            
            if (nargin < 2 || isempty(dbidx))                
                obj.featuredbindex = 1:data.FeatureDB.db().getndata;
            else
                obj.featuredbindex = dbidx;
            end                                        
        end
        
    end
    
end


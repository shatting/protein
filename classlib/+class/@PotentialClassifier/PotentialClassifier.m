classdef PotentialClassifier < class.IClassifier
    % +class@IClassifier
    %   abstract superclass for classifiers
    %
    %   known subclasses:
    %       +protein+HEFG@HEFGGammaClassifier
    %
    % PROPERTIES
    %   .potential (+suffclass@potential)
    %
    % CONSTRUCTORS
    %   obj = PotentialClassifier(potential)
        
    properties
        potential_values
    end
        
    methods
        function this = PotentialClassifier(featfun, potential_values)
            this = this@class.IClassifier(featfun);
            this.potential_values = potential_values;
        end
        
    end
    
    methods (Access=protected)
        function [ cl ] = getclasses( this, dataset )
            if (length(this.requiredfeaturenames) > 1)
                error('only one feature set allowed in potential classifiers');
            end
            feats = this.feature_fn(dataset.getdata(this.requiredfeaturenames{1}));

            cl = this.potential_values.calcvals(feats,this.potential_values.priorsORentropies,this.potential_values.ispriors);

        end
        
    end
    
end


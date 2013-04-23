classdef SClassifier < class.BoundClassifier
    % protein.SClassifier < class.BoundClassifier
    %   Abstract superclass for bound based s-classifiers
    %
    %   known subclasses: 
    %       protein.HEFG.HEFGSClassifier
    %       protein.HEL27.HEL27SClassifier
    %
    % PROPERTIES
    %
    % CONSTRUCTORS
    %   obj = SClassifier(bounds)
    %
    % ABSTRACT METHODS
    %   tc = gettorsioncenters(obj)
    %   distm = getdistm(obj)
    %
    % METHODS
    %
    
    properties        
    end
    
    methods
        function obj = SClassifier(bounds)
           obj = obj@class.BoundClassifier(bounds);           
        end                        
    end
    
    methods (Abstract)
        tc = gettorsioncenters(obj)
        distm = getdistm(obj)
    end
    
end


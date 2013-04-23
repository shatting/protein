classdef BoundClassifier < class.IClassifier    
    % +class@BoundClassifier
    %   abstract superclass for simple bound classifiers
    %
    %   known subclasses:
    %       +protein@SClassifier
    %
    % PROPERTIES
    %   .bounddef
    %       
    % CONSTRUCTORS
    %   obj = BoundClassifier(bounddef)
    %
    properties
        bounddef %
    end
    
    methods
        
        function obj = BoundClassifier(bounddef)
            obj.bounddef = bounddef;
        end                       
    end
    
end


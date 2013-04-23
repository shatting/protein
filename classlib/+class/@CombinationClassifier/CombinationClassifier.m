classdef CombinationClassifier < class.IClassifier
    % +class@CombinationClassifier
    %   combines two IClassifier instances (subclasses) by building the 
    %   product of their classes
    %
    % PROPERTIES
    %   .class1 (+class@IClassifier)
    %   .class2 (+class@IClassifier)
    %
    % CONSTRUCTORS
    %   obj = CombinationClassifier(class1,class2)
    %
    % ABSTRACT METHODS
    %   cl = getclasses(+data@Dataset)
    %   train()
    %   names = requiredfeaturenames();
    %   ncl = getncl()
    %
    % METHODS
    %   cl = get_classes(+data@Dataset)
    %       checks required features and calls getclasses()
    %   b = ishasreqfeatures(+data@Dataset)
    %
    
    properties
        class1;
        class2;
    end
    
    methods
        function ncl = getncl(obj)
            ncl = obj.class1.getncl*obj.class2.getncl;
        end
        
        function names = requiredfeaturenames(obj)
            % get all combinations
            f1 = obj.class1.requiredfeaturenames;
            f2 = obj.class2.requiredfeaturenames;
            names = {};
            for i=1:length(f1)
                for j=1:length(f2)
                    names = [names {union(f1{i},f2{j})}];
                end
            end
            
        end

        function obj = CombinationClassifier(class1,class2)
            obj.class1 = class1;
            obj.class2 = class2;
        end        
        
        function [combicl] = getcombicl(obj,cl1, cl2)
            combicl = obj.staticgetcombicl(cl1,cl2,obj.class1.getncl, obj.class2.getncl);
        end
    end
    
    methods (Access=protected)
        function cl = getclasses(obj,dataset)
            cl1 = obj.class1.get_classes(dataset);
            cl2 = obj.class2.get_classes(dataset);
            
            cl = obj.getcombicl(cl1, cl2);
            
        end      
        
    end
    
    methods (Static)
        function [cl1, cl2] = getcl1andcl2(combicl, ncl1, ncl2)
            cl = sysXto10(combicl,[ncl1 ncl2]);
            cl1 = cl(:,1);
            cl2 = cl(:,2);
        end
        
        function [combicl] = staticgetcombicl(cl1, cl2, ncl1, ncl2)
            combicl = sysXto10([cl1,cl2],[ncl1 ncl2])';
        end
    end
    
end


classdef IClassifier
    % +class@IClassifier
    %   abstract superclass for classifiers
    %   known subclasses: 
    %       +class@BoundClassifier
    %       +class@PotentialClassifier
    %       +class@CombinationClassifier
    %
    % PROPERTIES
    %   .feature_fn
    %       function handle (default @(x)x)
    %
    % CONSTRUCTORS
    %   obj = IClassifier()
    %       uses default feature_fn @(x)x
    %   obj = IClassifier(featfun)
    %       subclasses can call this at the end of their constructors    
    %
    % ABSTRACT METHODS
    %   cl = getclasses(+data@Dataset)
    %   %train()
    %   names = requiredfeaturenames();
    %   ncl = getncl()
    %
    % METHODS
    %   cl = get_classes(+data@Dataset)
    %       checks required features and calls getclasses()
    %   b = ishasreqfeatures(+data@Dataset)
    %
    
    properties
        feature_fn
    end
    
    methods
        
        function this = IClassifier(featfun)
            if (nargin==0), featfun = @(x)x; end
            this.feature_fn = featfun;
        end
        
        % checks if required features are present in dataset
        function cl = get_classes(this,dataset)
            if (~this.ishasreqfeatures(dataset))
                dprintf('not all features required for classification are present in dataset:\ndataset has {%s}.',stringcell2string(dataset.getfeaturenames));
                for i=1:length(this.requiredfeaturenames)
                    dprintf('classifier set %i: {%s}',i,stringcell2string(this.requiredfeaturenames{i}));
                end
                error('classification not possible.');
            end
            
            cl = this.getclasses(dataset);
        end
    end
    
    methods (Abstract)
        % names is a cell array of cell arrays containing required feature
        % sets         
        %train(obj);
        
        names = requiredfeaturenames(obj);
        
        ncl = getncl(obj);
    end

    methods (Access=protected, Abstract)
        cl = getclasses(obj,dataset)
    end
    
    methods
        
        function b = ishasreqfeatures(obj,dataset)
            names = obj.requiredfeaturenames;
            % comatibility
            if (~iscell(names{1}))
                names = {names};
            end
            
            b = 0;            
            for i=1:length(names),
                b = b || dataset.haveallfeatures(names{i});                    
            end
        end
        
    end
    
end


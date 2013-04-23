classdef HEFGTrainingRun < suffclass.TrainingRun
    % protein.HEFG.HEFGTrainingRun < suffclass.TrainingRun
    %
    %   Encapsulation of a single HEDG Gamma-class clssification run.
    %   Provides persistence by saving run data and via the .load(id) 
    %   method.
    %   
    % PROPERTIES
    %   .options
    %       .b_diagbds      (0) use Diagonal Bounds HEFG Classifier?
    %       .usediary       (0) save all console output?
    %       .usecomment     (0) ask for pre- and post-run comments?
    %       .docluster      (0) perform hierarchical clustering?
    %       .featfun        (@protein.HEFG.feat_pairs) data feature
    %                       function handle
    %       .initial_class_options 
    %                (protein.HEFG.HEFGTrainingRun.default_classoptions())
    %                       options to use for suffclass.classify()
    %                       subroutine    
    %   .info               
    %         .comments     pre-run comments
    %         .timedate     date of training
    %         .pr_comments  post-run comments
    %         .diary        console output
    %
    %   .class_result       return value of suffclass.classify()
    %
    %   .cluster_result     return value of suffclass.clusterh()
    % 
    % PUBLIC METHODS
    %   .run()              start training run, using current .options
    %
    % PUBLIC STATIC METHODS
    %   tr = .load(id)      load training run by id
    %   .assess()           perform assessment of the training run
    %   tr = .new()         factory method, uses default options
    
    
    properties        
        info
        class_result
        cluster_result
    end
    
    methods (Access=private)

    end
    
    methods
               
        function s = get_s_classifier(this)
            s =  protein.HEFG.HEFGSClassifier(this.options.b_diagbds);
        end
        
        function g = get_g_classifier(this)
            g = protein.HEFG.HEFGGammaClassifier(this);            
        end
        
        function pot = get_potential(this)
            pot = this.class_results.pot;
        end
        
        function this = run(this)

            % preliminaries
            clc;
            close all;
            diary off;
            if (exist('HEFG.dry','file')), delete('HEFG.dry'); end;

            if (this.options.usecomment)
                comments = sprintf(input('pre-run comments: ','s'));
            else
                comments = '';
            end
            if (this.options.usediary)
                diary HEFG.dry
            end

            sclassifier = this.get_s_classifier();
            hefg = sclassifier.get_classes(data.FeatureDB.db);

            seq = data.FeatureDB.db.getdata({'aa1','aa2','aa3','aa4'});
            fseq = this.options.featfun(seq);
            %fseq = fseq(1:1000,:);
            
            this.info = struct;
            this.info.comments = comments;
            this.info.timedate = datestr(now,'dd.mm.yy, HH:MM:SS');
            
            dprintf('################ hefg classes #####################');
            this.class_result = suffclass.classify(max(fseq,[],1), fseq,...
                hefg, this.options.initial_class_options);

            if (this.options.docluster)
                cluster_options=struct;    
                cluster_options.clust_minclfact = 0.2;
                cluster_options.clust_maxclfact = 0.8;
                cluster_options.clust_relcostbds = [0.10, 0.40];
                cluster_options.clust_maxdiff = 50;
                cluster_options.clust_minprev = 1e-3;
                cluster_options.clust_minstep = 0.02;    
                %hefg_class_info.options = options.class_options;
                this.cluster_result =...
                    class_info_cluster(this.class_result,cluster_options);   
            end
            
            % HEFG_assessment
            if (this.options.usecomment)
                pr_comments = sprintf(input('post-run comments: ','s'));
            else
                pr_comments = '';
            end
            this.info.pr_comments = pr_comments;

            % save diary
            if (this.options.usediary)
                diary off;    
                fid = fopen('HEFG.dry','r');
                r = [];
                while 1,
                    l = fgetl(fid);
                    if ~ischar(l), 
                        break;
                    end
                    s = sprintf('%s\n',l);
                    r = [r s];
                end
                fclose(fid);
                rundiary = r;
                delete('HEFG.dry');
            else
                rundiary = '';
            end
            this.info.diary = rundiary;

            % save new_deal_data
            k = 1;
            global datadir;
            fname = [datadir,filesep,'HEFG',filesep,'HEFG_gamma_',...
                datestr(now,'yymmdd'),'.mat'];    
            while exist(fname,'file')
                fname = [datadir,filesep,'HEFG',filesep,'HEFG_gamma_',...
                    datestr(now,'yymmdd'),'_',num2str(k),'.mat'];    
                k = k+1;    
            end
            save(fname,'this');

            % display results
            %new_deal_info('new_deal_data')
            this
            dprintf('saved to %s',fname);

            % trainingrun = protein.TrainingRun();
            % trainingrun.from_ndd(new_deal_data);
        end
        
        % HEFG_assesment
        function assess(this) 
            close all;

            %hefg = double(this.new_deal_data.hefg_info.cl(:,1));
            ts = data.FeatureDB.db.getdata({'t','tp'});
            cs = data.FeatureDB.db.getdata({'c','cp'});
            seq = data.FeatureDB.db.getdata({'aa1','aa2','aa3','aa4'});
            sclassifier = this.get_s_classifier();
            gclassifier = this.get_g_classifier();
            hefg = sclassifier.get_classes(data.FeatureDB.db);
            cl_end = gclassifier.get_classes(data.FeatureDB.db);

            if (~isfield(this,'cluster_result') ||...
            ask('clustering info found. assess non-clustered classes?')),
                clso = suffclass.utils.sortcl(cl_end);
                protein.classification_assessment([hefg clso], cs, ts, seq, 10)
                protein.new_deal_quality(hefg, clso, 1,90);
            end

%             if (isfield(new_deal_data.hefg_info,'cluster_result') && ask('clustering info found. assess cluster result classes?')),
%                 class_info_cluster_plot(new_deal_data.hefg_info);
% 
%                 cl = double(new_deal_data.hefg_info.cluster_result.cl);
%                 clso = sortcl(cl);
%                 classification_assessment([hefg clso], cs, ts, seq, 10);
%                 new_deal_quality(hefg, cl, 1,90);
%             end

        end
        
    end
    
    methods (Static)
        
        function tr = load(id)
            
            global datadir;
            load([datadir filesep 'HEFG' filesep 'HEFG_gamma_' id '.mat']);            
            tr = this;
        end
        
        function tr = new()
            tr = protein.HEFG.HEFGTrainingRun();
            tr.options = protein.HEFG.HEFGTrainingRun.default_options();            
        end
        
        function options = default_options(options)

            if (nargin==0 || isempty(options)), options = struct; end
            
            if ~isfield(options,'b_diagbds'),       options.b_diagbds=0; end;
            if ~isfield(options,'usediary'),        options.usediary=0; end;
            if ~isfield(options,'usecomment'),      options.usecomment=0; end;
            if ~isfield(options,'docluster'),       options.docluster=0; end; 
            if ~isfield(options,'featfun'),         options.featfun=@protein.HEFG.feat_pairs; end; %featfun = @(x) x;
            if ~isfield(options,'initial_class_options'),         
                options.initial_class_options = protein.HEFG.HEFGTrainingRun.default_classoptions();
            end            
        end
        
        function class_options = default_classoptions()
                 
            %if (nargin==0 || isempty(options)), options = struct; end

            class_options = struct(   'plots',0,...
                            'max_iterations',inf,...
                            'ask_timeout',0,...
                            'term_percentage',0.02,...
                            'new_classes',0.015,...
                            'remove_small',0.01,...                    
                            'stop_on_nonewclasses',0,...
                            'final_adjust',0,...
                            'savedata',0,...
                            'handleequalpotentials',1,...
                            'supervised',1,...
                            'autropy',1);
            end
        
        
    end
    
end


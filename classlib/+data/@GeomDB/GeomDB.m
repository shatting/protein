classdef GeomDB < handle
    % data.GEOMDB < handle
    %   Singleton class. Provides access to geom.Chain instances.
    % 
    % PROPERTIES
    %   rtfile      source data filename within /data directory. read only
    %   odata       {1xnch} original data from rtfile
    %   chains      {1xnch} geom.Chain objects
    %   afreq       
    %   naa;
    %   nfrag;
    %   res24;
    %   svdbb;
    %   geomversion; 
    %
    % STATIC METHODS
    %   outdb = data.GeomDB.db()
    %   outch = data.GeomDB.ch(i)
    %   outchs= data.GeomDB.chs()
    %   outnch= data.GeomDB.nch()
    %   data.GeomDB.reloadch(i)        reload i-th chain from disk
    %
    % METHODS
    %   [c, idx] = data.GeomDB.db.getsizencell(naamin,naamax) 
    %               returns all chains with length in [naamin,naamax] in a
    %               struct array c. idx corresponds to the index into
    %               data.GeomDB.ch(i)
    %   [ch, i] = data.GeomDB.db.getbyname(search)
    %               exact name search
    
    %% Singleton behaviour
    properties (GetAccess=private, Constant)
        dbprop = data.GeomDB('RT071127.mat');
        %dbprop = data.GeomDB('DSbond.mat');
    end
        
    methods (Static=true)
        function outdb = db
            outdb = data.GeomDB.dbprop;
        end
        
        function och = ch(i)
            och = data.GeomDB.dbprop.chains{i};                        
        end
        function ochains = chs
            ochains = data.GeomDB.dbprop.chains;
        end
        function onch = nch
            onch = length(data.GeomDB.dbprop.chains);
        end
        % reload a chain from original raw data file. useful for debugging.
        function reloadch(chidx)
            % we are assuming that all the other GeomDB fields did not change
            geomdb = data.GeomDB.dbprop; % get reference
            newgeom = geom.Chain(data.GeomDB.db.odata{chidx}); % create new Chain from raw data
            geomdb.chains{chidx} = newgeom; % overwrite old chain reference in .chains
        end
    end

    %% Data properties
    properties (SetAccess=private, GetAccess=private)
        
    end
    properties (SetAccess=private)
        rtfile;
        odata;
        chains;                
        afreq;        
        naa;
        nfrag;
        res24;
        svdbb;
        geomversion;        
    end
    
    methods 
        % bond = data.GeomDB.db{i}.bond;
%         function B = subsref(A,S)
%             B = A.chains{S.subs{1}};
%         end      
        function [c, idx] = getsizencell(obj,naamin,naamax)
            c = {};
            idx = [];
            for i=1:length(obj.chains),
                naai = length(obj.chains{i}.seq);
                if (naai <= naamax && naai >=naamin)
                    c = [c obj.chains{i}];
                    idx = [idx i];
                end
            end
        end
        
        function [ch, i] = getbyname(obj, search)
            ch = [];
            for i=1:length(obj.chains),
                namei = obj.chains{i}.name;
                if (strcmp(search,namei))
                    ch = obj.chains{i};
                    return;
                end
            end
        end
                
    end
    
    %% Methods
    methods (Access=private) % singleton: private constructor
        function obj = GeomDB(rtfile)  
            global datadir;
            load([datadir,filesep,rtfile]);    
            obj.odata = data;
            obj.rtfile = rtfile;
            
            obj = obj.fromdata();            
        end
        % [obj] = fromdata( obj )
        % augment data with geometry data
        function [obj] = fromdata( obj )
            tic
            data = obj.odata;
            %nch=length(data); % number of chains
            fieldnames(data{1})
            nch=length(data)

            % simple statistics
            naa=zeros(1,nch);
            res=zeros(1,nch);
            svdbb=zeros(3,nch);
            afreq=zeros(24,1);
            res24=inf;
            chains = cell(size(data));

            for ch=1:nch,
              if rem(ch,100)==0,
                s = showtime;
                dprintf('computing geometry for chain %i (elapsed %s)',ch,s);
              end;
              dat=data{ch};

              % amino acid frequencies
              aa=dat.seq;% amino acid sequence
              naa(ch)=length(aa);
              res(ch)=dat.res; 
              % if res(ch)==0, res(ch)=inf;data{ch}.res=inf; end;
              for t=1:naa(ch)-3,
                afreq(aa(t))=afreq(aa(t))+1;    
              end;
              if max(aa)>23, res24=min(res24,res(ch)); end;

              chains{ch} = geom.Chain(dat);

              svdbb(:,ch) = chains{ch}.geometry.svdbb;

            end;

            % readme=[
            % '% generated with data2geomdb.m                                 '
            % '% .geomdata{ch} contains (some of) the fields                  '
            % '%   .name    % name of subchain                                '
            % '%   .res     % resolution (-1 = NMR)                           '
            % '%   .seq     % amino acid sequence, uint8                      '
            % '%   .bond    % bond vectors, 3 columns, int8                   '
            % '%   and all fields from bond2geomstruct.m                      '
            % ];

            obj.chains = chains;
            %obj.nch = nch;
            obj.afreq = afreq;
            obj.naa = sum(naa);
            obj.nfrag = sum(naa-3);
            obj.res24 = res24;
            obj.svdbb = svdbb;
            obj.geomversion = geom.Geometry.version;

            end

        
    end
        
end


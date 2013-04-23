classdef suffstat < handle
% +suffclass@suffstat < handle
%
%   Main class of suffclass that provides methods for basic manipulation of 
%   sufficient statistics: adding or subtracting data, modifying the
%   number or merging of classes.
%   A sufficient statistic summarizes information about data vectors. 
%   Vectors may consist of any number of numerical (real numbers), 
%   categorical (integers with meaningless mean value) or ordinal (integers 
%   with meaningful mean value) features (nyi). suffstat is a handle class 
%   and is initialized with a call to the constructor
%
%       suff = suffstat(ncl,typ)
%   
%   where ncl, the number of classes, is an integer, and typ is a
%   1x(number of features) vector describing the type of each feature:
%   
%       typ(i)>0: categorical variable (meaningless mean) with typ(i) levels
%       typ(i)=0: numerical variable
%   
%   The following call adds information from the vectors contained in data to
%   the object suff, weighing each vector according to its weight given
%   by wt:
%
%       suff.adddata(data,cl,wt)
%
%   cl is a vector containing class labels for data. Please note, 
%   that any vectors with lables not falling into the range 1:suff.ncl will 
%   be ignored in the call and will not be added to the statistic. 
%
% PUBLIC PROPERTIES
%   .ncl       number of classes
%   .typ       feature types - see constructor
%   .nnum      number of numerical features
%   .ncat      number of numerical features
%   .nord      number of ordinal features (nyi)
%   .num       numerical features indexes
%   .cat       categorical features indexes
%   .ord       ordinal features indexes (nyi)
%   .moment    numerical moment matrix
%   .count     categorical moment matrix
%
% CONSTRUCTORS
%   obj = suffstat(ncl,type)
%       Initialize sufficient statistic for data of arbitrary type 
%       (>=1 feature) and ncl classes.
%       INPUT:
%           ncl         number of classes
%           typ(i)      type of feature i
%               typ(i)>0: categorical variable (meaningful mean) with 
%                         typ(i) levels
%               typ(i)=0: numerical variable
%       OUTPUT:
%           suff        SUFFSTAT object handle.
%
% PUBLIC METHODS    
%   obj = .clone()
%   .increasencl(ncl)            
%       Increase the number of classes to ncl.
%   .mergeclasses(mergelist) TODO: seems buggy
%       Merge classes from mergelist into one class.
%       INPUT:
%           mergelist   List of classes to be merged. New class will have
%                       label equal to minimal label of merged classes. All
%                       other class lables from mergelist are empty after
%                       call.                
%  .adddata(data,cl,wt);
%       Increment sufficient statistic using data with assigned 
%       classes cl and weights wt.
%  .subtractdata(data,cl,wt);
%       Decrement sufficient statistic using data with assigned 
%       classes cl and weights wt.
%  .reassigndata(data,clold,wtold,clnew,wtnew)
%       modify sufficient statistic using data when class assignments 
%       change from cold to clnew with new weights wtnew and old 
%       weights wtold.
%  .addsuff(insuff, intonewclasses)               
%       Add sufficient statistic suff to object.
%       INPUT:
%           intonewclasses     true:    create new class for each class in 
%                                       insuff
%                              false:   merge classes
    properties
        ncl
        typ
        nnum
        ncat
        nord
        num
        cat
        ord
        moment
        count
    end
    
    methods        
        
        function obj = suffstat(ncl,type)
          
          obj.ncl = ncl;
          obj.typ = type;
          
          obj.num=find(type==0);
          obj.nnum=length(obj.num);
          obj.cat=find(type~=0);
          obj.ncat=length(obj.cat);
          typc=type(obj.cat); % restrict type to categorical data
          obj.ord=find(typc<0);
          obj.nord=length(obj.ord);
          obj.moment=zeros(obj.nnum+1,obj.nnum+1,ncl);
          obj.count=zeros(obj.ncat,max(typc),ncl);
                    
        end
                        
        function obj = clone(inobj)
            obj = suffstat(inobj.ncl,inobj.type);
            obj.moment = inobj.moment;
            obj.count = inobj.count;
        end
             
        function increasencl(obj, ncl)            
            if ncl<=obj.ncl, 
                warning('there are already enough classes');
                return;
            end;
            obj.ncl=ncl;
            obj.moment(1,1,ncl)=0;
            obj.count(1,1,ncl)=0;
        end
        
        function mergeclasses(obj, mergelist)
            
            newcl = min(mergelist);
            
            obj.moment(:,:,newcl) = sum(obj.moment(:,:,mergelist),3);     
            % want the classes to stay in place, eg. for clustering
            obj.moment(:,:,mergelist(2:end)) = obj.moment(:,:,mergelist(2:end))*0;
            
            obj.count(:,:,newcl) = sum(s.count(:,:,mergelist),3);
            obj.count(:,:,mergelist(2:end)) = obj.count(:,:,mergelist(2:end))*0;  
            
        end

        function adddata(obj, data, cl, wt)
            
            m = min(data(:,obj.cat),[],1);
            if (any(m<=0)),
                error('categorical features must be > 0.');
            end
            
            if (nargin >= 4), withweights=1; else withweights=0;end
            % TODO: parfor
            for cli=1:obj.ncl,
                %if prt, disp(['increment class ',num2str(cli)]); end;
                clind=find(cl==cli);
                nclind=length(clind);
                % categorical part
                for i=1:obj.ncat,
                    %if prt, disp(['i=',num2str(i)]); end;
                    ii=obj.cat(i);
                    for a=1:obj.typ(ii),
                        aaidx = data(clind,ii)==a;
                        if withweights,
                            obj.count(i,a,cli)=obj.count(i,a,cli)+sum(wt(clind(aaidx)));
                        else
                            obj.count(i,a,cli)=obj.count(i,a,cli)+sum(aaidx);
                        end
                    end;
                end;
                % numerical part
                if obj.nnum>0,
                    if withweights,
                        X=[data(clind,obj.num), ones(nclind,1)];  % [x_i 1]
                        Y=X.*repmat(wt(clind),1,obj.nnum+1);      % [x_i*w_i w_i]
                        X=X'*Y;
                    else
                        X=[data(clind,obj.num), ones(nclind,1)];
                        X=X'*X;
                    end
                    obj.moment(:,:,cli)=obj.moment(:,:,cli)+X;        
                else
                    if withweights,
                        obj.moment(:,:,cli)=obj.moment(:,:,cli)+sum(wt(clind));
                    else
                        obj.moment(:,:,cli)=obj.moment(:,:,cli)+nclind;
                    end
                end % nnum > 0
            end % for
        end % adddata
              
        function subtractdata(obj,data,cl,wt)
            
            if (nargin >= 4), withweights=1; else withweights=0;end
            
            for cli=1:obj.ncl,
                %if prt, disp(['decrement class ',num2str(cli)]); end;
                clind=(cl==cli);
                nclind=length(clind);
                for i=1:obj.ncat,
                    ii=obj.cat(i);
                    for a=1:obj.typ(ii),
                        aaidx = data(clind,ii)==a;
                        if withweights,
                            obj.count(i,a,cli)=obj.count(i,a,cli)-sum(wt(clind(aaidx)));
                        else
                            obj.count(i,a,cli)=obj.count(i,a,cli)-sum(aaidx);
                        end
                    end;
                end;
                if obj.nnum>0,
                    if withweights,
                        X=[data(clind,obj.num), ones(nclind,1)];  % [x_i 1]
                        Y=X.*repmat(wt(clind),1,obj.nnum+1);      % [x_i*w_i w_i]
                        X=X'*Y;
                    else
                        X=[data(clind,obj.num), ones(nclind,1)];
                        X=X'*X;
                    end
                    obj.moment(:,:,cli)=obj.moment(:,:,cli)-X;
                else
                    if withweights,
                        obj.moment(:,:,cli)=obj.moment(:,:,cli)-sum(wt(clind));
                    else
                        obj.moment(:,:,cli)=obj.moment(:,:,cli)-nclind;
                    end
                end;
            end;
        end % subtractdata
        
        function reassigndata(obj,data,clold,wtold,clnew,wtnew)
            obj.adddata(data,clnew,wtnew);
            obj.subtractdata(data,clold,wtold);
        end
        
        function addsuff(obj, suff, intonewclasses)            
            if nargin < 3 || ~intonewclasses,                
                obj.moment=obj.moment+suff.moment;
                obj.count=obj.count+suff.count;
            else
                obj.ncl=obj.ncl+suff.ncl;                
                obj.moment(:,obj.ncl+1:suff.ncl)=suff.moment;
                obj.count(:,obj.ncl+1:suff.ncl)=suff.count;
            end
        end
        
    end
        
end


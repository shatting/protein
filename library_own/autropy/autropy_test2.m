
if ~exist('autropy_test2_loaded','var'),
    if ~exist('autropy_test2.mat','file'),
        getgeomseq;

        % options
        featfun = @feat_pairs;
        %featfun = @(x) x;

        data = featfun(seq);
        %fseq = fseq(1:1000,:);
        %HEL27_ccpt = HEL27_ccpt(1:1000,:);

        cl = HEL27;

        suff=suffstat(max(cl), max(data,[],1));
        suff=suffstat(suff,data,size(data,1),cl);
        save autropy_test2 data cl suff;
    else
        load autropy_test2;
    end
    autropy_test2_loaded = 1;
end

autropy_test_common
% formerly geomseq.m
% HEL27
% uses: RT071127.mat
%       db2fragdb.m
%       fragdb2linear.m
%       HEL27_anal08.m

if 0,
  if (exist('db')~=1),
    cd ..;startup;cd allesneu
    db = data2db( 'RT071127.mat' );    
  end
  if ~exist('fragdb'),
      fragdb = db2fragdb(db);
  end
elseif ~exist('fragdb'),
  tic
  disp('load data')
  load fragdb
  showtime
  %  nfrag: 221491
  %  freq(20,20,20,20)       frequency of primary sequence
  %  frags{20,20,20,20} 
  %      count               frequency
  %      data(count,nfeat)   features   
  %  nfeat: 6
  %  feats: {'cos(i)'  'cos(i+1)'  'tors(i,3)'  'len(i)'  'len(i+1)'  'res'}
  disp(fragdb)
end;

if ~exist('geom','var') || ~exist('seq','var'),
    [ geom,seq ] = fragdb2linear( fragdb );
end

HEL27 = HEL27_anal08(geom);

HEL27_ccpt = sys10toX(HEL27,[3 3 3]); % c_i, c_i+1 and t classes
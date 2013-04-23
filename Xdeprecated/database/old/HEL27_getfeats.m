% convenience script
loadgeomdb;

if ~exist('HEL27feats','var'),    
    [HEL27feats, fragaaseq] = geomdb2features(geomdb,{'c','cp','beta'});    
    hel27 = HEL27_classifier(HEL27feats);
    dprintf('loaded HEL27feats, hel27 and fragaaseq');
else
    dprintf('nothing loaded, clear HEL27feats to reload.');
end

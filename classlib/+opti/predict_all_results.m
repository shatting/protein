minrms = min(rrmsds,[],2);
idx= minrms > 0;
plot(len(idx),minrms(idx),'r+')
xlabel('Chain length');
ylabel('Minimum RMSD');
box off
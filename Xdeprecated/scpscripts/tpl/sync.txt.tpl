open #{username}#@login.mat.univie.ac.at
cd #{rprotdir}#
option exclude ".svn"
synchronize remote
chmod 755 chmodrec.sh
chmod 755 optimize.sh
call ./chmodrec.sh
close
exit
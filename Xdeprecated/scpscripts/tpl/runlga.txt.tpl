open virgy@login.mat.univie.ac.at
cd protein/stephan/lga/lga_bin/MOL2
lcd stephan
lcd lga
lcd lga_bin
lcd MOL2
put #{infile}#
cd ..
lcd ..
put runlga.sh
chmod 755 /users/virgy/protein/stephan/lga/lga_bin/runlga.sh
call /users/virgy/protein/stephan/lga/lga_bin/runlga.sh #{infile}#

cd TMP
get firsttry.pdb
cd .. 

close
exit
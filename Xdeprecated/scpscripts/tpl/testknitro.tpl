#!/bin/bash
CDIR=#{rampldir}#
AMPL=~herman/COCONUT/bin/ampl
KNITRODIR=/local/scratch/dferi/knitro/knitro-5.1.2-z/

# autoclosing ssh tunnel
# http://www.g-loaded.eu/2006/11/24/auto-closing-ssh-tunnels/
ssh -f -L 8349:gondor.mat.univie.ac.at:8349 -l #{username}# login.mat.univie.ac.at "sleep 0.5"
export ZIENA_LICENSE_NETWORK_ADDR=127.0.0.1

cd "$KNITRODIR"knitroampl/
"$AMPL" "$CDIR"/test1.run
cd "$CDIR"

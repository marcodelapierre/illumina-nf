#!/bin/bash

# python 3
export PATH=/software/pawsey0001/spack-tests/setonix/2022.01/software/cray-sles15-zen2/gcc-10.3.0/python-3.9.7-2wjic6fswcgdxvcpo3cwcnn5hti62t3n/bin:$PATH

# singularity
export PATH="/software/pawsey0001/spack-tests/setonix/2022.01/software/cray-sles15-zen2/gcc-10.3.0/singularity-3.8.5-eig5ikxq2yyrkklutfqpwq6khk2mcxah/bin:$PATH"
# maybe due to fs setup for home?
export SINGULARITY_HOME="/home/mdelapierre"
export SINGULARITY_BINDPATH="/home,/lus/joey/home"

# shpc
export PATH=/home/mdelapierre/containers-shpc-nf-tests/shpc/shpc/bin:$PATH
export PYTHONPATH=/home/mdelapierre/containers-shpc-nf-tests/shpc/shpc/lib/python3.9/site-packages:$PYTHONPATH

# shpc modules
module use /home/mdelapierre/containers-shpc-nf-tests/shpc/containers/modules

# spack
. /home/mdelapierre/Spack-sprints/setonix/spack/share/spack/setup-env.sh

spack load openjdk
spack load nextflow

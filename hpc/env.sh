# Shared environment for every job in hpc/. Sourced, not run.
T=/faststorage/project/DeiC-KU-L04/dduchene/tools
export RLIB_COLOUR=$T/Rlib_colour
export CMDSTAN_DIR=$T/cmdstan_colour
export CMDSTAN=$(ls -d $CMDSTAN_DIR/cmdstan-* 2>/dev/null | sort -V | tail -1)
IMG=/faststorage/project/DeiC-KU-L04/shared_software/containers/r-ver_4.4.1.sif
# Colour library first, then the existing library (brms, rstan, ape, MCMCglmm ...)
RSCRIPT="apptainer exec --env R_LIBS_USER=$RLIB_COLOUR:$T/Rlib,RLIB_COLOUR=$RLIB_COLOUR,CMDSTAN_DIR=$CMDSTAN_DIR,CMDSTAN=$CMDSTAN,OMP_NUM_THREADS=1 $IMG Rscript"

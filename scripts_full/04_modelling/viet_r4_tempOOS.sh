#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=48:00:0
#$ -pe smp 16
#$ -l mem=8G
#$ -l tmpfs=15G
#$ -N v_tempOOS
#$ -wd /home/ucbtgib/Scratch/dengue_viet_23/myriad_output

cd $TMPDIR

module unload compilers
module unload mpi
module load udunits/2.2.20/gnu-4.9.2
module load perl/5.22.0
module load python/2.7.12
module load proj.4/6.1.0
module load gdal/3.0.4/gnu-4.9.2
module load geos/3.5.0/gnu-4.9.2

/home/ucbtgib/Scratch/R-4.1.1/bin/Rscript /home/ucbtgib/Scratch/dengue_viet_23/03a_r4_vietnam_tempOOS.R